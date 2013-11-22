import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmap
plt.interactive(True)
from scipy import fftpack
import scipy.ndimage
import scipy
from spefiles import MyPrincetonSPEFile


def PSF( X,Y, xpos, ypos, sigma ):
    psf = 1/(2*np.pi*sigma**2) * np.exp( - .5 *( ((X-xpos)/sigma)**2 + ((Y-ypos)/sigma)**2 ))
    return psf

def sm_simulate_image( Nobjects=100, \
                           ellipsoid_coeffs=[1,1,1], \
                           Efield=[1,0,0], \
                           pixel_edge_length=350e-9, \
                           wavelength=6e-7):

    Npixelx = 512
    Npixely = 512

    # ccd grid
    X,Y   = np.meshgrid( np.arange(Npixelx), np.arange(Npixely) )
    image = np.zeros_like(X, dtype=np.float)

    # width of the psf is approximated as the wavelength (and
    # given in units of pixel edge length)
    psf_sigma = wavelength/pixel_edge_length

    # laser spot params
    ls_x = 250.0
    ls_y = 250.0
    ls_width   = 100.0

    for i in range(Nobjects):
        # pick random position
        xpos = np.random.random()*np.max(X)
        ypos = np.random.random()*np.max(Y)

        # get psf on fine grid
        psf = PSF(X,Y, xpos, ypos, psf_sigma)
        
        # now decide quantum yield (random for now)
        QY = np.random.random()

        # and query laser intensity in wide field
        I  = 1/(2*np.pi*ls_width**2) * np.exp( -.5 *( ((xpos-ls_x)/ls_width)**2 + ((ypos-ls_y)/ls_width)**2 ))

        # add that to the image
        image += I*QY*psf

    print np.max(image)

    # add laser spot background luminescence
    lsbg = 1e-1/(2*np.pi*ls_width**2) * np.exp( -.5 *( ((X-ls_x)/ls_width)**2 + ((Y-ls_y)/ls_width)**2 ))
    image += lsbg

    # scale up to get some counts
    image *= 1e9
    image = np.around( image )

    # add poisson noise
    image += np.random.poisson( lam=100, size=image.shape )

    # add broad some weird wobble to test frequency chopping closer to sm 
    image += np.around( 50*np.sin(2*np.pi*X/25.0) )
    image += np.around( 50*np.cos(2*np.pi*Y/25.0) )

    return X,Y,image

def gaussian_filter( image, sigma=2 ):
    return scipy.ndimage.filters.gaussian_filter( image, sigma=sigma )
    
def white_tophat_filter( image, size=2 ):
    return scipy.ndimage.white_tophat( image, size=(size,size) )

def fourier_filter(image, cutfreq1=12.0, cutfreq2=2.0, operation='normal'):
    fimage  = fftpack.fft2( image )
    fsimage = fftpack.fftshift( fimage )
    powerspec = np.abs(fsimage)**2
    
    if operation=='normal' or operation=='gaussian':
        mask = np.ones_like(powerspec)
    elif operation=='inverted':
        mask = np.zeros_like(powerspec)
    else:
        raise hell

    cutoff1 = np.argmin( np.abs( fftpack.fftfreq(512)-1/cutfreq1 ) )
    cutoff2 = np.argmin( np.abs( fftpack.fftfreq(512)-1/cutfreq2 ) )
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            if operation=='normal':
                if (i-mask.shape[0]/2.0)**2 + (j-mask.shape[1]/2.0)**2 < cutoff1**2:
                    mask[i,j] = 0
            elif operation=='inverted':                    
                if (i-mask.shape[0]/2.0)**2 + (j-mask.shape[1]/2.0)**2 < cutoff1**2:
                    mask[i,j] = 1
            elif operation=='gaussian':
                mask[i,j] -= np.exp( -.5*( \
                        (j-mask.shape[1]/2.0)**2/float(cutoff1)**2 \
                            + (i-mask.shape[0]/2.0)**2/float(cutoff1)**2 ))
            if (i-mask.shape[0]/2.0)**2 + (j-mask.shape[1]/2.0)**2 > cutoff2**2:
                if operation=='normal':
                    mask[i,j] = 0
                elif operation=='inverted':
                    mask[i,j] = 1
                elif operation=='gaussian':
                    mask[i,j] = 0


    mask = fftpack.ifftshift(mask)
    newimage = np.real( fftpack.ifft2( mask*fimage ) )
    return newimage, fftpack.fftshift(mask)


def mark_all_clusters( clusterimage ):
    # the image is supposed to be of dtype=int, where 
    # 0 means that a pixel did not meet the threshold,
    # -1 means that the pixel hasn't been assigned to a cluster yet,
    # 1...n means that the pixel belongs to the cluster with this number

    icluster = 0

    for ix in range(clusterimage.shape[1]):
        for iy in range(clusterimage.shape[0]):
            if clusterimage[iy,ix]==-1:
                icluster += 1
                clusterimage = mark_this_cluster_recursively(clusterimage,iy,ix,icluster)

    return clusterimage, icluster
    
def mark_this_cluster_recursively( clusterimage, iy, ix, icluster ):
    # set cluster marker
    clusterimage[iy,ix] = icluster

    # is the pixel to the left also above threshold?
    if ix!=0:
        if clusterimage[iy, ix-1]==-1:
            # call self to mark it
            clusterimage = mark_this_cluster_recursively( clusterimage, iy, ix-1, icluster )

    # is the pixel above also above threshold?
    if iy!=0:
        if clusterimage[iy-1, ix]==-1:
            # call self to mark it
            clusterimage = mark_this_cluster_recursively( clusterimage, iy-1, ix, icluster )

    # is the pixel to the right also above threshold?
    if ix!=clusterimage.shape[1]-1:
        if clusterimage[iy, ix+1]==-1:
            # call self to mark it
            clusterimage = mark_this_cluster_recursively( clusterimage, iy, ix+1, icluster )

    # is the pixel below also above threshold?            
    if iy!=clusterimage.shape[0]:
        if clusterimage[iy+1, ix]==-1:
            # call self to mark it
            clusterimage = mark_this_cluster_recursively( clusterimage, iy+1, ix, icluster )

    return clusterimage


def find_cluster_positions( clusterimage, min_size_threshold=4, max_size_threshold=20, originalimage=None ):
    # go through the image and add clusters to list (accumulating pixel positions and number) 
    c = []
    for yi in range(clusterimage.shape[0]):
        for xi in range(clusterimage.shape[1]):
            if not clusterimage[yi,xi]==0:
                c.append( [clusterimage[yi,xi], yi, xi] )
    c = np.array(c)
    Nclusters = np.max(c[:,0])

    if originalimage==None:
        originalimage = np.ones_like(clusterimage)

    # go through all clusters and:
    # add up pixel-value-weighted coordinates (rows and columns), 
    # also add up the pixel-values, and, lastly, keep track of the number of pixel we've added up.
    clusterpos = np.zeros( (Nclusters,4) )
    for ci in range(Nclusters):
        for a in c:
            if a[0]==ci+1:
                clusterpos[ci,0] += a[1] * originalimage[a[1],a[2]]   # weighted rows
                clusterpos[ci,1] += a[2] * originalimage[a[1],a[2]]   # weighted columns
                clusterpos[ci,2] += originalimage[a[1],a[2]]          # weights
                clusterpos[ci,3] += 1                                 # number of pixel

    # now get the positions by dividing the accumulated 
    # positions with the accumulated pixel values 
    pos = []
    for cp in clusterpos:
        # accept cluster only if it meets pixel-size criteria:
        if (cp[3] >= min_size_threshold) and (cp[3] <= max_size_threshold):
            pos.append( [cp[0]/cp[2], cp[1]/cp[2]] )
    pos = np.array(pos)

    return pos


if __name__=="__main__":
    X,Y,image = sm_simulate_image()

    image = MyPrincetonSPEFile('yuxis_sm_images/SNM0252.SPE').return_Array()
    image = np.squeeze( image.astype(np.float) )


    plt.figure()
    plt.imshow( np.log10(image), interpolation='nearest')
    plt.colorbar()

    newimage, mask = fourier_filter(image, cutfreq1=30.0, cutfreq2=1.0, operation='inverted')

    plt.figure()
    plt.imshow( newimage, interpolation='nearest')

    raise SystemExit
    clusterimage = newimage > np.std(newimage.flatten())*6
    clusterimage = -clusterimage.astype(np.int)
    plt.figure()
    plt.imshow( clusterimage, interpolation='nearest')


    cimage = clusterimage.copy()
    # clean up boundaries (got Fourier spills there)
    bskip=10
    cimage[:bskip,:]  = 0   # top edge
    cimage[:,:bskip]  = 0   # left edge
    cimage[-bskip:,:] = 0   # bottom edge
    cimage[:,-bskip:] = 0   # right edge

    cimage = mark_all_clusters(cimage)
    plt.figure()
    plt.imshow( cimage, interpolation='nearest')
    plt.colorbar()

