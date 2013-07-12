import os
import time
import numpy as np
import matplotlib.pyplot as plt
plt.interactive(True)
from datetime import datetime, timedelta
from mpl_toolkits.axes_grid1 import AxesGrid
import memory
import h5py

def show_mem():
    print "memory usage:"
    print memory.memory()/(1024*1024)
    print memory.resident()/(1024*1024)


def deal_with_date_time_string( motorobj, datetimestring ):        
    """This function converts the date+time string into the difference time 
    (in seconds) since the start of the experiment. The first value passed 
    is treated in a special way: it sets the zero time mark.

    It is the same function for both motors, that's why it's here,  
    outside the class definitions.
    """
    dt = datetime.strptime(datetimestring, "%m/%d/%Y %H:%M:%S.%f")
    if motorobj.experiment_start_datetime == None:
        motorobj.experiment_start_datetime = dt
        return 0.0
    else:
        td = dt - motorobj.experiment_start_datetime
        return td.total_seconds()
    # return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 1e6 



def grid_image_section_into_squares_and_define_spots( movie, res, bounds ):
    
    rb = bounds #[80,56,155,89]  # pixel indices (starting from zero!)

    leftedges = range(rb[1],rb[3],res)
    if leftedges[-1]+res > rb[3]:
        leftedges = leftedges[:-1]

    topedges = range(rb[0],rb[2],res)
    if topedges[-1]+res > rb[2]:
        topedges = topedges[:-1]

    # Note: Spot definitions will include all edges, ie define_spot( [0,2,10,11] ) will
    # create a spot covering 11 pixel in x (columns) and 10 in y (rows).
    # At the same time, the leftedges/topedges generated here will not include the last
    # edge (rb[3] in columns

    for xi in leftedges:
        for yi in topedges:
            movie.define_spot( [yi,xi,yi-1+res,xi-1+res] )

    return 


def trim_noisy_data( movie, what='M_ex', threshold=None ):

    if what=='M_ex' or what=='M_em':
        if threshold==None:
            threshold = 3*movie.bg_spot.std
        quantity  = 'mean_intensity'

    elif what=='phase_ex':
        if threshold==None:
            threshold = 0.1
        quantity  = 'M_ex'

    elif what=='phase_em':
        if threshold==None:
            threshold = 0.1
        quantity  = 'M_em'

    elif what=='LS':   # LS is a bit special here
        if not hasattr(movie.validspots[0],'phase_ex_trimmed'):
            trim_noisy_data( movie, what='phase_ex', threshold=None )
        if not hasattr(movie.validspots[0],'phase_em_trimmed'):
            trim_noisy_data( movie, what='phase_em', threshold=None ) 
        for si,s in enumerate(movie.validspots):
            setattr(s,'LS_trimmed', getattr(s,'phase_ex')-getattr(s,'phase_em') )
        return  # early return from this one

    # all others go through this one
    for si,s in enumerate(movie.validspots):
        if getattr(s,quantity) <= threshold:
            setattr(s,what+'_trimmed', np.nan )
        else:
            setattr(s,what+'_trimmed', getattr(s,what) )


def update_image_files( movie, what, fileprefix ):
    filename = fileprefix + what + '_output_data.txt'
    # if the file exists, try loading and updating it
    if os.path.isfile(filename):
        try: 
            sc = np.loadtxt(filename)
        except IOError:
            sc = getattr(movie, what+'_image')
        # where the current image is not nan, update the old image
        new_value_indices = ~np.isnan( getattr(movie, what+'_image') ) 
        sc[new_value_indices] = getattr(movie, what+'_image')[new_value_indices]
        np.savetxt( filename, sc )
    else:
        np.savetxt( filename, getattr(movie, what+'_image') )


def save_hdf5( movie, fileprefix, proc, images=True, spots=True ):

    # first we grab all the data we have now

    if images:
        # compile dict of images
        imagedict = {}
        for d in dir(movie):
            if d.endswith('_image') and type(getattr(movie,d))==np.ndarray:
                imagedict[d] = getattr(movie,d)

    # now let's see if there's an existing file
    filename = fileprefix + movie.data_basename + '_output'+'_proc'+str(proc)+'.hdf5'
    print 'Looking for output file with name %s ...' % filename
    if os.path.isfile(filename):
        try: 
            # read the images into a dict
            rfid = h5py.File(filename,'r')
            rfidimages = rfid['images']
            readimagedict = {}
            for item in rfidimages.items():
                readimagedict[item[0]] = np.array( item[1] )
            rfid.close()
            print 'Found existing output file, updating data'

            # update read images with current images
            for item in imagedict.items():
                new_value_indices = ~np.isnan( imagedict[item[0]] )
                readimagedict[item[0]][new_value_indices] = imagedict[item[0]][new_value_indices]
                # now readimagedict has the latest values

        except IOError:
            raise IOError('Something blocking file access?')

    else:
        print 'No output file found; creating it now'
        readimagedict = imagedict

    fid = h5py.File(filename,'w')

    fid.create_group('images')
    for imagename in readimagedict:
        fid.create_dataset('images/'+imagename, data=readimagedict[imagename])
    
    # if spots:
    #     # basically save spot data here: first, we need the coordinates, which
    #     # are stored as a single [Nspots x 4] matrix
    #     Nspots = len(movie.validspots)
    #     coordmatrix = np.zeros( (Nspots,4) )
    #     for si,s in enumerate(movie.validspots):
    #         coordmatrix[si,:] = s.coords
        
    #     # now we need the fit data (I0, M0, phase) for each spot, for each 
    #     # portrait, for each line...  [Nspots x Nportraits x Nlines x 3]
    #     Nportraits = len( movie.validspots[0].portraits )
    #     Nlines     = len( movie.validspots[0].portraits[0].lines )
    #     fitmatrix = np.zeros( (Nspots,Nportraits,Nlines,3) )
    #     for si,s in enumerate(movie.validspots):
    #         for pi,p in enumerate(s.portraits):
    #             for li,l in enumerate(p.lines):
    #                 fitmatrix[si,pi,li,:] = [l.I0, l.M0, l.phase]

    # fid.create_group('spots')
    # fid.create_dataset('spots/coordinates', data=coordmatrix)
    # fid.create_dataset('spots/fits',data=fitmatrix)

    fid.close()


def combine_outputs( basename, fileprefix ):
    # collect list of files    
    filelist = []
    for file in os.listdir("."):
        if file.endswith('.hdf5') and file.startswith(basename+'_output'):
            filelist.append(file)

    # now grab the images from the first file (whichever that one is; doesn't matter)
    rfid = h5py.File(filelist[0],'r')
    rfidimages = rfid['images']
    readimagedict = {}
    for item in rfidimages.items():
        readimagedict[item[0]] = np.array( item[1] )
    rfid.close()

    # now do the same with the remaining files, updating the previous images
    for f in filelist[1:]:
        rfid = h5py.File(f,'r')
        rfidimages = rfid['images']
        imagedict = {}
        for item in rfidimages.items():
            imagedict[item[0]] = np.array( item[1] )
        # update first images with current images
        for item in imagedict.items():
            new_value_indices = ~np.isnan( imagedict[item[0]] )
            readimagedict[item[0]][new_value_indices] = imagedict[item[0]][new_value_indices]
        rfid.close()

    # now save those to a new file
    fid = h5py.File(basename+'_output.hdf5','w')
    fid.create_group('images')
    for imagename in readimagedict:
        fid.create_dataset('images/'+imagename, data=readimagedict[imagename])
    fid.close()


def save_spot_data( movie, what='M_ex', whole_image=True, fileprefix='' ):
    # We assume in the following that the first spot (movie.spots[0]) is
    # the upper left spot in the image, with the origin (0,0) of the
    # image also being in the upper left corner of the plot.
    # xinit = movie.spots[0].coords[0]
    # yinit = movie.spots[0].coords[1]

    xdim=movie.spots[-1].coords[2]-movie.spots[0].coords[0]+1
    ydim=movie.spots[-1].coords[3]-movie.spots[0].coords[1]+1
    xinit = movie.spots[0].coords[0]
    yinit = movie.spots[0].coords[1]
    
    if whole_image:
        fs = np.ones( (movie.camera_data.datasize[1],movie.camera_data.datasize[2]) ) * np.nan
        xinit = 0
        yinit = 0
    else:
        fs = np.ones((ydim,xdim)) * np.nan
    
    for si,s in enumerate(movie.validspots):
        xi = s.coords[0]-xinit
        xf = s.coords[2]-xinit+1  # edges...
        yi = s.coords[1]-yinit
        yf = s.coords[3]-yinit+1
        fs[yi:yf,xi:xf] = getattr(s,what)

    np.savetxt(fileprefix + what + '_output_data.txt', fs)

        

def show_spot_data( movie, what='M_ex', which_cmap=None, show_bg_spot=True ):

    import matplotlib.pyplot as plt
    import matplotlib.cm as cmap
    from matplotlib.patches import Rectangle

    if which_cmap==None:
        colormap = cmap.jet

    plt.figure()
    
    # draw average as background, use colormap gray
#    plt.imshow( np.mean( movie.camera_data.rawdata, axis=0 ), cmap=cmap.gray )
    plt.imshow( movie.camera_data.average_image, cmap=cmap.gray )

    ax = plt.gca()

    if show_bg_spot:
        p = Rectangle( (movie.bg_spot.coords[0],movie.bg_spot.coords[1]), \
                           movie.bg_spot.width, movie.bg_spot.height, \
                           facecolor=[1,0,0,0.1], edgecolor=[1,0,0,.75])
        ax.add_patch( p )

    # prepare intensities, etc
    mean_intensities = []
    Ms_ex     = []
    Ms_em     = []
    phases_ex = []
    phases_em = []
    LSs       = []
    ET_rulers = []
    residuals = []

    # collect the data from all spots
    for ss in movie.validspots:
        mean_intensities.append( ss.mean_intensity )
        Ms_ex.append( ss.M_ex )
        Ms_em.append( ss.M_em )
        phases_ex.append( ss.phase_ex )
        phases_em.append( ss.phase_em )
        LSs.append( ss.LS )
        ET_rulers.append( ss.ET_ruler )
        residuals.append( ss.residual )
    # turn the lists into numpy arrays
    intensity = np.array(mean_intensities)
    M_ex = np.array(Ms_ex)
    M_em = np.array(Ms_em)
    phase_ex = np.array(phases_ex)
    phase_em = np.array(phases_em)
    LS = np.array(LSs)
    ET_ruler = np.array(ET_rulers)
    residual = np.array(residuals)

    # rescale values to color range 
    intensity_color = (intensity-np.min(intensity))/np.max(intensity)                
    M_ex_color = (M_ex-np.min(M_ex))/np.max(M_ex)
    M_em_color = (M_em-np.min(M_em))/np.max(M_em)
    phase_ex_color = (phase_ex-np.min(phase_ex))/np.max(phase_ex)
    phase_em_color = (phase_em-np.min(phase_em))/np.max(phase_em)
    LS_color = (LS-np.min(LS))/np.max(LS)
    ET_ruler_color = (ET_ruler-np.min(ET_ruler))/np.max(ET_ruler)
    residual_color = (residual-np.min(residual))/np.max(residual)

    # edges! the +1 accounts for the fact that a spot includes its edges (really? does it?)
    xdim=movie.spots[-1].coords[2]-movie.spots[0].coords[0]+1
    ydim=movie.spots[-1].coords[3]-movie.spots[0].coords[1]+1
    fs=np.zeros((ydim,xdim))
    xinit = movie.spots[0].coords[0]
    yinit = movie.spots[0].coords[1]
     
    for si,s in enumerate(movie.validspots):
        xi = s.coords[0]-xinit
        xf = s.coords[2]-xinit+1  # edges...
        yi = s.coords[1]-yinit
        yf = s.coords[3]-yinit+1

        # determine color (color axes)
        if what=='M_ex':
            col = colormap( M_ex_color[si] )
            fs[yi:yf,xi:xf] = s.M_ex
        elif what=='M_em':
            col = colormap( M_em_color[si] )
            fs[yi:yf,xi:xf] = s.M_em
        elif what=='phase_ex':
            col = colormap( phase_ex_color[si] )
            fs[yi:yf,xi:xf] = s.phase_ex
        elif what=='phase_em':
            col = colormap( phase_em_color[si] )
            fs[yi:yf,xi:xf] = s.phase_em
        elif what=='LS':
            col = colormap( LS_color[si] )
            fs[yi:yf,xi:xf] = s.LS
        elif what=='ET_ruler':
            col = colormap( ET_ruler_color[si] )
            fs[yi:yf,xi:xf] = s.ET_ruler
        elif what=='mean intensity':
            col = colormap( intensity_color[si] )
            fs[yi:yf,xi:xf] = intensity[si]
        elif what=='residual':
            col = colormap( residual_color[si] )
            fs[yi:yf,xi:xf] = residual[si]
            # print yi, yf
            # print xi, xf
            # print s.intensity[0]
            # import time
            # time.sleep(.5)

        else:
            print "Not sure how to interpret what=%s" % (what)
            return

        p = Rectangle((s.coords[0],s.coords[1]), s.width, s.height, \
                          facecolor=col, edgecolor=None, linewidth=0, alpha=1)
        ax.add_patch( p )
    
    np.savetxt(what+'data.txt', fs)

    ax.figure.canvas.draw()

    #plt.figure()
    #plt.hist( intensity )
    #plt.draw()





    