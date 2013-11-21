import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import os

plt.interactive(True)

prefix = '/home/kiwimatto/Desktop/130925 - MEHPPV YUXI/TDM3/'
basename = 'TDM3-488-OD1-02'

print 'Requested file is: %s' % prefix+basename+'_output.hdf5'
f = h5py.File( prefix+basename+'_output.hdf5' )
os.chdir( prefix )

images     = ['mean_intensity','M_ex','M_em','phase_ex','phase_em','ET_ruler']
multiplier = [1, 1, 1, 180.0/np.pi, 180.0/np.pi, 1 ]
colormap   = [cm.jet, cm.jet, cm.jet, cm.hsv, cm.hsv, cm.jet]

bounds = [0,511,0,511]
if bounds==None:
    bounds=[0,512,0,512]


for i,im in enumerate(images):

    print '--> ',i
    meanint = np.array( f['images/mean_intensity_image'] )[ bounds[0]:bounds[1], bounds[2]:bounds[3] ]
    meanint /= np.max(meanint)
    bla = np.array( f['images/'+im+'_image'] )[ bounds[0]:bounds[1], bounds[2]:bounds[3] ]

    plt.figure()
    abc = plt.imshow( bla * multiplier[i], cmap=colormap[i], alpha=1.0 )
    # print abc._A.shape
    # abc._A[:,:,3] *= meanint
    plt.draw()

    plt.colorbar()
    plt.title( im )
    plt.savefig( im+'_image.png')

