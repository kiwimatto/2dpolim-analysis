import numpy as np
import h5py



filename = '/home/kiwimatto/Desktop/2D Jagg/Area 1 Blue 2D_spot_output.hdf5'



f = h5py.File(filename,'r')

for x in f.keys():
    if x.startswith('spot_'):
        spotnumber = np.int(x[5:])
        if f[x]['ok_by_eye'].value:
            center = np.array( f[x]['setup/center'] )
            print 'spot %d --- coordinates: %d,%d' % (spotnumber, center[0],center[1])

f.close()
