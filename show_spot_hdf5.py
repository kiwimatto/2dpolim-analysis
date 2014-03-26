import os, os.path, time
import numpy as np
import matplotlib.pyplot as plt
plt.interactive(True)
import h5py


prefix = os.path.normpath( '/home/kiwimatto/Desktop/2D Jagg' )+os.path.sep
basename = 'Area 1 Red 2D'


# read the stuff...
filename = prefix+basename+'_spot_output.hdf5'
assert os.path.isfile( filename ), 'File not found!'
fid = h5py.File( filename, 'r' )


# find how many spots we got
keys = fid.keys()
Nspots = 0
for k in fid.keys():
    if k.startswith('spot_'):
        Nspots += 1

print "Nspots = %d" % Nspots

# go through all spots and collect what we need into nice lists
M_ex = []
M_em = []
for k in fid.keys():
    if k.startswith('spot_'):
        M_ex.append( np.array(fid[k+'/contrasts/M_ex']) )
        M_em.append( np.array(fid[k+'/contrasts/M_em']) )

M_ex = np.array(M_ex)
M_em = np.array(M_em)
  
plt.plot( M_ex, M_em, 'ro' )
plt.plot( [0,1], [0,1], 'k:' )
plt.xlim( 0,1.1 )
plt.ylim( 0,1.1 )

