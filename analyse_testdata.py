import sys
import numpy as np
from util2dpolim.movie import Movie
from util2dpolim.misc import grid_image_section_into_squares_and_define_spots, \
    save_hdf5, \
    combine_outputs, \
    show_mem, \
    import_spot_positions
import time as stopwatch

#import_spot_positions( movie, coords_filename, boxedgelength=5 ):

from mpi4py import MPI
comm = MPI.COMM_WORLD
myrank = comm.Get_rank()
nprocs = comm.Get_size()

show_mem()
tstart = stopwatch.time()

prefix = '/home/kiwimatto/Desktop/Lund/2D/2dpolim-analysis/olle_simu/'
basename = 'testdata_stretched'

m = Movie( prefix, basename )
m.find_portraits( frameoffset=0 )
m.find_lines()

m.define_background_spot( [16,15,16,15] )

for i in np.arange(21)*31+15:
    spot = m.define_spot( [16,i,16,i] )
    print spot
    print spot.pixel," --> ",spot.intensity

m.correct_excitation_intensities()
m.correct_emission_intensities()   # .5, 45.0/180*np.pi )
m.are_spots_valid( SNR=3 )
    
if not len(m.validspots)==0:
    myspots = np.array_split( np.arange(len(m.validspots)), nprocs )
        
    m.fit_all_portraits_spot_parallel_selective( myspots[myrank] )
    m.find_modulation_depths_and_phases_selective( myspots[myrank] )
    m.ETrulerFFT_selective( myspots[myrank] )
    m.ETmodel_selective( myspots[myrank] )
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    for si,s in enumerate(m.validspots):
        et = s.ET_ruler
        if si==3:
            s.values_for_ETruler( newdatalength=1024, showfft=True )
        else:
            s.values_for_ETruler( newdatalength=1024, showfft=False )
        
        if not et==0:
            diff = (et-s.ET_ruler)/et
        else:
            diff = np.inf
        print "mex: %f, mem: %f, ruler: %f, model: %f " % (s.M_ex, s.M_em, s.ET_ruler, s.ETmodel_et )
        
#    save_hdf5( m, myspots[myrank], prefix, myrank )
else:
    print "-------------->  no valid spots..."  

print 'p=',myrank,': done. ',(stopwatch.time()-tstart)
    
raise SystemExit


comm.barrier()

# first process gets to combine them into a single file
#if myrank==0: combine_outputs( m.data_basename, prefix )

# all done



journal = np.loadtxt('journal.txt')
stretch = journal[:,0]
M   = journal[:,1]
Mr  = journal[:,2:-3]
Mex = journal[:,-3]
Mem = journal[:,-2]
r   = journal[:,-1]

print "stre\tM_ex\tM_em\tph_ex\tph_em\tLS\tr\tETr"

dat = []
for i,s in enumerate(m.validspots):
    print "%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t" % \
        (stretch[i], s.M_ex, s.M_em, s.phase_ex, s.phase_em, s.LS, s.anisotropy, s.ET_ruler)
    dat.append( [stretch[i], s.M_ex, s.M_em, s.phase_ex, s.phase_em, s.LS, s.anisotropy, s.ET_ruler] )

dat = np.array(dat)

print M
print Mr
print Mex
print Mem


import matplotlib.pyplot as plt
plt.interactive(True)

mmax = []
mmin = []
msum = []
for i,s in enumerate(stretch):
    mmax.append( np.max( m.validspots[i].recover_average_portrait_matrix() ) )
    mmin.append( np.min( m.validspots[i].recover_average_portrait_matrix() ) )
    print np.sum( m.validspots[i].recover_average_portrait_matrix() )

# for i,s in enumerate(stretch):
#     plt.imshow( m.validspots[i].recover_average_portrait_matrix(), \
#                     interpolation='nearest', origin='lower', vmax=np.max(mmax), vmin=np.min(mmin) )
#     plt.colorbar()
#     plt.savefig('portrait_%d.png' % i)
#     plt.clf()
# raise SystemExit


figwidth  = 12
figheight = 4
fig = plt.figure(figsize=(figwidth,figheight))

b = np.linspace( 0, 2*np.pi, 100 )
radius = .33
for i,s in enumerate(stretch):
    ax1 = fig.add_subplot( 1,3,1, position=[.05, .1, .1, .8] )
    ax2 = fig.add_subplot( 1,3,2, position=[.2, .1, .4, .8] )
    ax3 = fig.add_subplot( 1,3,3, position=[.6, .1, .4, .8] )

    ax1.plot( radius*np.cos(b), s*radius*np.sin(b), 'r' )
    ax1.set_xlim(-2.03, 2.03)
    ax1.set_ylim(-5.50, 5.50)

    ax2.plot( dat[:,0],dat[:,1],'b' )
    ax2.plot( dat[:,0],dat[:,2],'g' )
    ax2.plot( dat[:,0],dat[:,6],'k--' )
    ax2.plot( dat[:,0],dat[:,7],'r' )
#    ax2.plot( stretch,M,'b--' )
#    ax2.plot( stretch,r,'k--' )
    ax2.axvline( dat[i,0], color='red', ls=':' )
    ax2.set_ylim( -.05, 1.05 )
    ax2.set_xlabel('stretch factor')
    ax2.legend(('Mex','Mem','r','ET'), loc='lower right')

    im = ax3.imshow( m.validspots[i].recover_average_portrait_matrix(), \
                    interpolation='nearest', origin='lower', vmax=np.max(mmax), vmin=np.min(mmin) )
    plt.colorbar( im, ax=ax3 )

    plt.savefig('fig%03d.png' % i)
    plt.clf()

import os
os.system('convert -delay 50 fig*.png fig.gif')

raise SystemExit

figwidth  = 8
figheight = 4
fig = plt.figure()#figsize=(figwidth,figheight))

ax1 = fig.add_subplot( 1,1,1 )

ax1.plot( dat[:,0],dat[:,1],'b' )
ax1.plot( dat[:,0],dat[:,2],'g' )
ax1.plot( dat[:,0],dat[:,6],'k--' )
ax1.plot( dat[:,0],dat[:,7],'r' )

#ax1.set_ylim( -.05, 1.05 )
#ax1.set_xlabel('stretch factor')
ax1.legend(('Mex','Mem','r','ET'), loc='right')

ax1.set_xlabel('number of energy hops')
ax1.set_ylabel('Mex, Mem, r, ET')

# im = ax2.imshow( m.validspots[i].recover_average_portrait_matrix(), \
#                      interpolation='nearest', origin='lower', vmax=np.max(mmax), vmin=np.min(mmin) )
# plt.colorbar( im, ax=ax2 )

plt.savefig('fig.png')
#plt.clf()

# import os
# os.system('convert -delay 50 fig*.png fig.gif')


# fig.clf()
# ax = fig.add_subplot(111)
# ax.plot( dat[:,6], dat[:,1], 'o' )
# ax.plot( dat[:,6], M, 'ro' )
# plt.show()


