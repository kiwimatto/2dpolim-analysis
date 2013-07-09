import sys
import numpy as np
from util2dpolim.movie import Movie
from util2dpolim.misc import grid_image_section_into_squares_and_define_spots, update_image_files, show_mem
import time as stopwatch

from mpi4py import MPI
comm = MPI.COMM_WORLD
myrank = comm.Get_rank()
nprocs = comm.Get_size()

show_mem()

starttime = stopwatch.time()

prefix = '/home/kiwimatto/Desktop/Lund/2D/2dpolim-analysis/test_data__cool_new_setup_header/'

tstart = stopwatch.time()

# scan whole image, roughly 10 rows at a time
Nrowsatatime = 10
for r in np.arange(0,2,Nrowsatatime):
    m = Movie( prefix, 'moviename' )
    m.define_background_spot( [0,100,50,300] )

    bounds = [ 0, r, 512, r+Nrowsatatime ]
    if myrank==0: print 'current bounds: ',bounds

    grid_image_section_into_squares_and_define_spots( m, res=1, bounds=bounds )
    if myrank==0: print 'nspots: ',len(m.spots)

    m.collect_data()
    m.startstop()
    m.assign_portrait_data()
    #m.are_spots_valid(SNR=4)

print 'time taken: %fs' % (stopwatch.time()-starttime)

raise SystemExit

m.chew_a_bit(SNR=0)





# for si,s in enumerate(m.validspots):
#     print "si=%d\tLS=%f\tM_ex=%f\tM_em=%f" % (si, s.LS, s.M_ex, s.M_em)

print 'p=',myrank,': done. ',(stopwatch.time()-tstart)

comm.barrier()

stopwatch.sleep(myrank)
update_image_files(m, 'mean_intensity', fileprefix=prefix )
update_image_files(m, 'SNR', fileprefix=prefix )
update_image_files(m, 'M_ex', fileprefix=prefix )
update_image_files(m, 'M_em', fileprefix=prefix )
update_image_files(m, 'phase_ex', fileprefix=prefix )
update_image_files(m, 'phase_em', fileprefix=prefix )
update_image_files(m, 'LS', fileprefix=prefix )
update_image_files(m, 'r', fileprefix=prefix )
update_image_files(m, 'ET_ruler', fileprefix=prefix )



print "time taken: %fs" % (stopwatch.time()-tstart)





