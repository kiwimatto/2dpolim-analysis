import sys
from util2dpolim.movie import Movie
from util2dpolim.misc import grid_image_section_into_squares_and_define_spots, update_image_files, show_mem
import time as stopwatch

from mpi4py import MPI
comm = MPI.COMM_WORLD
myrank = comm.Get_rank()
nprocs = comm.Get_size()

show_mem()

prefix = '/home/kiwimatto/Desktop/Lund/2D/2dpolim-analysis/'

m = Movie( '.', 'S1' )         

tstart = stopwatch.time()
m.define_background_spot( [0,100,50,300] )
fullbounds = [ 300,160,400,440 ]

if myrank==0: print 'fullbounds=',fullbounds
mybounds = [ fullbounds[0] +myrank*(fullbounds[2]-fullbounds[0])/nprocs, \
                 fullbounds[1], \
                 fullbounds[0] +(myrank+1)*(fullbounds[2]-fullbounds[0])/nprocs, \
                 fullbounds[3] ]
print 'p=',myrank,': mybounds=',mybounds

grid_image_section_into_squares_and_define_spots( m, res=10, bounds=mybounds )
print 'p=',myrank,': nspots=',len(m.spots)

use_alternate_path=False
if use_alternate_path:
    m.collect_data()
    m.startstop()
    m.assign_portrait_data()
    for s in m.spots:
        s.check_if_valid(SNR=4)
        if s.isvalid:
#            s.collect_and_assign()
            s.cos_fit()
            s.findModDepths()
            print '.',
            sys.stdout.flush()
else:
    m.chew_a_bit(SNR=0)

# for si,s in enumerate(m.validspots):
#     print "si=%d\tLS=%f\tM_ex=%f\tM_em=%f" % (si, s.LS, s.M_ex, s.M_em)

print 'p=',myrank,': done. ',(stopwatch.time()-tstart)

comm.barrier()

stopwatch.sleep(myrank)
update_image_files(m, 'mean_intensity', fileprefix=prefix )
update_image_files(m, 'M_ex', fileprefix=prefix )
update_image_files(m, 'M_em', fileprefix=prefix )
update_image_files(m, 'phase_ex', fileprefix=prefix )
update_image_files(m, 'phase_em', fileprefix=prefix )
update_image_files(m, 'ET_ruler', fileprefix=prefix )
update_image_files(m, 'LS', fileprefix=prefix )
update_image_files(m, 'LS', fileprefix=prefix )

print "time taken: %fs" % (stopwatch.time()-tstart)





