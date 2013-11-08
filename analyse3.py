import sys
import numpy as np
from util2dpolim.movie import Movie
from util2dpolim.misc import grid_image_section_into_squares_and_define_spots, update_image_files, show_mem
import time as stopwatch

from mpi4py import MPI
comm = MPI.COMM_WORLD
myrank = comm.Get_rank()
nprocs = comm.Get_size()


prefix = '/home/rafael/Desktop/Win/130627-CT/S1/'

tstart = stopwatch.time()

m = Movie( prefix, 'S1-633ex-OD1-gain1-675LP-01' )
m.define_background_spot( [85,430,400,460] )

fullbounds = [ 85, 47, 400, 420 ]
if myrank==0: print 'fullbounds=',fullbounds

mybounds = [ fullbounds[0] +myrank*(fullbounds[2]-fullbounds[0])/nprocs, \
                 fullbounds[1], \
                 fullbounds[0] +(myrank+1)*(fullbounds[2]-fullbounds[0])/nprocs, \
                 fullbounds[3] ]
print 'p=',myrank,': mybounds=',mybounds

grid_image_section_into_squares_and_define_spots( m, res=2, bounds=mybounds )
print 'p=',myrank,': nspots=',len(m.spots)

m.chew_a_bit(SNR=5)

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





