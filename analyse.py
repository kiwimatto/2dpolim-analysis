from util_2d import *
from util_misc import grid_image_section_into_squares_and_define_spots, show_spot_data, save_spot_data, update_image_files
import time as stopwatch

from mpi4py import MPI
comm = MPI.COMM_WORLD
myrank = comm.Get_rank()
nprocs = comm.Get_size()

def show_mem():
    print "memory usage:"
    import memory
    print memory.memory()/(1024*1024)
    print memory.resident()/(1024*1024)
show_mem()

##movie = Movie( '01/AM01.SPE', '01/AMtest1_ex.txt', '01/AMtest1_em.txt' )
# movie.define_spot( [60,76,85,95], label='transmission blob' )
# movie.chew()
# movie.portraits[0].show_picture()


#prefix = '/media/sf_shared_with_VM/S3/S3-TQ1Filter/'
#prefix = '/home/kiwimatto/Desktop/130422/S3/'
#prefix = '/home/kiwimatto/Desktop/130514/a-MEH/'
#prefix = '/home/rafael/Desktop/Win/LC/130527-LC/Ink/'
prefix = '/home/kiwimatto/Desktop/Lund/Experimental/130530/'
          
global_phase = 1.0 * np.pi/180.0   # must be in radians!!!

#AMdegrees = [90]#,120,150,180]

tstart = stopwatch.time()

m = Movie( prefix+"MEHPPV-SA1-OD2-488nm-03.SPE", prefix+"MS-MEHPPV-SA1-OD2-488nm-03.txt", \
               phase_offset_excitation=global_phase, which_setup='new setup', \
               excitation_optical_element='Polarizer')

m.define_background_spot( [0,450,330,450] )

fullbounds = [ 50,120,330,430 ]

if myrank==0: print 'fullbounds=',fullbounds
mybounds = [ fullbounds[0] +myrank*(fullbounds[2]-fullbounds[0])/nprocs, \
                 fullbounds[1], \
                 fullbounds[0] +(myrank+1)*(fullbounds[2]-fullbounds[0])/nprocs, \
                 fullbounds[3] ]
print 'p=',myrank,': mybounds=',mybounds

grid_image_section_into_squares_and_define_spots( m, res=8, bounds=mybounds )    
print 'p=',myrank,': nspots=',len(m.spots)


m.chew_a_bit(SNR=4)
#print "LS=%f\tM_ex=%f\tM_em=%f" % (m.spots[0].LS, m.spots[0].M_ex, m.spots[0].M_em)

print 'p=',myrank,': done. ',(stopwatch.time()-tstart)

comm.barrier()

import time
time.sleep(myrank)
#update_image_files(m, 'intensity_time_average', fileprefix=prefix )
update_image_files(m, 'M_ex', fileprefix=prefix )
update_image_files(m, 'M_em', fileprefix=prefix )
update_image_files(m, 'phase_ex', fileprefix=prefix )
update_image_files(m, 'phase_em', fileprefix=prefix )
update_image_files(m, 'ET_ruler', fileprefix=prefix )
update_image_files(m, 'LS', fileprefix=prefix )

save_spot_data(m, 'intensity_time_average', fileprefix=prefix )
# # save_spot_data(m, 'M_ex', fileprefix=prefix )
# # save_spot_data(m, 'M_em', fileprefix=prefix )
# # save_spot_data(m, 'phase_ex', fileprefix=prefix )
# # save_spot_data(m, 'phase_em', fileprefix=prefix )
# # save_spot_data(m, 'ET_ruler', fileprefix=prefix )
# # save_spot_data(m, 'LS', fileprefix=prefix )


# show_mem()

# print "time taken: %fs" % (stopwatch.time()-tstart)





