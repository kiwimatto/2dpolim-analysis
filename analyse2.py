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

prefix = '/home/rafael/Desktop/Win/TDM3/'
basename = 'TDM3-488-OD1-02'
bgbounds   = [10,155,95,480]
fullbounds = [100,155,480,480]
resolution = 2
Nsplit     = 15
rafaSNR    = 0


# bounds in x,y format: (left column, upper row, right column, lower row) -- where 'upper' and 'lower' 
# correspond to the way the image is plotted (matrix-style, origin in the top left corner of the picture)
topedges = np.arange(fullbounds[1], fullbounds[3], resolution )  
splittopedges = np.array_split( topedges, Nsplit )


### for single-molecule samples you can use this:

# m = Movie( prefix, basename )
# m.startstop()
# m.define_background_spot( bgbounds )
# import_spot_positions( m, 'coords.txt', boxedgelength=5 )



#for r in np.arange(fullbounds[1], fullbounds[3]-Nrowsatatime, Nrowsatatime):
for ste in splittopedges:
    print '-->Rafa you changed the BG substraction to have the right shape<--'
    m = Movie( prefix, basename )
    m.startstop()

    m.define_background_spot( bgbounds )
    
    for line in ste:
        for col in range( fullbounds[0], fullbounds[2], resolution ):
            bounds = [ col, line, col+resolution-1, line+resolution-1 ]
            m.define_spot( bounds )

#    if myrank==0: print 'current bounds: ',bounds

#    grid_image_section_into_squares_and_define_spots( m, res=resolution, bounds=bounds )
    if myrank==0: print 'nspots: ',len(m.spots)

    m.correct_excitation_intensities()
    m.correct_emission_intensities()   # .5, 45.0/180*np.pi )
    m.are_spots_valid( SNR=rafaSNR )

    # the rest is done only if we actually have any valid spots here
    if not len(m.validspots)==0:
        # m.fit_all_portraits_spot_parallel()
        # m.find_modulation_depths_and_phases()

        myspots = np.array_split( np.arange(len(m.validspots)), nprocs )
        m.fit_all_portraits_spot_parallel_selective( myspots[myrank] )
        m.find_modulation_depths_and_phases_selective( myspots[myrank] )
        m.ETrulerFFT_selective( myspots[myrank] )
        #    m.ETmodel_selective( myspots[myrank] )

        # all processes save their contributions separately
        save_hdf5( m, myspots[myrank], prefix, myrank )

print 'p=',myrank,': done. ',(stopwatch.time()-tstart)

comm.barrier()
stopwatch.sleep(1)

# first process gets to combine them into a single file
if myrank==0: 
    print "proc 0 here: combining output"
    combine_outputs( m.data_basename, prefix )

# all done







