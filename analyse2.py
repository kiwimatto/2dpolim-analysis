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


### for single-molecule samples you can use this:

# m = Movie( prefix, basename )
# m.startstop()
# m.define_background_spot( bgbounds )
# import_spot_positions( m, 'coords.txt', boxedgelength=5 )

# prefix = '/home/kiwimatto/Desktop/Lund/2D/2dpolim-analysis/test_emission_correction/'
# basename = 'S3-633ex-OD1-gain1-675LP-01'
prefix = '/home/kiwimatto/Desktop/Lund/Experimental/130719__new_fredrik_samples/'
basename = 'edge1'

prefix = '/home/kiwimatto/Desktop/130925 - MEHPPV YUXI/TDM5/'
basename = 'TDM5-488-OD1-01'

# bounds in x,y format: (left column, upper row, right column, lower row) -- where 'upper' and 'lower' 
# correspond to the way the image is plotted (matrix-style, origin in the top left corner of the picture)

bgbounds   = [20,30,100,500]         #[110,405,400,450] 
fullbounds = [130,30,400,500]        #[110, 80,400,360]
resolution = 2
Nsplit     = 3

topedges = np.arange(fullbounds[1], fullbounds[3], resolution )  
splittopedges = np.array_split( topedges, Nsplit )

#for r in np.arange(fullbounds[1], fullbounds[3]-Nrowsatatime, Nrowsatatime):
for ste in splittopedges:

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

    m.are_spots_valid( SNR=3 )

    m.test_mp( Nprocs=3 )
    stopwatch.sleep(1)
    raise SystemExit


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
        save_hdf5( m, myspots[myrank], myrank )

print 'p=',myrank,': done. ',(stopwatch.time()-tstart)

comm.barrier()

# first process gets to combine them into a single file
if myrank==0: combine_outputs( m.data_basename, prefix )

# all done







