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

prefix = '/home/rafael/Desktop/Win/Well01/'
basename = 'S1-W1-Area02'

# bounds in x,y format: (left column, upper row, right column, lower row) -- where 'upper' and 'lower' 
# correspond to the way the image is plotted (matrix-style, origin in the top left corner of the picture)
bgbounds   = [120,190,180,400]
fullbounds = [190,190,340,400]
resolution = 1
Nsplit     = 10
rafaSNR    = 7
VFRrafa    = .5
testrun    = False   #False/True

topedges = np.arange(fullbounds[1], fullbounds[3], resolution )  
splittopedges = np.array_split( topedges, Nsplit )


### for single-molecule samples you can use this:

# m = Movie( prefix, basename )
# m.startstop()
# m.define_background_spot( bgbounds )
# import_spot_positions( m, 'coords.txt', boxedgelength=5 )


#for r in np.arange(fullbounds[1], fullbounds[3]-Nrowsatatime, Nrowsatatime):
for iste,ste in enumerate(splittopedges):

    print "=== Now dealing with slice %d of %d ===" % (iste, Nsplit)
    m = Movie( prefix, basename )
    m.startstop()
    m.define_background_spot( bgbounds )
    
    for line in ste:
        for col in range( fullbounds[0], fullbounds[2], resolution ):
            bounds = [ col, line, col+resolution-1, line+resolution-1 ]
            m.define_spot( bounds )

    if myrank==0: print 'nspots: ',len(m.spots)

    m.correct_excitation_intensities()
    m.correct_emission_intensities()   # .5, 45.0/180*np.pi )
    m.are_spots_valid( SNR=rafaSNR, validframesratio=VFRrafa )

    # the rest is done only if we actually have any valid spots here
    if not len(m.validspots)==0:
        # m.fit_all_portraits_spot_parallel()
        # m.find_modulation_depths_and_phases()

        myspots = np.array_split( np.arange(len(m.validspots)), nprocs )

        if not testrun:
            m.fit_all_portraits_spot_parallel_selective( myspots[myrank] )
            m.find_modulation_depths_and_phases_selective( myspots[myrank] )
            m.ETrulerFFT_selective( myspots[myrank] )
        #    m.ETmodel_selective( myspots[myrank] )

        # all processes save their contributions separately
        save_hdf5( m, myspots[myrank], myrank )

print 'p=',myrank,': done. ',(stopwatch.time()-tstart)

comm.barrier()
stopwatch.sleep(1)

# first process gets to combine them into a single file
if myrank==0: 
    print "proc 0 here: combining output"
    print m.data_directory
    combine_outputs( m )

# all done







