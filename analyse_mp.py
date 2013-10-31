import sys
import numpy as np
from util2dpolim.movie import Movie
from util2dpolim.misc import save_hdf5, combine_outputs
import time as stopwatch



prefix = '/home/rafael/Desktop/Win/Well01/'
basename = 'S1-W1-Area01'

# bounds in x,y format: (left column, upper row, right column, lower row) -- where 'upper' and 'lower' 
# correspond to the way the image is plotted (matrix-style, origin in the top left corner of the picture)
bgbounds   = [125,185,175,415]
fullbounds = [180,185,340,415]
resolution = 1
Nsplit     = 10
rafaSNR    = 6
Nprocs     = 4

topedges = np.arange(fullbounds[1], fullbounds[3], resolution )  
splittopedges = np.array_split( topedges, Nsplit )

for ste in splittopedges:

    m = Movie( prefix, basename )
    m.startstop()

    m.define_background_spot( bgbounds )
    
    for line in ste:
        for col in range( fullbounds[0], fullbounds[2], resolution ):
            bounds = [ col, line, col+resolution-1, line+resolution-1 ]
            m.define_spot( bounds )
    print 'nspots: ',len(m.spots)

    m.correct_excitation_intensities()
    m.correct_emission_intensities()

    m.are_spots_valid( SNR=rafaSNR )

    if not len(m.validspots)==0:
        tstart = stopwatch.time()
        m.test_mp( Nprocs=Nprocs )
        print stopwatch.time()-tstart

        save_hdf5( m, myspots=np.arange(len(m.validspots)) )

