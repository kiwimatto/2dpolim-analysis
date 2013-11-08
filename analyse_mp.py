import sys
import numpy as np
from util2dpolim.movie import Movie
from util2dpolim.misc import save_hdf5, combine_outputs
import time as stopwatch

prefix = '/home/kiwimatto/Desktop/130925 - MEHPPV YUXI/TDM5/'
basename = 'TDM5-488-OD1-01'

# bounds in x,y format: (left column, upper row, right column, lower row) -- where 'upper' and 'lower' 
# correspond to the way the image is plotted (matrix-style, origin in the top left corner of the picture)

bgbounds   = [20,30,100,500]         #[110,405,400,450] 
fullbounds = [130,30,400,500]        #[110, 80,400,360]
resolution = 8
Nsplit     = 1
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

    m.are_spots_valid( SNR=1 )

    if not len(m.validspots)==0:
        tstart = stopwatch.time()
        m.run_mp( Nprocs=Nprocs, fits=True, mods=True, ETruler=True )
        print stopwatch.time()-tstart

        save_hdf5( m, myspots=np.arange(len(m.validspots)) )

