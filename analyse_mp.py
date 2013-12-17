import sys
import numpy as np
from util2dpolim.movie import Movie
from util2dpolim.misc import save_hdf5, combine_outputs, pixel_list
import time as stopwatch



prefix = 'x:\TDM3'
basename = 'TDM3-488-OD1-02'

# bounds in x,y format: (left column, upper row, right column, lower row) -- where 'upper' and 'lower' 
# correspond to the way the image is plotted (matrix-style, origin in the top left corner of the picture)
bgbounds   = [1,200,50,500]         #[110,405,400,450] 
fullbounds = [150,200,450,500]        #[110, 80,400,360]
resolution = 2
Nsplit     = 1
SNR    = 0
VFR    = .5
Nprocs = 2

topedges = np.arange(fullbounds[1], fullbounds[3], resolution )  
splittopedges = np.array_split( topedges, Nsplit )

for ste in splittopedges:

    m = Movie( prefix, basename )

    boolimage = np.ones( (m.sample_data.rawdata.shape[1],m.sample_data.rawdata.shape[2]), dtype=np.bool )*True
    blankfitexclusion = {'left':140, 'right':450, 'upper':180, 'lower':480, 'op':'exclude'}
    boolimage = pixel_list( m, blankfitexclusion, boolimage )
    m.fit_blank_image( boolimage, verbosity=0 )

    m.find_portraits( frameoffset=0 )
    m.find_lines()

    m.define_background_spot( bgbounds )
    
    for line in ste:
        for col in range( fullbounds[0], fullbounds[2], resolution ):
            bounds = [ col, line, col+resolution-1, line+resolution-1 ]
            m.define_spot( bounds )
    print 'nspots: ',len(m.spots)

    m.correct_excitation_intensities()
    m.correct_emission_intensities()

    m.are_spots_valid( SNR=SNR, validframesratio=VFR )
    #m.validspots = range(len(m.spots))

    if not len(m.validspots)==0:
        tstart = stopwatch.time()
        m.run_mp( Nprocs=Nprocs, fits=True, mods=True, ETruler=True )
        print stopwatch.time()-tstart

#        save_hdf5( m, myspots=np.arange(len(m.validspots)) )

