import sys
import numpy as np
from util2dpolim.movie import Movie
from util2dpolim.misc import save_hdf5, pixel_list, import_spot_positions
import time as stopwatch


prefix = '/home/kiwimatto/Desktop/130925 - MEHPPV YUXI/TDM3/'
basename = 'TDM3-488-OD1-02'

# bounds in x,y format: (left column, upper row, right column, lower row) -- where 'upper' and 'lower' 
# correspond to the way the image is plotted (matrix-style, origin in the top left corner of the picture)
bgbounds   = [1,200,50,500]         #[110,405,400,450] 
SNR    = 1
VFR    = .5
Nprocs = 4
blankfitexclusion = {'left':150, 'right':450, 'upper':200, 'lower':500, 'op':'exclude'}


m = Movie( prefix, basename )
m.find_portraits( frameoffset=1 )
m.find_lines()

#### blank fitting ####
boolimage = np.ones( (m.sample_data.rawdata.shape[1],m.sample_data.rawdata.shape[2]), dtype=np.bool )*True
boolimage = pixel_list( m, blankfitexclusion, boolimage )
print m.portrait_indices
m.fit_blank_image( boolimage, verbosity=0 )

m.define_background_spot( bgbounds )

import_spot_positions( m, 'coords.txt', 4, 'circle' )

m.correct_excitation_intensities()
m.correct_emission_intensities()

m.are_spots_valid( SNR=SNR, validframesratio=VFR )

if not len(m.validspots)==0:
    tstart = stopwatch.time()
    m.run_mp( Nprocs=Nprocs, fits=True, mods=True, ETruler=True, ETmodel=True )
    print stopwatch.time()-tstart

save_hdf5( m, myspots=np.arange(len(m.validspots)) )

