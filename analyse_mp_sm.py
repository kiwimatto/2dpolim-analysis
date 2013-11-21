import sys
import numpy as np
from util2dpolim.movie import Movie
from util2dpolim.misc import save_hdf5, pixel_list, import_spot_positions
import time as stopwatch


prefix = '/home/rafael/Desktop/Win/TDM5/'
basename = 'TDM5-488-OD106-02'


# bounds in x,y format: (left column, upper row, right column, lower row) -- where 'upper' and 'lower' 
# correspond to the way the image is plotted (matrix-style, origin in the top left corner of the picture)

bgbounds   = [10,154,60,502]
SNR        = 3
VFR        = .5
Nprocs = 2

m = Movie( prefix, basename )
m.find_portraits( frameoffset=0 )
m.find_lines()
m.define_background_spot( bgbounds )

#### blank fitting ####
print m.sample_data.rawdata.shape
boolimage = np.ones( (m.sample_data.rawdata.shape[1],m.sample_data.rawdata.shape[2]), dtype=np.bool )*True
blankfitexclusion = {'left':140, 'right':450, 'upper':180, 'lower':450, 'op':'exclude'}
boolimage = pixel_list( blankfitexclusion, boolimage )
m.fit_blank_image( verbosity=0 )

import_spot_positions( m, 'coords.txt', 4, 'circle' )

m.correct_excitation_intensities()
m.correct_emission_intensities()
m.are_spots_valid( SNR=SNR, validframesratio=VFR )

if not len(m.validspots)==0:
    tstart = stopwatch.time()
    m.run_mp( Nprocs=Nprocs, fits=True, mods=True, ETruler=True )
    print stopwatch.time()-tstart

save_hdf5( m, myspots=np.arange(len(m.validspots)) )

