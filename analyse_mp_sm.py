import sys
import numpy as np
from util2dpolim.movie import Movie
from util2dpolim.misc import save_hdf5, save_spot_hdf5, pixel_list, import_spot_positions
import time as stopwatch


#### data directory and file name ####

#prefix = '/home/kiwimatto/Desktop/Lund/Aboma_Jaggs/Full 2D/exp1'
prefix = '/home/rafael/Desktop/Win/Crystal'
basename = 'S01-P01-01'

SNR        = 2
VFR        = .5
Nprocs     = 3

m = Movie( prefix, basename )
m.find_portraits( frameoffset=0 ) # it was 1 before
m.find_lines()


#### blank fitting ####

boolimage = np.ones( (m.sample_data.rawdata.shape[1],m.sample_data.rawdata.shape[2]), dtype=np.bool )*True
blankfitexclusion = {'left':140, 'right':450, 'upper':180, 'lower':480, 'op':'exclude'}
boolimage = pixel_list( m, blankfitexclusion, boolimage )
m.fit_blank_image( boolimage, verbosity=0 )


#### background spot definition ####

# bounds in x,y format: (left column, upper row, right column, lower row) -- where 'upper' and 'lower' 
# correspond to the way the image is plotted (matrix-style, origin in the top left corner of the picture)
# bgbounds   = [10,100,60,500]         #[110,405,400,450] 
# m.define_background_spot( bgbounds )



#### main part of the analysis ####

import_spot_positions( m, basename, 6, 'circle', use_borderbg=True  )

m.correct_excitation_intensities()
m.correct_emission_intensities()
m.are_spots_valid( SNR=SNR, validframesratio=VFR )

if not len(m.validspots)==0:
    print len(m.validspots)
    tstart = stopwatch.time()

    #### decide which part of the analysis you want done here ###

    m.run_mp( Nprocs=Nprocs, fits=True, mods=True, ETruler=True, ETmodel=True )
    print stopwatch.time()-tstart


#### write output contrasts into a text file --- this is only needed for the combination of multiple data sets ####
allp = []
for i,s in enumerate(m.validspots):
    allp.append( [s.mean_intensity, s.M_ex,s.M_em, s.phase_ex,s.phase_em,s.LS, s.ET_ruler, s.ETmodel_et,s.ETmodel_md_fu,s.ETmodel_th_fu] )
allp = np.array(allp)
np.savetxt('bla',allp)

# for d in dir(m.validspots[0]):
#     x = getattr(m.validspots[0],d)
#     if type(x)==np.ndarray:
#         print d+"\t"+str(x.shape)+"\t"+str(x.dtype)
#     else:
#         print d+str(type(x))


save_spot_hdf5( m )

#### save output to hdf5 ####
save_hdf5( m )


