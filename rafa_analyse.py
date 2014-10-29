import os
import sys
import numpy as np
from util2dpolim.movie import Movie
from util2dpolim.misc import save_hdf5, combine_outputs, pixel_list#, import_spot_positions
import time as stopwatch

prefix = '/home/rafa/Desktop/share/OB 3month'
filelist = os.listdir( prefix )

# bounds in x,y format: (left column, upper row, right column, lower row) -- where 'upper' and 'lower' 
# correspond to the way the image is plotted (matrix-style, origin in the top left corner of the picture)
bgbounds   = [460,25,505,500]         #[110,405,400,450] 
fullbounds = [25,30,455,500]        #[110, 80,400,360]
resolution = 1
Nsplit     = 20
SNR    = 3
VFR    = .6
Nprocs = 3
blankfit   = False#True#
exclusion_bf = [28,85,442,500] # fullbounds # for most cases we can use fullbounds
DoETruler = True#False#
DoETmodel = True#False# 

topedges = np.arange(fullbounds[1], fullbounds[3], resolution )  
splittopedges = np.array_split( topedges, Nsplit )


for x in filelist:
    a = not x.startswith('blank')
    print x.endswith('.SPE')
    print a
    if (x.endswith('.SPE') and a):
        print x[:-4]
        basename = x[:-4]
        for ste in splittopedges:
            print '-------------current y been considered--------------'
            print ste
#            print basename
            m = Movie( prefix, basename )

            m.find_portraits( frameoffset=0 )
            m.find_lines()
            #### blank fitting ####
	    if blankfit:
    		boolimage = np.ones( (m.sample_data.rawdata.shape[1],m.sample_data.rawdata.shape[2]),
                                      dtype=np.bool )*True
    		blankfitexclusion = {'left':exclusion_bf[0], 'right':exclusion_bf[2],
                                     'upper':exclusion_bf[1], 'lower':exclusion_bf[3], 'op':'exclude'}
    		boolimage = pixel_list( m, blankfitexclusion, boolimage )
    		m.fit_blank_image( boolimage, verbosity=0 )
		m.define_background_spot( bgbounds )
		print 'done with BG fit and substraction'
            else:
                m.define_background_spot( bgbounds )
                print 'done with BG substraction'
            for line in ste:
                for col in range( fullbounds[0], fullbounds[2], resolution ):
                    bounds = [ col, line, col+resolution-1, line+resolution-1 ]
                    m.define_spot( bounds )
            print 'nspots: ',len(m.spots)

            m.correct_excitation_intensities()
            m.correct_emission_intensities()

            m.are_spots_valid( SNR=SNR, validframesratio=VFR )
            
            if not len(m.validspots)==0:
                tstart = stopwatch.time()
                m.run_mp( Nprocs=Nprocs, fits=True, mods=True, ETruler=DoETruler, ETmodel=DoETmodel )
                print stopwatch.time()-tstart

            save_hdf5( m )


