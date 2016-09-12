import sys
import numpy as np
from util2dpolim.movie import Movie
from util2dpolim.misc import * # save_hdf5, save_spot_hdf5, remove_hdf5_files, pixel_list, import_spot_positions
import time as stopwatch

######### USER INPUTS ######### USER INPUTS ######### USER INPUTS #########
### data directory and file name ####
## prefix is the data directory
prefix = '/home/rafa/Desktop/share/160205-TQ1'
#basename is the name of the SPE movie without the .SPE
basename = '1028-B2'

### Spot valudation parameters ###
SNR               = 1.0
VFR               = .4

### Background related input
##is there a blank to fit?
blankfit          = False    # True 
# exclusion_bf      = [45,85,420,495]  # only needed when blankfit is True
##are you going to use border as background
border_background = False
##if not then tell me what to use as background area
bgbounds          = [485,20,505,450] #only needed when border_background is False
                                     # values are: [left upper right lower]

######### END OF USER INPUTS ######### END OF USER INPUTS #########

m = Movie( prefix, basename )
m.find_portraits( frameoffset=0 )
m.find_lines()

#### blank fitting ####
if blankfit:
    boolimage = np.ones( (m.sample_data.rawdata.shape[1],m.sample_data.rawdata.shape[2]), dtype=np.bool )*True
    blankfitexclusion = {'left':exclusion_bf[0], 'right':exclusion_bf[2],
                                     'upper':exclusion_bf[1], 'lower':exclusion_bf[3], 'op':'exclude'}
    boolimage = pixel_list( m, blankfitexclusion, boolimage )
    m.fit_blank_image( boolimage, verbosity=0 )

#### now we have to import the spots. In here is very important to know if
#### the spots uses its borders as BG or not.

if border_background:
    print 'BG is defined automatically around the spot'
    import_userdefinedspot( m, prefix,  basename, 'UserDefined', use_borderbg=True )

else:
    assert type(bgbounds) == list, 'bgbound has to be a list'
    assert len(bgbounds) == 4, 'bgbound has to be a list of 4 numbers [left upper right lower]'
    m.define_background_spot( bgbounds )
    print 'BG was defined by user'
    import_userdefinedspot( m, prefix,  basename, 'UserDefined', use_borderbg=False )
	

#### main part of the analysis #### This next line takes care of BG definition due to use_borderbg

m.correct_excitation_intensities()
m.correct_emission_intensities()
m.are_spots_valid( SNR=SNR, validframesratio=VFR )

if not len(m.validspots)==0:
    print len(m.validspots)
    tstart = stopwatch.time()

    #### decide which part of the analysis you want done here ###

    # multi-proc (Linux only):
    # m.run_mp( Nprocs=Nprocs, fits=True, mods=True, ETruler=True, ETmodel=True )

    # single-core:
    m.fit_all_portraits_spot_parallel_selective( range(len(m.validspots)) )
    m.find_modulation_depths_and_phases_selective( range(len(m.validspots)) )
    for s in m.validspots:
        s.values_for_ETruler( newdatalength=1024 )
    m.ETmodel_selective( range(len(m.validspots)) )

#### save output to hdf5 ####
# For sm data, we don't update, but replace hdf5 files.
# So we delete first,
remove_hdf5_files( m )    
# and then write the new ones:
save_spot_hdf5( m )
save_hdf5( m )

# we're done here







