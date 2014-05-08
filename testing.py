import numpy as np
import matplotlib.pyplot as plt
plt.interactive(True)
from scipy.optimize import fmin


prefix   = '/home/kiwimatto/Desktop/Lund/2D/'
basename = 'ETmodeltest'


#######################################################
#########  Simple simulation to test ET model #########
#######################################################

# We create two absorbing dipoles, orthogonal to each other. We can 
# control the modulation in excitation by scaling one dipole between 0 and 1.
# We also control the main absorption axis (angle of the unscaled dipole).

def orthogonal_dipoles( phi0, gr ):
    """This function returns two orthogonal dipoles: the
    first has unit length and points along phi0 [rad], the
    second is orthogonal to this (phi0-pi/2) and is scaled
    by the argument gr (expected to be between 0 and 1)."""

    dipoleA = np.array( [np.cos(phi0), np.sin(phi0), 0] )
    dipoleB = gr * np.cross( dipoleA, [0,0,1] )

    return dipoleA, dipoleB

# this in-line function converts the desired modulation depth into a scaler for
# the second dipole---which is passed as an argument to orthogonal_dipoles()
ratio_from_mod = lambda mod: np.sqrt((1-mod)/(1+mod))

# we fix our 'measurement' conditions
exangles = np.linspace(0, np.pi, 6, endpoint=False)
emangles = np.linspace(0, np.pi, 4, endpoint=False)
print exangles * 180/np.pi
print emangles * 180/np.pi

# create absorbing dipoles
M_ex     = .5
phase_ex = 30 * np.pi/180.0
absD1, absD2 = orthogonal_dipoles( phase_ex, gr=ratio_from_mod(M_ex) )

# create funnel dipoles
M_fu     = .9
phase_fu = 45 * np.pi/180.0
funD1, funD2 = orthogonal_dipoles( phase_fu, gr=ratio_from_mod(M_fu) )

# declare the transition matrix
absD1_to_funD1 = 0.1
absD1_to_funD2 = 0.9
absD2_to_funD1 = 0.0
absD2_to_funD2 = 0.3
T = np.array( [ [absD1_to_funD1, absD2_to_funD1], [absD1_to_funD2, absD2_to_funD2] ] )
# what's left? One minus the sum down the columns is the fraction that the absorbing 
# dipoles _didn't_ give away:
noT = 1-np.sum(T,axis=0)
# make sure we're not being silly
assert np.all(noT<1), "Dipoles can't give away more than they have!"

# calculate intensity frame-by-frame
intensities = []
for ema in emangles:
    for exa in exangles:
        # abs probabilites
        Pabs1 = np.dot( absD1, [np.cos(exa),np.sin(exa),0] )**2
        Pabs2 = np.dot( absD2, [np.cos(exa),np.sin(exa),0] )**2
        # emission probs of abs dipoles
        Pemi1 = np.dot( absD1, [np.cos(ema),np.sin(ema),0] )**2
        Pemi2 = np.dot( absD2, [np.cos(ema),np.sin(ema),0] )**2
        # emission probs of funnel dipoles
        Pfun1 = np.dot( funD1, [np.cos(ema),np.sin(ema),0] )**2
        Pfun2 = np.dot( funD2, [np.cos(ema),np.sin(ema),0] )**2

        # intensity from funnel
        ET   = np.dot( [Pfun1, Pfun2], np.dot( T, [Pabs1,Pabs2] ) )
        # intensity from abs dipoles
        noET = np.dot( noT, [Pemi1, Pemi2] )
        # total
        intensities.append( ET+noET )

print intensities



###################################################
#########  output for analysis software  ##########
###################################################

# this function writes a bogus motor file
def write_motor_file( prefix, basename, exangles, emangles ):
    header = """optical element in excitation: L/2 Plate
rotation mode in excitation: stepwise
rotation mode in emission: stepwise
phase offset in deg: 0.000000
full excitation power in mW: 1.0
optical density of filter: 1.000000
excitation power after filter: 1.0
real gain: 1.000000
camera: Cascade
laser wavelength in nm: 488.000000
objective: 60X oil
em correction modulation depth: 0.0
em correction phase: 0.000000
NA: 0.650000
user notes: 
This is a bogus header intended for use with test data.
END-OF-HEADER 
Frame	Excitation Physical Angle	Excitation Polarization Angle	Excitation Power [A.U.]	Emission Angle
"""
    fid = open( prefix+'MS-'+basename+'.txt', 'w' )
    fid.write(header)
    frame = 0
    for ema in emangles:
        for exa in exangles:
            fid.write('%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % \
                          (frame, exa*180/np.pi, exa*180/np.pi, 1.0, ema*180/np.pi))
            frame += 1
    fid.close()


# write the bogus motor file with the angles we've used
write_motor_file( prefix, basename, exangles, emangles )

# create a 2-pixel movie (one pixel for background)
data = np.zeros( (exangles.size*emangles.size,2,1) )
data[:,1,0] = intensities
np.save( prefix+basename, data )



###############################################################
#########  now we run the analysis software on this  ##########
###############################################################

from util2dpolim.movie import Movie
from util2dpolim.misc import *

m = Movie( prefix, basename )
m.find_portraits( frameoffset=0 )
m.find_lines()

m.define_background_spot( [0,0,0,0] )
m.define_spot( [0,1,0,1] )

m.correct_excitation_intensities()
m.correct_emission_intensities()
m.are_spots_valid( SNR=0, validframesratio=1 )

m.fit_all_portraits_spot_parallel_selective( range(len(m.validspots)) )
m.find_modulation_depths_and_phases_selective( range(len(m.validspots)) )
for s in m.validspots:
    s.values_for_ETruler( newdatalength=1024 )
m.ETmodel_selective( range(len(m.validspots)) )

remove_hdf5_files( m )    
save_spot_hdf5( m )
save_hdf5( m )







