import numpy as np
import matplotlib.pyplot as plt
import time as stopwatch
plt.interactive(True)
from scipy.optimize import fmin

# 2d polim imports
from util2dpolim.movie import Movie
from util2dpolim.fitting import SFA_full_func
from util2dpolim.misc import write_motor_file, remove_hdf5_files, save_spot_hdf5, save_hdf5

if __name__ == '__main__':

    ######################################################
    ####### USER INPUTS ##################################
    # information about file to use for saving information, NOTE folder must exist !!!!
    prefix = 'C:\Users\MIPP\Desktop\Test\\'
    basename = 'ETmodeltest'

    # information for fit
    Nprocs = 1
    DoETruler = True  # False#
    DoETmodel = True  # False#

    # information for simulated spot
    # portrait
    ex_ang_number = 6
    em_ang_number = 4
    # absorbing dipoles
    M_ex = .0
    phase_ex = 20 * np.pi / 180.0
    # create funnel dipoles
    M_fu = .9
    phase_fu = phase_ex  # 45 * np.pi/180.0
    gr = 1.0
    et = 0.0
    #######################################################
    #### END of user inputs ###############################
    #######################################################



    #######################################################
    #########  Simple simulation to test ET model #########
    #######################################################

    phase_offset_in_deg = 0.0
    phase_offset_in_rad = phase_offset_in_deg * np.pi/180.0

    # we fix our 'measurement' conditions
    exangles_motorfile = np.linspace(0, np.pi, ex_ang_number, endpoint=False)
    exangles = exangles_motorfile + phase_offset_in_rad
    emangles = np.linspace(0, np.pi, em_ang_number, endpoint=False)
    print exangles * 180/np.pi
    print emangles * 180/np.pi

    EXA,EMA = np.meshgrid( exangles, emangles )
    intensity = SFA_full_func( [M_fu, phase_fu, gr, et], EXA, EMA, M_ex, phase_ex ).reshape((EXA.shape))
    print intensity.shape

    ###################################################
    #########  output for analysis software  ##########
    ###################################################

    # write the bogus motor file with the angles we've used --- and possibly cheat with the phase offset
    write_motor_file( prefix, basename, exangles_motorfile, emangles, phase_offset=0 ) # phase_offset_in_deg )

    # create a 2-pixel movie (one pixel for background)
    data = np.zeros( (exangles.size*emangles.size,2,1) )
    data[:,1,0] = intensity.flatten()
    np.save( prefix+basename, data )


    ###############################################################
    #########  now we run the analysis software on this  ##########
    ###############################################################

    bgbounds = [0,0,0,0]
    bounds = [0,1,0,1]
    m = Movie( prefix, basename )
    m.find_portraits( frameoffset=0 )
    m.find_lines()

    m.define_background_spot( bgbounds )
    m.define_spot( bounds )
    print 'nspots: ', len(m.spots)

    m.correct_excitation_intensities()
    m.correct_emission_intensities()

    m.are_spots_valid( SNR=0, validframesratio=1 )

    if not len(m.validspots) == 0:
        tstart = stopwatch.time()
        m.fit_all_portraits_spot_parallel_selective(range(len(m.validspots)))
        m.find_modulation_depths_and_phases_selective(range(len(m.validspots)))
        for s in m.validspots:
            s.values_for_ETruler(newdatalength=1024)
        m.ETmodel_selective(range(len(m.validspots)))

        # m.run_mp(Nprocs=Nprocs, fits=True, mods=True, ETruler=DoETruler, ETmodel=DoETmodel)
        print stopwatch.time() - tstart
    print 'Done with analysis, now we save the hdf5 file'

    save_spot_hdf5(m)
    save_hdf5( m )

    print "Mex  = %f" % (m.validspots[0].M_ex)
    print "phex = %f" % (m.validspots[0].phase_ex * 180./np.pi)
    print "Mem  = %f" % (m.validspots[0].M_em)
    print "phem = %f" % (m.validspots[0].phase_em * 180./np.pi)
    print "LS   = %f" % (m.validspots[0].LS * 180./np.pi)
    print "et   = %f" % (m.validspots[0].ETmodel_et)
    print "anisotropy 2D          = %f" % m.validspots[0].anisotropy
    print "anisotropy classical   = %f" % m.validspots[0].anisotropy_n