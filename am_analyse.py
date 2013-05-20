import sys
from util_2d import *

spefilename   = sys.argv[1]
motorfilename = sys.argv[2]
global_phase  = np.float( sys.argv[3] )
bg_c_1        = np.int( sys.argv[4] )
bg_c_2        = np.int( sys.argv[5] )
bg_c_3        = np.int( sys.argv[6] )
bg_c_4        = np.int( sys.argv[7] )
sig_c_1       = np.int( sys.argv[8] )
sig_c_2       = np.int( sys.argv[9] )
sig_c_3       = np.int( sys.argv[10] )
sig_c_4       = np.int( sys.argv[11] )
SNR           = np.int( sys.argv[12] )

bg_coords  = [bg_c_1, bg_c_2, bg_c_3, bg_c_4]
sig_coords = [sig_c_1, sig_c_2, sig_c_3, sig_c_4]

#raise SystemExit
m = Movie( spefilename, motorfilename, phase_offset_excitation=global_phase*np.pi/180.0, use_new_fitter=True )

m.define_background_spot( bg_coords )
m.define_spot( sig_coords )
m.chew_AM(SNR=SNR)
m.validspots[0].export_averagematrix('spotmatrix.npy')
