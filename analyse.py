from util_2d import *
from util_misc import grid_image_section_into_squares_and_define_spots, show_spot_data

def show_mem():
    print "memory usage:"
    import memory
    print memory.memory()/(1024*1024)
    print memory.resident()/(1024*1024)

show_mem()

##movie = Movie( '01/AM01.SPE', '01/AMtest1_ex.txt', '01/AMtest1_em.txt' )
# movie.define_spot( [60,76,85,95], label='transmission blob' )
# movie.chew()
# movie.portraits[0].show_picture()


#prefix = '/media/sf_shared_with_VM/S3/S3-TQ1Filter/'
#prefix = '/home/kiwimatto/Desktop/130422/S3/'
prefix = '/home/kiwimatto/Desktop/Lund/2D/2dpolim-analysis/bla/'
          
global_phase = 9.0 * np.pi/180.0   # must be in radians!!!

#AMdegrees = [90]#,120,150,180]

#for am in AMdegrees:
for f in ['090','120','150','180']:
    m = Movie( prefix+"AM633-"+f+".SPE", prefix+"MS-AM633-"+f+".txt", \
                   phase_offset_excitation=global_phase, use_new_fitter=True )

    m.define_background_spot( [0,0,89,511] )
#m.define_background_spot( [0,6,7,7] )
    
    grid_image_section_into_squares_and_define_spots( m, res=100, bounds=[200,250,300,350] )
#grid_image_section_into_squares_and_define_spots( m, res=1, bounds=[0,0,1,1] )


    m.chew_a_bit()
#    print "LS=%f\tM_ex=%f\tM_em=%f" % (m.spots[0].LS, m.spots[0].M_ex, m.spots[0].M_em)

    show_mem()




raise SystemExit








# import matplotlib.pyplot as plt

# for i in range(4):
#     movies[i].show_average_picture()


prefix = '/home/kiwimatto/Desktop/Matthias Sample/'
global_phase = -59.0

movie = Movie( prefix+'Mov-02.SPE', prefix+'ex.txt', prefix+'em.txt', \
                   phase_offset_excitation=global_phase )
movie.define_background_spot( [0,95,248,103] )


class Storage:
    def __init__( self, movie ):
        self.pic=movie.averageportrait.copy() 
        self.M_ex=movie.M_ex 
        self.M_em=movie.M_em 
        self.phase_ex=movie.phase_ex 
        self.phase_em=movie.phase_em 
        self.LS=movie.LS 
        self.int = movie.spots[0].mean_I

import time as stopwatch
start = stopwatch.time()
Bigstore = []
for xi in leftedges:
    smallstore = []
    for yi in topedges:
        movie.define_spot( [yi,xi,yi-1+res,xi-1+res] )
        movie.chew(quiet=True)
        smallstore.append( Storage( movie ) )
        movie.spots = []
    Bigstore.append( smallstore )

print "time taken: %fs" % (stopwatch.time()-start)

# Bigstore can be indexed via x and y edges indices (in that order!)
mod_ex_matrix = np.zeros( (len(Bigstore[0]),len(Bigstore)) )
mod_em_matrix = np.zeros( (len(Bigstore[0]),len(Bigstore)) )
mod_LS_matrix = np.zeros( (len(Bigstore[0]),len(Bigstore)) )

for i in range( len(Bigstore) ):
    for j in range( len(Bigstore[0]) ):
        mod_ex_matrix[j,i] = Bigstore[i][j].M_ex
        mod_em_matrix[j,i] = Bigstore[i][j].M_em
        mod_LS_matrix[j,i] = Bigstore[i][j].LS

plt.matshow( np.mean( movie.camera_data.rawdata, axis=0 )[rb[0]:rb[2],rb[1]:rb[3]] )
