import time
import numpy as np
from datetime import datetime, timedelta

def deal_with_date_time_string( motorobj, datetimestring ):        
    """This function converts the date+time string into the difference time 
    (in seconds) since the start of the experiment. The first value passed 
    is treated in a special way: it sets the zero time mark.

    It is the same function for both motors, that's why it's here,  
    outside the class definitions.
    """
    dt = datetime.strptime(datetimestring, "%m/%d/%Y %H:%M:%S.%f")
    if motorobj.experiment_start_datetime == None:
        motorobj.experiment_start_datetime = dt
        return 0.0
    else:
        td = dt - motorobj.experiment_start_datetime
        return td.total_seconds()
    # return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 1e6 



def grid_image_section_into_squares_and_define_spots( movie, res, bounds ):
    
    rb = rectangular_blob = bounds #[80,56,155,89]  # pixel indices (starting from zero!)

    leftedges = range(rb[1],rb[3],res)
    if leftedges[-1]+res > rb[3]:
        leftedges = leftedges[:-1]

    topedges = range(rb[0],rb[2],res)
    if topedges[-1]+res > rb[2]:
        topedges = topedges[:-1]

    for xi in leftedges:
        for yi in topedges:
            movie.define_spot( [yi,xi,yi-1+res,xi-1+res] )

    return 



def show_spot_data( movie, what='M_ex', which_cmap=None, show_bg_spot=True ):

    import matplotlib.pyplot as plt
    import matplotlib.cm as cmap
    from matplotlib.patches import Rectangle

    if which_cmap==None:
        colormap = cmap.jet

    plt.figure()
    
    # draw average as background, use colormap gray
#    plt.imshow( np.mean( movie.camera_data.rawdata, axis=0 ), cmap=cmap.gray )
    plt.imshow( movie.camera_data.average_image, cmap=cmap.gray )

    ax = plt.gca()

    if show_bg_spot:
        p = Rectangle( (movie.bg_spot.coords[0],movie.bg_spot.coords[1]), \
                           movie.bg_spot.width, movie.bg_spot.height, \
                           facecolor=[1,0,0,0.1], edgecolor=[1,0,0,.75])
        ax.add_patch( p )

    # prepare intensities, etc
    mean_intensities = []
    Ms_ex     = []
    Ms_em     = []
    phases_ex = []
    phases_em = []
    LSs       = []
    ET_rulers = []
    
    for ss in movie.spots:
        mean_intensities.append( ss.intensity_time_average )
        Ms_ex.append( ss.M_ex )
        Ms_em.append( ss.M_em )
        phases_ex.append( ss.phase_ex )
        phases_em.append( ss.phase_em )
        LSs.append( ss.LS )
        ET_rulers.append( ss.ET_ruler )

    intensity = np.array(mean_intensities)
    M_ex = np.array(Ms_ex)
    M_em = np.array(Ms_em)
    phase_ex = np.array(phases_ex)
    phase_em = np.array(phases_em)
    LS = np.array(LSs)
    ET_ruler = np.array(ET_rulers)

    # rescale values to color range 
    intensity_color = (intensity-np.min(intensity))/np.max(intensity)                
    M_ex_color = (M_ex-np.min(M_ex))/np.max(M_ex)
    M_em_color = (M_em-np.min(M_em))/np.max(M_em)
    phase_ex_color = (phase_ex-np.min(phase_ex))/np.max(phase_ex)
    phase_em_color = (phase_em-np.min(phase_em))/np.max(phase_em)
    LS_color = (LS-np.min(LS))/np.max(LS)
    ET_ruler_color = (ET_ruler-np.min(ET_ruler))/np.max(ET_ruler)

    # edges! the +1 accounts for the fact that a spot includes its edges
    xdim=movie.spots[-1].coords[2]-movie.spots[0].coords[0]+1
    ydim=movie.spots[-1].coords[3]-movie.spots[0].coords[1]+1
    fs=np.zeros((ydim,xdim))
    xinit = movie.spots[0].coords[0]
    yinit = movie.spots[0].coords[1]
     
    for si,s in enumerate(movie.spots):
        xi = s.coords[0]-xinit
        xf = s.coords[2]-xinit+1  # edges...
        yi = s.coords[1]-yinit
        yf = s.coords[3]-yinit+1

        # determine color (color axes 
        if what=='M_ex':
            col = colormap( M_ex_color[si] )
            fs[yi:yf,xi:xf] = s.M_ex
        elif what=='M_em':
            col = colormap( M_em_color[si] )
            fs[yi:yf,xi:xf] = s.M_em
        elif what=='phase_ex':
            col = colormap( phase_ex_color[si] )
            fs[yi:yf,xi:xf] = s.phase_ex
        elif what=='phase_em':
            col = colormap( phase_em_color[si] )
            fs[yi:yf,xi:xf] = s.phase_em
        elif what=='LS':
            col = colormap( LS_color[si] )
            fs[yi:yf,xi:xf] = s.LS
        elif what=='ET_ruler':
            col = colormap( ET_ruler_color[si] )
            fs[yi:yf,xi:xf] = s.ET_ruler
        elif what=='mean intensity':
            col = colormap( intensity_color[si] )
            fs[yi:yf,xi:xf] = intensity[si]
            # print yi, yf
            # print xi, xf
            # print s.intensity[0]
            # import time
            # time.sleep(.5)

        else:
            print "Not sure how to interpret what=%s" % (what)
            return

        p = Rectangle((s.coords[0],s.coords[1]), s.width, s.height, \
                          facecolor=col, edgecolor=None, linewidth=0, alpha=1)
        ax.add_patch( p )
    
    np.savetxt(what+'data.txt', fs)

    
    ax.figure.canvas.draw()



def create_test_data_set():

    Npixel_x = 8
    Npixel_y = 8

    # incident intensity distribution as a simple 2d gaussian
    pos_x = 3
    pos_y = 0
    sigma_x = sigma_y = 2
    X,Y = np.meshgrid( np.arange(Npixel_x,dtype=float), np.arange(Npixel_y,dtype=float) )
    laserspot = .5/(np.pi*sigma_x*sigma_y) * \
        np.exp( -.5*(X-pos_x)**2/sigma_x**2 -.5*(Y-pos_y)**2/sigma_y**2 )

    laserspot /= np.max(laserspot)

    maxsignal  = 1000
    laserspot *= maxsignal
    # bg         = 2
    # laserspot += bg

    ex_angle_increment_per_sec  = 100.0  # deg/s
    em_angle_change_every_N_sec = 4.0    # s
    em_increment                = 22.5   # deg
    shutter_off_time            = .1     # s

    Nframes = 500
    integration_time = .1                # s
    timer_step = .05                     # s

    timer      = np.arange( 0, Nframes*integration_time, timer_step )
    frametimes = np.arange( 0, Nframes*integration_time, integration_time )

    data = np.zeros( (Nframes, Npixel_x, Npixel_y) )

    md_ex = np.random.random(size=(Npixel_y,Npixel_x))
    md_fu = np.random.random(size=(Npixel_y,Npixel_x))
    phase_ex = (np.random.random(size=(Npixel_y,Npixel_x))-.5) * np.pi   # rad
    phase_fu = (np.random.random(size=(Npixel_y,Npixel_x))-.5) * np.pi   # rad
    gr = np.random.random(size=(Npixel_y,Npixel_x))
    et = np.random.random(size=(Npixel_y,Npixel_x))

    alpha = 0.5 * np.arccos( .5*(((gr+2)*md_ex)-gr) )   # rad

    ph_ii_minus = phase_ex - alpha   # rad
    ph_ii_plus  = phase_ex + alpha   # rad
    
    exaframe = np.zeros( (Nframes,) )
    emaframe = np.zeros( (Nframes,) )
    for i in range(len(frametimes)):
        ex = ex_angle_increment_per_sec*frametimes[i]   # deg
        em = np.floor(frametimes[i]/em_angle_change_every_N_sec) * em_increment  # deg
        exaframe[i] = ex   # deg
        emaframe[i] = em   # deg
        ex *= np.pi/180.0    # rad
        em *= np.pi/180.0    # rad

        Fnoet  =    np.cos( ex-ph_ii_minus )**2 * np.cos( em-ph_ii_minus )**2
        Fnoet += gr*np.cos( ex-phase_ex )**2    * np.cos( em-phase_ex )**2
        Fnoet +=    np.cos( ex-ph_ii_plus )**2  * np.cos( em-ph_ii_plus )**2
        
#        print Fnoet

        Fnoet /= (2+gr)
#        Fnoet /= np.max(Fnoet)

        Fet   = .25 * (1+md_ex*np.cos(2*(ex-phase_ex))) \
            * (1+md_fu*np.cos(2*(em-phase_fu-phase_ex)))
#        Fet  /= np.sum(Fet)
#        print et*Fet
#        print (1-et)*Fnoet

        data[i,:,:] = (et*Fet + (1-et)*Fnoet) * laserspot
#        print et*Fet + (1-et)*Fnoet
        print 'frame number %d done' % i

    exa = np.zeros_like(timer)
    ema = np.zeros_like(timer)
    shutter = np.ones_like(timer, dtype=np.bool)
    for i in range(len(timer)):
        print '.',
        exa[i] = ex_angle_increment_per_sec*timer[i]  # assuming we start at 0deg   # deg
        ema[i] = np.floor(timer[i]/em_angle_change_every_N_sec) * em_increment      # deg

        if np.mod(timer[i],em_angle_change_every_N_sec) <= shutter_off_time: 
            if not timer[i] <= shutter_off_time:  # stay on at start
                shutter[i] = 0


    # print data.shape
    # print timer.shape
    # print exa.shape
    # print ema.shape
    # print exaframe.shape
    # print emaframe.shape

    writeTestDataMotorFile(timer,exa,ema,shutter)
    writeTestDataFile(data)
    writeTestDataParameters( md_ex, md_fu, phase_ex, phase_fu, gr, et )

    import matplotlib.pyplot as plt
    plt.interactive(True)    

    plt.plot( range(data.shape[0]), data[:,0,0] )
    plt.plot( range(data.shape[0]), exaframe/1000.0 )
    plt.plot( range(data.shape[0]), emaframe/1000.0 )

    plt.figure()
    plt.imshow( laserspot, interpolation='none' )

    return data, exaframe, emaframe


def writeTestDataMotorFile(timer,exa,ema,shutter):
    towrite = ['Date       Time         Motor Em        Motor Ex        Shutter Status\n']
    starttime = time.time()
    for i in range(len(timer)):
        # construct line to print, first: data-time string
        line  = time.strftime( '%d/%m/%Y %H:%M:%S', time.localtime(starttime+timer[i]) )
        line += '.%02d\t' % np.int(100*np.mod(starttime+timer[i],1))   # add some fractional seconds
        line += '%E\t' % ema[i]       #
        line += '%E\t' % exa[i]       # 
        if shutter[i]:
            line += 'open' 
        else:
            line += 'close'
        line += '\n'
        towrite.append( line )

    f = open('testmotordata.txt','w')
    f.writelines( towrite )
    f.close()

def writeTestDataFile(data):
    np.save( 'testdata.npy', data )

def writeTestDataParameters( md_ex, md_fu, phase_ex, phase_fu, gr, et ):
    arr = np.dstack( [md_ex, md_fu, phase_ex, phase_fu, gr, et] )
    np.save( 'testdataparams.npy', arr )

def compareTestParamsWithOutput( movie, paramfilename ):
    
    import matplotlib.pyplot as plt
    import matplotlib.cm as cmap
    from matplotlib.patches import Rectangle

    params = np.load( paramfilename )

    plt.imshow( params[:,:,0], cmap=cmap.gray )
    ax = plt.gca()

#    mexs = []
    for s in movie.spots:
#        mexs.append( s.M_ex )
        md_ex = params[s.coords[0],s.coords[1],0]
        s.M_ex_diff = s.M_ex - md_ex
        print 's.M_ex=%f\tmd_ex=%f\tdiff=%f' % (s.M_ex, md_ex, s.M_ex_diff)

        col = cmap.jet(s.M_ex)
        print col
        p = Rectangle((s.coords[0],s.coords[1]), s.width, s.height, \
                          facecolor=col, edgecolor=None, linewidth=0, alpha=.3)
        
        ax.add_patch( p )

    plt.draw()




def generate_single_funnel_test_data( excitation_angles, emission_angles, \
                                          md_ex=0, md_fu=1, \
                                          phase_ex=0, phase_fu=0, \
                                          gr=1.0, et=1.0 ):

    ex, em = np.meshgrid( excitation_angles, emission_angles )

    alpha = 0.5 * np.arccos( .5*(((gr+2)*md_ex)-gr) )

    ph_ii_minus = phase_ex - alpha
    ph_ii_plus  = phase_ex + alpha
    
    print ph_ii_minus
    print ph_ii_plus

    Fnoet  =    np.cos( ex-ph_ii_minus )**2 * np.cos( em-ph_ii_minus )**2
    Fnoet += gr*np.cos( ex-phase_ex )**2    * np.cos( em-phase_ex )**2
    Fnoet +=    np.cos( ex-ph_ii_plus )**2  * np.cos( em-ph_ii_plus )**2
        
    Fnoet /= (2+gr)
    
    Fet   = .25 * (1+md_ex*np.cos(2*(ex-phase_ex))) \
        * (1+md_fu*np.cos(2*(em-phase_fu-phase_ex)))
    
    
    Fem = et*Fet + (1-et)*Fnoet


    import matplotlib.pyplot as plt
    plt.interactive(True)
    plt.matshow( Fem, origin='bottom' )
    plt.colorbar()


    
