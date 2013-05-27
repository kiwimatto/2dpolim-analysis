import time
import numpy as np
import matplotlib.pyplot as plt
plt.interactive(True)
from datetime import datetime, timedelta
from mpl_toolkits.axes_grid1 import AxesGrid

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


def trim_noisy_data( movie, what='M_ex', threshold=None ):

    if what=='M_ex' or what=='M_em':
        if threshold==None:
            threshold = 3*movie.bg_spot.std
        quantity  = 'intensity_time_average'

    elif what=='phase_ex':
        if threshold==None:
            threshold = 0.1
        quantity  = 'M_ex'

    elif what=='phase_em':
        if threshold==None:
            threshold = 0.1
        quantity  = 'M_em'

    elif what=='LS':   # LS is a bit special here
        if not hasattr(movie.validspots[0],'phase_ex_trimmed'):
            trim_noisy_data( movie, what='phase_ex', threshold=None )
        if not hasattr(movie.validspots[0],'phase_em_trimmed'):
            trim_noisy_data( movie, what='phase_em', threshold=None ) 
        for si,s in enumerate(movie.validspots):
            setattr(s,'LS_trimmed', getattr(s,'phase_ex')-getattr(s,'phase_em') )
        return  # early return from this one

    # all others go through this one
    for si,s in enumerate(movie.validspots):
        if getattr(s,quantity) <= threshold:
            setattr(s,what+'_trimmed', np.nan )
        else:
            setattr(s,what+'_trimmed', getattr(s,what) )


def save_spot_data( movie, what='M_ex', whole_image=True, fileprefix='' ):
    # We assume in the following that the first spot (movie.spots[0]) is
    # the upper left spot in the image, with the origin (0,0) of the
    # image also being in the upper left corner of the plot.
    # xinit = movie.spots[0].coords[0]
    # yinit = movie.spots[0].coords[1]

    xdim=movie.spots[-1].coords[2]-movie.spots[0].coords[0]+1
    ydim=movie.spots[-1].coords[3]-movie.spots[0].coords[1]+1
    xinit = movie.spots[0].coords[0]
    yinit = movie.spots[0].coords[1]
    
    if whole_image:
        fs = np.ones( (movie.camera_data.datasize[1],movie.camera_data.datasize[2]) ) * np.nan
        xinit = 0
        yinit = 0
    else:
        fs = np.ones((ydim,xdim)) * np.nan
    
    for si,s in enumerate(movie.validspots):
        xi = s.coords[0]-xinit
        xf = s.coords[2]-xinit+1  # edges...
        yi = s.coords[1]-yinit
        yf = s.coords[3]-yinit+1
        fs[yi:yf,xi:xf] = getattr(s,what)

    np.savetxt(fileprefix + what + '_output_data.txt', fs)

        

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
    residuals = []

    # collect the data from all spots
    for ss in movie.validspots:
        mean_intensities.append( ss.intensity_time_average )
        Ms_ex.append( ss.M_ex )
        Ms_em.append( ss.M_em )
        phases_ex.append( ss.phase_ex )
        phases_em.append( ss.phase_em )
        LSs.append( ss.LS )
        ET_rulers.append( ss.ET_ruler )
        residuals.append( ss.residual )
    # turn the lists into numpy arrays
    intensity = np.array(mean_intensities)
    M_ex = np.array(Ms_ex)
    M_em = np.array(Ms_em)
    phase_ex = np.array(phases_ex)
    phase_em = np.array(phases_em)
    LS = np.array(LSs)
    ET_ruler = np.array(ET_rulers)
    residual = np.array(residuals)

    # rescale values to color range 
    intensity_color = (intensity-np.min(intensity))/np.max(intensity)                
    M_ex_color = (M_ex-np.min(M_ex))/np.max(M_ex)
    M_em_color = (M_em-np.min(M_em))/np.max(M_em)
    phase_ex_color = (phase_ex-np.min(phase_ex))/np.max(phase_ex)
    phase_em_color = (phase_em-np.min(phase_em))/np.max(phase_em)
    LS_color = (LS-np.min(LS))/np.max(LS)
    ET_ruler_color = (ET_ruler-np.min(ET_ruler))/np.max(ET_ruler)
    residual_color = (residual-np.min(residual))/np.max(residual)

    # edges! the +1 accounts for the fact that a spot includes its edges (really? does it?)
    xdim=movie.spots[-1].coords[2]-movie.spots[0].coords[0]+1
    ydim=movie.spots[-1].coords[3]-movie.spots[0].coords[1]+1
    fs=np.zeros((ydim,xdim))
    xinit = movie.spots[0].coords[0]
    yinit = movie.spots[0].coords[1]
     
    for si,s in enumerate(movie.validspots):
        xi = s.coords[0]-xinit
        xf = s.coords[2]-xinit+1  # edges...
        yi = s.coords[1]-yinit
        yf = s.coords[3]-yinit+1

        # determine color (color axes)
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
        elif what=='residual':
            col = colormap( residual_color[si] )
            fs[yi:yf,xi:xf] = residual[si]
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

    #plt.figure()
    #plt.hist( intensity )
    #plt.draw()



def run_self_test():
    import util_2d

    print "Creating test data set"
    create_test_data_set()
    print "Loading parameter file"
    p = np.load('testdataparams.npy')
    print "Analysing test data"
    m = util_2d.Movie( "testdata.npy", "testmotordata.txt", \
               phase_offset_excitation=0, use_new_fitter=True )
    grid_image_section_into_squares_and_define_spots( m, res=1, bounds=[0,0,16,16] )
    m.chew_a_bit()

    save_spot_data( m, what='M_ex', whole_image=True )
    save_spot_data( m, what='M_em', whole_image=True )
    save_spot_data( m, what='phase_ex', whole_image=True )
    save_spot_data( m, what='phase_em', whole_image=True )
    save_spot_data( m, what='LS', whole_image=True )

    mex = np.loadtxt('M_ex_output_data.txt')
    mem = np.loadtxt('M_em_output_data.txt')
    pex = np.loadtxt('phase_ex_output_data.txt')
    pem = np.loadtxt('phase_em_output_data.txt')
    ls  = np.loadtxt('LS_output_data.txt')

    fig = plt.figure(figsize=(18,12))
    grid = AxesGrid(fig, [.05,.05,.9,.9], # similar to subplot(132)
                    nrows_ncols = (3, 5),
                    axes_pad = 0.5,
                    share_all=True,
                    label_mode = "L",
                    cbar_location = "right",
                    cbar_mode="each",
                    cbar_pad=.04
                    )

    # M_ex
    im = grid[0].imshow(p[:,:,0], interpolation='nearest')
    grid.cbar_axes[0].colorbar(im)    
    im = grid[5].imshow(mex, interpolation='nearest')
    grid.cbar_axes[5].colorbar(im)
    im = grid[10].imshow(p[:,:,0]-mex, interpolation='nearest')
    grid.cbar_axes[10].colorbar(im)

    # M_em
    im = grid[1].imshow(p[:,:,1], interpolation='nearest')
    grid.cbar_axes[1].colorbar(im)
    im = grid[6].imshow(mem, interpolation='nearest')
    grid.cbar_axes[6].colorbar(im)
    im = grid[11].imshow(p[:,:,1]-mem, interpolation='nearest')
    grid.cbar_axes[11].colorbar(im)

    # phase_ex
    im = grid[2].imshow(p[:,:,2], interpolation='nearest')
    grid.cbar_axes[2].colorbar(im)
    im = grid[7].imshow(pex, interpolation='nearest')
    grid.cbar_axes[7].colorbar(im)
    im = grid[12].imshow(p[:,:,2]-pex, interpolation='nearest')
    grid.cbar_axes[12].colorbar(im)

    # phase_em
    im = grid[3].imshow(p[:,:,3], interpolation='nearest')
    grid.cbar_axes[3].colorbar(im)
    im = grid[8].imshow(pem, interpolation='nearest')
    grid.cbar_axes[8].colorbar(im)
    actual_emission_phase = p[:,:,2]+p[:,:,3] 
    im = grid[13].imshow((actual_emission_phase)-pem, interpolation='nearest')
    grid.cbar_axes[13].colorbar(im)

    # LS
    im = grid[4].imshow(p[:,:,2]-p[:,:,3], interpolation='nearest')
    grid.cbar_axes[4].colorbar(im)
    im = grid[9].imshow(pex-pem, interpolation='nearest')
    grid.cbar_axes[9].colorbar(im)
    actual_emission_phase = p[:,:,2]+p[:,:,3]
    im = grid[14].imshow(p[:,:,2]-actual_emission_phase-(pex-pem), interpolation='nearest')
    grid.cbar_axes[14].colorbar(im)
    
    grid[0].set_title('M_ex')
    grid[1].set_title('M_em')
    grid[2].set_title('phase_ex')
    grid[3].set_title('phase_em')
    grid[4].set_title('LS')

    grid[0].set_ylabel('test input')
    grid[5].set_ylabel('analysis output')
    grid[10].set_ylabel('difference')


    # a11 = f.add_subplot(3,5,1)
    # a12 = f.add_subplot(3,5,6)
    # a13 = f.add_subplot(3,5,11)
    # plt.colorbar()
    # a11.imshow(p[:,:,0])
    # a12.imshow(mex)
    # a13.imshow(p[:,:,0]-mex)
    plt.draw()


def create_test_data_set( illumination='flat', peakphotons=1000, noise=False, SNR=100, flat_bg=0, debug=False):

    Npixel_x = 16
    Npixel_y = 16

    X,Y = np.meshgrid( np.arange(Npixel_x,dtype=float), np.arange(Npixel_y,dtype=float) )
    if illumination=='flat':
        laserspot = np.ones_like(X)
    elif illumination=='2dgauss':
        # incident intensity distribution as a simple 2d gaussian
        pos_x = 3
        pos_y = 0
        sigma_x = sigma_y = 2
        laserspot = .5/(np.pi*sigma_x*sigma_y) * \
            np.exp( -.5*(X-pos_x)**2/sigma_x**2 -.5*(Y-pos_y)**2/sigma_y**2 )

        laserspot /= np.max(laserspot)
    else:
        raise ValueError('Not sure how to interpret illumination=%s ...' % illumination)

    laserspot *= peakphotons

    # angle-vs-time: slopes and steps and intervals etc..
    ex_angle_increment_per_sec  = 100.0  # deg/s
    em_angle_change_every_N_sec = 4.0    # s
    em_increment                = 22.5   # deg
    shutter_off_time            = .1     # s

    # time defs for movie
    Nframes = 500
    integration_time = .1                # s
    timer_step = .05                     # s

    # timer running at labview-output freq
    timer      = np.arange( 0, Nframes*integration_time, timer_step )
    # timer for the movie
    frametimes = np.arange( 0, Nframes*integration_time, integration_time )

    # init movie 
    data = np.zeros( (Nframes, Npixel_y, Npixel_x) )

    # init all parameters of the ET model
    md_ex    = np.outer( np.ones((Npixel_y,)), np.linspace(0,1,Npixel_x) )
    md_fu    = np.outer( np.linspace(0,1,Npixel_y), np.ones((Npixel_x,)) )
    phase_ex = (np.random.random(size=(Npixel_y,Npixel_x))-.5) * np.pi   # rad
    phase_fu = (np.random.random(size=(Npixel_y,Npixel_x))-.5) * np.pi   # rad
    gr       = np.random.random(size=(Npixel_y,Npixel_x))
    et       = np.ones((Npixel_y,Npixel_x))

    # ET model calculations
    alpha = 0.5 * np.arccos( .5*(((gr+2)*md_ex)-gr) )   # rad

    ph_ii_minus = phase_ex - alpha   # rad
    ph_ii_plus  = phase_ex + alpha   # rad
    
    # we generate the movie
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

        # store into data array
        data[i,:,:] = (et*Fet + (1-et)*Fnoet) * laserspot

        # # if flat illumination, then make room for a background spot in the top left 4x4 pixel
        # if illumination=='flat':
        #     data[i,:4,:4] = 0

        # add flat bg
        data[i,:,:] += flat_bg

        # add noise
        if noise:
            data[i,:,:] += np.random.normal( scale=peakphotons/SNR, size=(Npixel_y,Npixel_x) )

#        print et*Fet + (1-et)*Fnoet
#        print 'frame number %d done' % i

    # now we generate the motor data
    exa = np.zeros_like(timer)
    ema = np.zeros_like(timer)
    shutter = np.ones_like(timer, dtype=np.bool)
    for i in range(len(timer)):
#        print '.',
        exa[i] = ex_angle_increment_per_sec*timer[i]  # assuming we start at 0deg   # deg
        ema[i] = np.floor(timer[i]/em_angle_change_every_N_sec) * em_increment      # deg

        if np.mod(timer[i],em_angle_change_every_N_sec) <= shutter_off_time: 
            if not timer[i] <= shutter_off_time:  # stay on at start
                shutter[i] = 0

    if debug:
        print data.shape
        print timer.shape
        print exa.shape
        print ema.shape
        print exaframe.shape
        print emaframe.shape

    # write all this into files!
    writeTestDataMotorFile(timer,exa,ema,shutter)
    writeTestDataFile(data)
    writeTestDataParameters( md_ex, md_fu, phase_ex, phase_fu, gr, et )

    if debug:
        import matplotlib.pyplot as plt
        plt.interactive(True)    

        plt.plot( range(data.shape[0]), data[:,0,0] )
        plt.plot( range(data.shape[0]), exaframe/1000.0 )
        plt.plot( range(data.shape[0]), emaframe/1000.0 )

        plt.figure()
        plt.imshow( laserspot, interpolation='none' )

    return


def writeTestDataMotorFile(timer,exa,ema,shutter):
    towrite = ['Date       Time         Motor Em        Motor Ex        Shutter Status\n']
    starttime = time.time()
    for i in range(len(timer)):
        # construct line to print, first: data-time string
        line  = time.strftime( '%m/%d/%Y %H:%M:%S', time.localtime(starttime+timer[i]) )
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
    for s in movie.validspots:
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


    
