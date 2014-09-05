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
