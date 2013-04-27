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
        p = Rectangle((movie.bg.coords[0],movie.bg.coords[1]), movie.bg.width, movie.bg.height, \
                          facecolor=[1,0,0,0.1], edgecolor=[1,0,0,.75])
        ax.add_patch( p )

    # prepare intensities, etc
    inten = []
    Mex   = []
    for ss in movie.spots:
        inten.append( np.mean(ss.intensity) )
        Mex.append( ss.M_ex[0] )

    intensity = np.array(inten)
    Mex = np.array(Mex)

    intensity = (intensity-np.min(intensity))/np.max(intensity)                
    Mex = (Mex-np.min(Mex))/np.max(Mex)

    # make a patch for each spot
    for si,s in enumerate(movie.spots):

        # determine color (color axes 
        if what=='M_ex':
            col = colormap( Mex[si] )
        elif what=='M_em':
            col = colormap( s.M_em[0] )
        elif what=='phase_ex':
            # get all phases
            phases = []
            for ss in movie.spots:
                phases.append( ss.phase_ex[0] )
            pha = np.array(phases)
            # scale to range 0--1
            pha = (pha-np.min(pha))/np.max(pha)
            col = colormap( pha[si] )
        elif what=='phase_em':
            # get all phases
            phases = []
            for ss in movie.spots:
                phases.append( ss.phase_em[0] )
            pha = np.array(phases)
            # scale to range 0--1
            pha = (pha-np.min(pha))/np.max(pha)
            col = colormap( pha[si] )
        elif what=='LS':
            col = colormap( s.LS[0] )
        elif what=='ET_ruler':
            col = colormap( s.ET_ruler )
        elif what=='mean intensity':
            col = colormap( intensity[si] )
        else:
            print "Not sure how to interpret what=%s" % (what)
            return

        p = Rectangle((s.coords[0],s.coords[1]), s.width, s.height, facecolor=col, edgecolor=None, alpha=1)
        ax.add_patch( p )

    ax.figure.canvas.draw()

        
