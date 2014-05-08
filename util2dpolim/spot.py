import numpy as np
import scipy.optimize as so
import scipy.ndimage as sn
from itertools import chain
#from scipy.spatial import ConvexHull
from matplotlib.patches import Rectangle, Polygon, Circle
# from matplotlib.collections import PatchCollection
# import matplotlib.cm as cmap


class Spot:
    def __init__( self, \
                      shape, \
                      int_type, \
                      label, \
                      parent, \
                      is_bg_spot=False, \
                      which_bg=None, \
                      use_blank=True, \
                      use_exspot=False, \
                      use_borderbg=False ):
                  
        """
        There's something noteworthy (speak: important) about the coordinates
        which define the box that is the 'spot'. First of all, the convention
        is [left, bottom, right, top], where the coordinates (0,0) define the
        _lower left corner_ of the image. If you plot an image with matplotlib's 
        imshow() or matshow(), then the origin is in the _top left_ (the way 
        you would write a matrix on paper) --- so there's potential for confusion
        here. Secondly, the coordinates given are interpreted as inclusive 
        boundaries, so that [3,3,5,5] gives a 3x3 spot from which intensities 
        are computed. This explains the occurence of various +1s in the code 
        below...



        We have potentially three data sources which determine the spot intensity:
        sample data, blank data, and excitation spot data.
        The data flow from raw inputs to spot intensity is supposed to be as follows:

        1) Define a background area. This should be outside the exposure slit, 
           capturing the dark noise only, since all other background sources will
           typically have a spatial dependence and need to be corrected through 
           blank data.
        2) For all three data sets, sample, blank, and excitation spot, subtract the 
           background mean intensity from the rest of the data. For sample and blank
           this happens on a per-frame basis, while for excitation spot data this is done
           with the average frame.
        3) Subtract bg-corrected blank from bg-corrected sample.
        4) Apply bg-corrected excitation spot data to blank-corrected sample:
           a) use for threshold-based selection of region of interest, and/or
           b) divide by the excitation spot intensity to yield brightness information

        """

        self.shape    = shape
        self.label    = label
        self.parent   = parent
        self.int_type = int_type
        self.set_spot_center()
        self.create_pixel_list( shape )
        self.is_bg_spot = is_bg_spot
        self.which_bg   = which_bg
        self.use_blank  = use_blank
        self.use_exspot = use_exspot
        self.use_borderbg = use_borderbg

        self.have_bg_sample  = None
        self.have_bg_blank   = None
        self.have_bg_exspot  = None
        self.have_blank      = None
        self.have_exspot     = None

        self.value_for_coverage_image = 1
        self.value_for_coverage_valid_image = 2
        self.value_for_coverage_ownbg_image = 100

#        self.width  = coords[2] -coords[0] +1
#        self.height = coords[3] -coords[1] +1


        # Fit parameter storage
        # dimensions: Nportraits x Nlines x 4
        # Last dimension encodes [phase, I0, M, residual]
        self.linefitparams     = np.nan*np.ones( (self.parent.Nportraits, self.parent.Nlines, 4) )
        self.verticalfitparams = np.nan*np.ones( (self.parent.Nportraits, \
                                                      self.parent.excitation_angles_grid.size, 4) )


        # if this spot is a background spot, then this is a special case
        if is_bg_spot:
            # if which_bg=='sample':
            #     data = self.parent.sample_data.rawdata[:, coords[1]:coords[3]+1, coords[0]:coords[2]+1 ]
            # elif which_bg=='blank':
            #     data = self.parent.blank_data.rawdata[:, coords[1]:coords[3]+1, coords[0]:coords[2]+1 ]
            # elif which_bg=='exspot':
            #     data = self.parent.exspot_data.rawdata[:, coords[1]:coords[3]+1, coords[0]:coords[2]+1 ]
            # else:
            #     raise ValueError('Spot class: do not understand which_bg="%s"' % which_bg)

            self.intensity = self.calculate_spot_intensity( which_data=which_bg, int_type=int_type )
            self.std       = self.calculate_spot_std( which_data=which_bg, int_type=int_type )
            # print "BG spot std:"
            # print self.std
            
        else:
            # here is all the normal stuff (when this spot is _not_ a bg spot)

            self.have_bg_sample    = (not self.parent.bg_spot_sample==None)
            self.have_bg_blank     = (not self.parent.bg_spot_blank==None)
            self.have_bg_exspot    = (not self.parent.bg_spot_exspot==None)
            self.have_blank  = (not self.parent.blank_data==None)
            self.have_exspot = (not self.parent.exspot_data==None)

            # spot knows its own background!
            if use_borderbg:
                Isample = self.calculate_spot_intensity( which_data='sample', int_type=int_type )
                self.dilate( borderwidth=2 )
                # print Isample
                # print self.borderbg
                # print Isample-self.borderbg
                # print "--------------------------------------"
                Isample -= self.borderbg
                self.intensity_type = int_type
                self.intensity      = Isample
                self.mean_intensity = np.mean(Isample)

                # record the spot in the mean intensity image
                for p in self.pixel:
                    self.parent.mean_intensity_image[ p[0], p[1] ] = self.mean_intensity
                return

            # grab blank data
            if self.have_blank and use_blank: 
#                blank  = self.parent.blank_data.rawdata[:, coords[1]:coords[3]+1, coords[0]:coords[2]+1 ]
                Iblank = self.calculate_spot_intensity( which_data='blank', int_type=int_type )
            # and exspot data
            if self.have_exspot and use_exspot: 
#                exspot  = self.parent.exspot_data.rawdata[:, coords[1]:coords[3]+1, coords[0]:coords[2]+1 ]
                Iexspot = self.calculate_spot_intensity( which_data='exspot', int_type=int_type )
            # and sample data
#            sample  = self.parent.sample_data.rawdata[:, coords[1]:coords[3]+1, coords[0]:coords[2]+1 ]
            Isample = self.calculate_spot_intensity( which_data='sample', int_type=int_type )

            if self.have_bg_blank:
                Iblank  -= self.parent.bg_spot_blank.intensity
            if self.have_bg_exspot:
                Iexspot -= self.parent.bg_spot_exspot.intensity
            if self.have_bg_sample:
                Isample -= self.parent.bg_spot_sample.intensity

            # at this point, all available data (sample, blank, exspot) has been collected and
            # and bg-corrected (if a bg-spot was defined in the parent class). 
            # subtract blank from sample
            if self.have_blank and use_blank:
#                Iblank2 =  np.concatenate((Iblank,Iblank))
#                print Iblank2
                Isample -= Iblank

            # divide by exspot power if that is required
            if self.have_exspot and use_exspot:
                Isample /= Iexspot
                
            self.intensity_type = int_type
            self.intensity      = Isample
            self.mean_intensity = np.mean(Isample)

            # record the spot in the mean intensity image
            for p in self.pixel:
                self.parent.mean_intensity_image[ p[0], p[1] ] = self.mean_intensity


    def enclosing_polygon( self ):
        # find edges (Note that the edge 'orientation' is important! Top edges point from left to right,
        # bottom edges from right to left!)
        edges = []
        for p in self.pixel:
            if not self.has_pixel( (p[0]-1,p[1]) ):
                edges.append( ((p[0]-.5,p[1]-.5), (p[0]-.5,p[1]+.5)) )   # top edge
            if not self.has_pixel( (p[0]+1,p[1]) ):
                edges.append( ((p[0]+.5,p[1]+.5), (p[0]+.5,p[1]-.5)) )   # bottom edge
            if not self.has_pixel( (p[0],p[1]-1) ):
                edges.append( ((p[0]+.5,p[1]-.5), (p[0]-.5,p[1]-.5)) )   # left edge
            if not self.has_pixel( (p[0],p[1]+1) ):
                edges.append( ((p[0]-.5,p[1]+.5), (p[0]+.5,p[1]+.5)) )   # right edge
        # print edges

        # extract vertices
        vertices = []
        for e in edges:
            vertices.append( e[0] )
            vertices.append( e[1] )

        vertices = np.array( list( set(vertices) ) )
        #print vertices

        # build connection information (here is where the edge orientation becomes important)
        conn = []
        # start at (arbitrary) first vertex
        conn.append( (vertices[0][0],vertices[0][1]) )
        #print conn

        while len(edges) > 0:
            for ei,e in enumerate( edges ):             # go through all remaining edges (again)
                if e[0]==(conn[-1][0],conn[-1][1]):     # if last connection vertex is same as edge start
                    conn.append( e[1] )                 # append edge ending vertex to conn
                    edges.pop(ei)                       # remove this edge from list of edges

        # we've been thinking in x- and y-coordinates, but what we need is row-column:
        conn = np.array(conn)
        conn[:,[0, 1]] = conn[:,[1, 0]]

        # create polygon and return
        polygon = Polygon( conn, closed=True, facecolor='red', edgecolor='red', alpha=.3, zorder=7 )
        return polygon


    def dilate( self, borderwidth=1 ):
        borderwidth = np.int(borderwidth)

        # we first find the corner coordinates (we could look them up in the shape dict for Rectangles,
        # but not generally for all shapes)
        minr = np.inf
        minc = np.inf
        maxr = -np.inf
        maxc = -np.inf
        for p in self.pixel:
            if p[0]<minr: minr=p[0]
            if p[0]>maxr: maxr=p[0]
            if p[1]<minc: minc=p[1]
            if p[1]>maxc: maxc=p[1]
                
        # now prepare a binary bitmap which fits the spot, plus
        # a border of borderwidth pixel surrounding it
        b = np.zeros( (maxr-minr+1 +2*borderwidth, maxc-minc+1 +2*borderwidth), dtype=np.bool )
        # mark the spot pixel in the bitmap
        for p in self.pixel:
            b[ p[0]-minr+borderwidth, p[1]-minc+borderwidth ] = True

        # now dilate the spot
        if self.shape['type']=='Rectangle':
            structure = sn.generate_binary_structure(2, 2)
        else:
            structure = None
        c = sn.binary_dilation(b, structure=structure, iterations=borderwidth)
        
        # bg is difference between original and dilated spot
        bgbitmap = c-b

        # find bg pixel positions, and add them to bgpixel list (unless 
        # we have a bg exclusion map, and that pixel is on it)
        bgpixel = []
        for ri in range(maxr-minr+1 +2*borderwidth):
            for ci in range(maxc-minc+1 +2*borderwidth):
                if bgbitmap[ri,ci]:
                    # if hasattr(self.parent,'bgexclusionmap'):
                    #     if self.parent.bgexclusionmap[ri+minr-1, ci+minc-1]==False:
                    #         bgpixel.append( (ri+minr-1, ci+minc-1) )
                    # else:
                    bgpixel.append( (ri+minr-2, ci+minc-2) )



        self.bgpixel = bgpixel
#        print 'Number of border pixel: %d' % len(bgpixel)

        data  = np.array( [self.parent.sample_data.rawdata[:,p[0],p[1]] for p in self.bgpixel] ).T
        self.borderbg    = np.sum( data, axis=1 )/float(len(self.bgpixel))
        self.borderbgstd = np.std( data, axis=1 )
        
        return bgbitmap


    def check_exclusion_zones( self, share_bg_pixel=True ):
        # # make a flat list of all pixel of all spots
        # allpixel   = list( chain.from_iterable([s.pixel for s in self.parent.spots]) )
        # allbgpixel = list( chain.from_iterable([s.bgpixel for s in self.parent.spots]) )

        newbgps = []
        for bgp in self.bgpixel:
            isalreadytaken = False
            for s in self.parent.spots:
                if s.has_pixel( bgp ):
                    isalreadytaken = True
                    break
                if not share_bg_pixel:
                    if s.has_pixel( bgp, checkbgpixel=True ):
                        isalreadytaken = True
                        break
            if not isalreadytaken:
                newbgps.append( bgp )

        self.bgpixel = newbgps

        Isample = self.calculate_spot_intensity( which_data='sample', int_type=int_type )
        # print Isample
        # print self.borderbg
        # print Isample-self.borderbg
        # print "--------------------------------------"
        data  = np.array( [self.parent.sample_data.rawdata[:,p[0],p[1]] for p in self.bgpixel] ).T
        self.borderbg    = np.sum( data, axis=1 )/float(len(self.pixel))
        self.borderbgstd = np.std( data, axis=1 )
        Isample -= self.borderbg
        self.intensity_type = int_type
        self.intensity      = Isample
        self.mean_intensity = np.mean(Isample)

        # record the spot in the mean intensity image
        for p in self.pixel:
            self.parent.mean_intensity_image[ p[0], p[1] ] = self.mean_intensity
        return


    def graphical_representation( self, color='red' ):
        r = None
        if self.shape['type']=='Rectangle':
            r = Rectangle( (self.shape['left']-.5, self.shape['upper']-.5), \
                               self.shape['right'] - self.shape['left'] + 1, \
                               self.shape['lower'] - self.shape['upper'] + 1, \
                               facecolor=color, edgecolor=color, alpha=.3, zorder=7 )
#            r = self.enclosing_polygon()

        elif self.shape['type']=='Circle':
            r = self.enclosing_polygon()
#            r = Circle( self.shape['center'], radius=self.shape['radius'], zorder=7 )
#            center = shape['center']
#            radius = shape['radius']
#            r = PatchCollection()

        else:
            raise Shitstorm        

        r.set_gid( "spot_%s" % self.shape['type'] )

        return r

    def has_pixel( self, pixel, checkbgpixel=False ):
        if checkbgpixel:
            for p in self.bgpixel:
                if p==pixel: return True
        else:
            for p in self.pixel:
                if p==pixel: return True
        return False

    def set_spot_center( self ):
        if self.shape['type']=='Rectangle':
            self.center = ( (self.shape['lower'] + self.shape['upper'])/2.0, \
                                (self.shape['right'] + self.shape['left'])/2.0 )
        elif self.shape['type']=='Circle':
            self.center = self.shape['center']
        else:
            raise fuckOff

    def create_pixel_list( self, shape ):
        if shape['type']=='Rectangle':
            left  = shape['left']
            right = shape['right']
            lower = shape['lower']
            upper = shape['upper']
            # validate coords:
            assert lower >= upper   # lower in the image means larger row number!
            assert (left >= 0) and (upper >= 0)
            assert (right<self.parent.sample_data.datasize[2]) and (lower<self.parent.sample_data.datasize[1])

            pixel = []
            for col in range(left,right+1):
                for row in range(upper,lower+1):
                    pixel.append( (row,col) )

        elif shape['type']=='Circle':
            center = shape['center']
            radius = shape['radius']
            # through all pixel within radius and see if they are in the circle
            pixel = []
            radceilint = np.int(np.ceil( radius ))
            for row in range( center[0]-radceilint, center[0]+radceilint+3 ):
                for col in range( center[1]-radceilint, center[1]+radceilint+3 ):
                    if (row-center[0])**2+(col-center[1])**2 <= radius**2:
                        pixel.append( (row,col) )

        else:
            raise ValueError('Spot class, shape key %s not understood.' % shape['type'])
        
        self.pixel = pixel

    def calculate_spot_std( self, which_data, int_type ):
        if which_data=='sample':
            rawdata = self.parent.sample_data.rawdata
        elif which_data=='blank':
            rawdata = self.parent.blank_data.rawdata
        elif which_data=='exspot':
            rawdata = self.parent.exspot_data.rawdata
        else:
            raise hell

        data = np.array( [rawdata[:,p[0],p[1]] for p in self.pixel] ).T
        # now we have a 2d array, one row per frame, one column per pixel

        std = np.std( data, axis=1 )
        # print 'spot std is:'
        # print std

        return std


    def calculate_spot_intensity( self, which_data, int_type ):
        if which_data=='sample':
            rawdata = self.parent.sample_data.rawdata
        elif which_data=='blank':
            rawdata = self.parent.blank_data.rawdata
        elif which_data=='exspot':
            rawdata = self.parent.exspot_data.rawdata
        else:
            raise hell

        data = np.array( [rawdata[:,p[0],p[1]] for p in self.pixel] ).T
        # now we have a 2d array, one row per frame, one column per pixel

        # for i,p in enumerate(self.pixel):
        #     print "(%f,%f) %f" % (p[0],p[1],data[0,i])

        if int_type=='mean':
            I = np.sum( data, axis=1 )/float(len(self.pixel))
            # I = np.array( [np.sum(d).astype(np.float) for d in sample] )
            # I /= self.width*self.height
        elif int_type=='max':
            I = np.max( data, axis=1 )
#            I = np.array( [np.max(d).astype(np.float) for d in sample] )
        elif int_type=='min':
            I = np.min( data, axis=1 )
#            I = np.array( [np.min(d).astype(np.float) for d in sample] )
        elif int_type=='gaussian':
            self.gaussianfit( data )       ## broken, for testing only
        else:
            raise ValueError("Spot class: bad int_type='%s'" % (int_type))
        return I


    def gaussianfit(self, rawdata):        
        x = np.arange( rawdata.shape[1] )
        y = np.arange( rawdata.shape[2] )
        X,Y = np.meshgrid( x,y )
        
        def gaussian2d( p, X, Y ):
            return p[0]*np.exp( -.5*((X-p[1])**2-(Y-p[2])**2)/p[3]**2 ) + p[4]

        def gaussian2d_err( p, X, Y, data ):
            err = np.sum( (gaussian2d(p,X,Y)-data)**2 )
            # make it really bad if parameters are outside reasonable bounds
#            if p[1]<0 or p[2]<0 or p[1]>data.shape[1] or p[2]>data.shape[0]

        initpars = [1.0, 1.0, 1.0, 1.0, 0.0]
        pars = []
        for i,data in enumerate(rawdata):
            pars.append( so.fmin( gaussian2d_err, initpars, (X,Y,data) ) )

        print pars

        

    def __str__(self):        
        return "spot: "+str(self.shape)

        # return "Spot object: int_type=%s\tlowerleft=[%d,%d]\twidth=%d\theight=%d\tlabel=%s" % \
        #     (self.intensity_type, self.coords[0], self.coords[1], \
        #          self.coords[2]-self.coords[0]+1, self.coords[3]-self.coords[1]+1, self.label)

#     def store_property_in_image(self, image, prop):
#         # record the spot in the mean intensity image
#         for p in self.pixel:
#             image[ p[0], p[1] ] = getattr(self,prop)
# #        image[ self.coords[1]:self.coords[3]+1, self.coords[0]:self.coords[2]+1 ] = getattr(self,prop)


    def export_averagematrix(self,filename):
        # np.save(filename, self.averagematrix)
        np.save(filename, self.recover_average_portrait_matrix() )


    def recover_average_portrait_matrix( self ):
        mycos = lambda a, ph, I, M: I*( 1+M*( np.cos(2*(a-ph)) ) )
        
        pic = np.zeros( (self.parent.emission_angles_grid.size, \
                             self.parent.excitation_angles_grid.size) )

        for pi in range(self.parent.Nportraits):
            for exi in range( self.parent.excitation_angles_grid.size ):
                pic[:,exi] += mycos( self.parent.emission_angles_grid, \
                                        self.verticalfitparams[pi,exi,0], \
                                        self.verticalfitparams[pi,exi,1], \
                                        self.verticalfitparams[pi,exi,2] )
        pic /= float(self.parent.Nportraits)
        return pic


    def portrait_residuals( self, iportrait ):
        mycos = lambda a, ph, I, M: I*( 1+M*( np.cos(2*(a-ph)) ) )

        # indices of the angles which were measured in the portrait
        #iuexa = self.parent.uexa_portrait_indices
        #iuema = self.parent.uema_portrait_indices
      
        experimental_portrait = np.zeros( (self.parent.Nlines, len(self.parent.unique_exangles)) )
        experimental_exangles = np.zeros( (self.parent.Nlines, len(self.parent.unique_exangles)) )
        experimental_emangles = np.zeros( (self.parent.Nlines, len(self.parent.unique_exangles)) )

        line_fitted_portrait  = np.zeros( (self.parent.Nlines, len(self.parent.unique_exangles)) )
        line_fitted_exangles  = np.zeros( (self.parent.Nlines, len(self.parent.unique_exangles)) )
        line_fitted_emangles  = np.zeros( (self.parent.Nlines, len(self.parent.unique_exangles)) )
        
        full_fitted_portrait  = np.zeros( (self.parent.Nlines, len(self.parent.unique_exangles)) )
        full_fitted_exangles  = np.zeros( (self.parent.Nlines, len(self.parent.unique_exangles)) )
        full_fitted_emangles  = np.zeros( (self.parent.Nlines, len(self.parent.unique_exangles)) )
        # go through all lines
        for iline in range(self.parent.Nlines):
            # get correct angles here
            exangles,emangle = self.parent.retrieve_angles( iportrait, iline )
            # arrange lines of spot intensities into rows in this array
            experimental_portrait[iline,:] = self.retrieve_intensity( iportrait, iline )
            experimental_exangles[iline,:] = exangles
            experimental_emangles[iline,:] = np.ones_like(exangles)*emangle
            # arrange fitted lines in the same way
            line_fitted_portrait[iline,:]  = self.retrieve_line_fit( iportrait, iline, exangles )
            line_fitted_exangles[iline,:]  = exangles
            line_fitted_emangles[iline,:] = np.ones_like(exangles)*emangle
            # now the fully fitted portraits (ie vertical fits)
            for iexa,exa in enumerate(exangles):
                full_fitted_exangles[ iline,iexa ] = exa
                full_fitted_emangles[ iline,iexa ] = emangle
                # what's the index of this exangle in the portrait?
                piexa = np.argmin( np.abs(exa-self.parent.excitation_angles_grid) )
                assert np.min( np.abs(exa-self.parent.excitation_angles_grid) ) < 10*np.finfo(np.float).eps
                # got it? ok, then use that vertical fit 
                full_fitted_portrait[iline,iexa] = mycos( emangle, \
                                        self.verticalfitparams[iportrait,piexa,0], \
                                        self.verticalfitparams[iportrait,piexa,1], \
                                        self.verticalfitparams[iportrait,piexa,2] )

        assert np.all(experimental_exangles == line_fitted_exangles) 
        assert np.all(experimental_exangles == full_fitted_exangles)
        assert np.all(experimental_emangles == line_fitted_emangles)
        assert np.all(experimental_emangles == full_fitted_emangles)

        residuals_line = experimental_portrait - line_fitted_portrait
        residuals_full = experimental_portrait - full_fitted_portrait

        # print experimental_portrait
        # print line_fitted_portrait
        # print full_fitted_portrait
        # print residuals_line
        # print residuals_full
        print 'residual of line fits    : ',np.linalg.norm(residuals_line)
        print 'residual of vertical fits: ',np.linalg.norm(residuals_full)
        self.residual_line = np.linalg.norm(residuals_line)
        self.residual_full = np.linalg.norm(residuals_full)

        return self.residual_full

    def retrieve_line_fit( self, iportrait, iline, angles ):
        I0 = self.linefitparams[iportrait, iline, 1]
        M0 = self.linefitparams[iportrait, iline, 2]
        ph = self.linefitparams[iportrait, iline, 0]
        return  I0 *( 1+M0*np.cos( 2*(angles-ph) ) )


    def retrieve_intensity( self, iportrait, iline ):
        pstart = self.parent.portrait_indices[ iportrait ]
        pstop  = self.parent.portrait_indices[ iportrait+1 ]
        # lstart = self.parent.line_edges[ iline ]
        # lstop  = self.parent.line_edges[ iline+1 ]
        # print 'portrait %d line %d:' % (iportrait, iline)
        # print self.intensity[ pstart:pstop ][ lstart:lstop ]
        # print self.intensity
        # print '==============================================='
#        return self.intensity[ pstart:pstop ][ lstart:lstop ]
        return self.intensity[ pstart:pstop ][ self.parent.line_indices[iportrait][iline] ]


    def values_for_ETruler( self, newdatalength=2048, emainc=10*np.exp(1), showfft=False ):

        def look_up_cosine( self, exa, emas, Ndigits ):
            mycos = lambda a, ph, I, M: I*( 1+M*( np.cos(2*(a-ph)) ) )

            which_portrait = 0
            iexa = np.argmin( np.abs(exa-self.parent.excitation_angles_grid) )
            assert np.abs(exa-self.parent.excitation_angles_grid[iexa]) < 1e-5

            return mycos( emas, \
                              self.verticalfitparams[which_portrait,iexa,0], \
                              self.verticalfitparams[which_portrait,iexa,1], \
                              self.verticalfitparams[which_portrait,iexa,2] )
        

        # uexa   = np.array([0.0, 30.0, 60.0, 90.0, 120.0, 150.0])/180.0*np.pi
        uexa   = self.parent.unique_exangles
        exdiff = np.pi/len(uexa)

        slope = emainc/180.0        # go up emainc degree in emission for 180 deg in excitation
        exa = np.linspace(0, newdatalength-1, newdatalength)*exdiff
        ema = slope * exa

        # now round and wrap back to the [0,180)[0,180) portrait
        Ndigits = 6
        exa = np.round( np.mod( exa, np.pi ), Ndigits )
        ema = np.round( np.mod( ema, np.pi ), Ndigits )

        # import matplotlib.pyplot as plt
        # plt.plot( exa, ema, '.', alpha=.5 )
        # plt.show()

#        print np.unique( exa ).size
#        print np.unique( ema ).size

        intensities = np.zeros_like(ema)
        for i in range(len(uexa)):
            intensities[ i::len(uexa) ] = look_up_cosine( self, uexa[i], ema[i::len(uexa)], Ndigits )
            # print ema[i::len(uexa)].size
            # print ema[i::len(uexa)][-1]

        Fint  = np.fft.fft(intensities)
        Pint  = np.real( Fint*np.conj(Fint) )/float(exa.size)
        nPint = Pint[1:newdatalength/2]/np.sum(Pint[1:newdatalength/2])

        kappa = .5 * np.arccos( .5*(3*self.M_ex-1) )
        alpha = np.array( [-kappa, 0, kappa] )

        modelNoET = np.zeros_like( exa )
        for n in range(3):
            modelNoET[:] += np.cos(exa-alpha[n])**2 * np.cos(ema-alpha[n])**2
        modelNoET /= 3.0
        FmodelNoET  = np.fft.fft( modelNoET )
        PmodelNoET  = np.real( FmodelNoET*np.conj(FmodelNoET) )/float(exa.size)
        nPmodelNoET = PmodelNoET[1:newdatalength/2]/np.sum(PmodelNoET[1:newdatalength/2])

        # We are still awesome, and this time we think we know why:
        #
        # Both frequencies depend on the data length (because we never translate into actual fft freqs).
        # How many steps per 180deg do we do? --> That is our sampling frequency in excitation and emission.
        # The number of steps can be calculated from:
        #  - the number of angles in excitation (which goes from 0..180, excluding 180), 
        #    represented by the grid size
        #  - the slope, as given by emainc/180, which defines the 'vertical' sampling frequency in the portrait

        # horizontal freq (that is, sampling in excitation):
        number_of_samples_in_ex = float(uexa.size)
        hf = np.round( newdatalength/number_of_samples_in_ex )
#        print newdatalength
#        print number_of_samples_in_ex
#        print hf
        # vertical freq (in emission)
        exstep_deg = 180/number_of_samples_in_ex
        emstep_deg = slope*exstep_deg
        number_of_samples_in_em = 180/emstep_deg
        vf = np.round( newdatalength/number_of_samples_in_em )

        # satellites
        satdiff = np.abs( hf-vf )
        satsum  = np.abs( hf+vf )

        # everything has been shifted by one because we excluded the first element when normalizing
        hf -= 1
        vf -= 1
        satdiff -= 1
        satsum  -= 1

        # freq-difference for the summation over peaks: 1/3rd of the difference between the largest two
        # frequencies. We sort the frequencies first because we do not know if hf>vf (depends on chosen
        # slope (emainc)). 
        sortedfreqs = np.sort(np.array([hf,vf,satdiff,satsum]))
        df = np.diff( sortedfreqs )[-1]
        df = np.round( df/3.0 )

        # summation bounds:
        sb = np.array( [[f-df,f+df] for f in sortedfreqs] )
        sb[sb<0] = 0

        # peak intensities:
        peakints = np.array( [ np.sum( nPint[int(sb[i,0]):int(sb[i,1])] ) for i in range(sortedfreqs.size) ] )

        # peak intensities for the no-ET model:
        modelNoETpeakints = np.array( [ np.sum( nPmodelNoET[int(sb[i,0]):int(sb[i,1])] ) \
                                            for i in range(sortedfreqs.size) ] )

        # if peak sums make no sense, then say something and set result to NaN
        if (np.abs(np.sum(peakints)-1) > 0.08):
            print "ET ruler: peakints sum is bad..."
            self.ET_ruler = np.nan
            return
        if (np.abs(np.sum(modelNoETpeakints)-1) > 0.08):
            print "ET ruler: MODELpeakints sum is bad..."
            self.ET_ruler = np.nan
            return

        crossdiff      = np.abs( peakints[1]-peakints[3] )
        modelcrossdiff = np.abs( modelNoETpeakints[1]-modelNoETpeakints[3] )

        ruler = 1-(crossdiff/modelcrossdiff)

        if (ruler < -.2) or (ruler > 1.2):
            print "\t\tRuler has gone bonkers (ruler=%f)." % (ruler)
            print "\t\tWill continue anyways and set ruler to NaN."
            ruler = np.nan
        # else:
        #     if ruler < 0:
        #         ruler = 0
        #     if ruler > 1:
        #         ruler = 1
            
        self.ET_ruler = ruler
        self.parent.store_property_in_image( self, 'ET_ruler_image', 'ET_ruler' )
#        print 'ruler= ',ruler

        # Here we've reached our peaks in life. Only downhill from here...

        if showfft:
            import matplotlib.pyplot as plt
            plt.plot( nPint )
            plt.plot( nPmodelNoET, 'r' )
            plt.show()

        return exa,ema,intensities
