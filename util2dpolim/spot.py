import numpy as np

class Spot:
    def __init__(self, rawdata, coords, bg, int_type, label, parent, is_bg_spot=False, \
                     blankdata=False):
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
        """

        self.coords = coords    #[left, bottom, right, top]
        self.label  = label
        self.width  = coords[2] -coords[0] +1
        self.height = coords[3] -coords[1] +1
        self.parent = parent

        # Fit parameter storage
        # dimensions: Nportraits x Nlines x 4
        # Last dimension encodes [phase, I0, M, residual]
        self.linefitparams     = np.nan*np.ones( (self.parent.Nportraits, self.parent.Nlines, 4) )
        self.verticalfitparams = np.nan*np.ones( (self.parent.Nportraits, \
                                                      self.parent.excitation_angles_grid.size, 4) )

        # work out frame-dependent intensities for that spot
        # - mean:
        if int_type=='mean':
            I  = np.sum( np.sum( \
                    rawdata[:, coords[1]:coords[3]+1, coords[0]:coords[2]+1 ], \
                        axis=2), axis=1 ).astype( np.float )
            I /= self.width*self.height
            # work out blank signal if present
            if blankdata:
                Iblank  = np.sum( self.blank_image[ coords[1]:coords[3]+1, coords[0]:coords[2]+1 ] )
                Iblank /= self.width*self.height

        elif int_type=='max':
            # - maximum:
            I = np.max( np.max( \
                    rawdata[:, coords[1]:coords[3]+1, coords[0]:coords[2]+1 ], \
                        axis=2), axis=1 ).astype( np.float )
            # work out blank signal if present
            if blankdata:
                Iblank  = np.max( self.blank_image[ coords[1]:coords[3]+1, coords[0]:coords[2]+1 ] )
                Iblank /= self.width*self.height

        elif int_type=='min':
            # - minimum:
            I = np.min( np.min( \
                    self.camera_data.rawdata[:, coords[1]:coords[3]+1, coords[0]:coords[2]+1 ], \
                        axis=2), axis=1 ).astype( np.float )
            # work out blank signal if present
            if blankdata:
                Iblank  = np.min( self.blank_image[ coords[1]:coords[3]+1, coords[0]:coords[2]+1 ] )
                Iblank /= self.width*self.height

        else:
            raise ValueError("Spot __init__ did not understand int_type='%s' (should be mean|max|min)" % (int_type))

        # special: take standard deviation if this is the background spot
        if is_bg_spot:
            self.std  = np.std( rawdata[:, coords[1]:coords[3]+1, coords[0]:coords[2]+1 ] )

        # remove background
        I -= bg
        # remove blank
        if blankdata:
            I -= Iblank

        self.intensity_type = int_type
        self.intensity      = I
        self.mean_intensity = np.mean(I)
        self.bg_correction  = bg
        if not is_bg_spot:
            self.parent.mean_intensity_image[ \
                self.coords[1]:self.coords[3]+1, self.coords[0]:self.coords[2]+1 \
                    ] = self.mean_intensity


    def __str__(self):
        return "Spot object: int_type=%s\tlowerleft=[%d,%d]\twidth=%d\theight=%d\tlabel=%s" % \
            (self.intensity_type, self.coords[0], self.coords[1], \
                 self.coords[2]-self.coords[0]+1, self.coords[3]-self.coords[1]+1, self.label)

    def store_property_in_image(self, image, prop):
        image[ self.coords[1]:self.coords[3]+1, self.coords[0]:self.coords[2]+1 ] = getattr(self,prop)

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
        pic /= self.parent.Nportraits
        return pic


    def retrieve_line_fit( self, iportrait, iline, angles ):
        I0 = self.linefitparams[iportrait, iline, 1]
        M0 = self.linefitparams[iportrait, iline, 2]
        ph = self.linefitparams[iportrait, iline, 0]
        return  I0 *( 1+M0*np.cos( 2*(angles-ph) ) )


    def retrieve_intensity( self, iportrait, iline ):
        pstart = self.parent.portrait_indices[ iportrait, 0 ]
        pstop  = self.parent.portrait_indices[ iportrait, 1 ]+1
        lstart = self.parent.line_edges[ iline ]
        lstop  = self.parent.line_edges[ iline+1 ]
        return self.intensity[ pstart:pstop ][ lstart:lstop ]
