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
        np.save(filename,self.averagematrix)

    def recover_average_portrait_matrix(self):
        pic = self.portraits[0].recover_portrait_matrix()
        n = 1
        while n < len(self.portraits):
            pic += self.portraits[n].recover_portrait_matrix()
            n += 1
        pic /= n
        return pic


    def collect_and_assign(self):
        # grab angles
        if self.parent.which_setup=='cool new setup':
            exangles = self.parent.motors.excitation_angles
            emangles = self.parent.motors.emission_angles
        else:
            exptime = self.parent.camera_data.exposuretime
            exangles = np.array( [self.parent.excitation_motor.angle(t,exposuretime=exptime) \
                                      for t in self.parent.timeaxis] )
            emangles = np.array( [self.parent.emission_motor.angle(t,exposuretime=exptime) \
                                      for t in self.parent.timeaxis] )
        # truth value array for frame validity
        validframes = emangles != -1
        self.parent.Nvalidframes = np.sum(validframes)

        # pick angles of valid frames and round them
        emangles = np.round(emangles[validframes], decimals=2)
        exangles = np.round(exangles[validframes], decimals=2)

        # indices of valid frames
        trueindices = validframes.nonzero()[0]

        # store valid angles and frames temporarily
        ex = []
        em = []
        I  = []
#        for i in range( self.timeaxis[validframes].size ):
        for i in range( np.sum(validframes) ):
            ex.append( exangles[i] )
            em.append( emangles[i] )
            I.append( self.intensity[trueindices[i]] )
        ex = np.array(ex)
        em = np.array(em)
        I  = np.array(I)
                
        # grab portrait indices (generated by Movie.startstop())
        pind = self.parent.portrait_indices
        Nportraits = pind.shape[0]

        # make a list of portraits, using stored angles and frames
        portraitlist = []
        for n in range(Nportraits):
            portraitlist.append( Portrait( ex[pind[n,0]:pind[n,1]+1], em[pind[n,0]:pind[n,1]+1], \
                                               I[pind[n,0]:pind[n,1]+1], self ) )
        # store list in spot object
        self.portraits = portraitlist

        # the spot.intensity is no longer needed
        del(self.intensity)

        # we're done here



    def get_portrait_data(self):
        pind = self.parent.portrait_indices
        Nportraits = pind.shape[0]

        # create empty list
        portraitlist = []
        # go through all portraits
        for n in range(Nportraits):
            # grab nth portrait
            d = self.parent.data[pind[n,0]:pind[n,1]+1, :]
            if self.parent.datamode=='truedata':
                # use only valid rows
                d = d[d[:,1]==1]
                # from the remaining (valid) rows, grab
                # columns ex, em and Int
                ex = d[:, 2]
                em = d[:, 3]
                I  = d[:, 4+si]
            elif self.parent.datamode=='validdata':            
                # simply grab columns ex, em and Int
                ex = d[:, 1]
                em = d[:, 2]
                I  = d[:, 3+si]
            # store portrait object in list
            portraitlist.append( Portrait( ex, em, I, spot ) )

        # store list in spot object
        self.portraits = portraitlist


    def check_if_valid(self,SNR):
        # do we actually have the background std
        bgstd = 0
        if hasattr( self.parent, 'bg_spot' ):
            bgstd = self.parent.bg_spot.std
        else:
            print "Dude --- no background spot defined, therefore no standard deviation. Will treat all spots as valid (i.e. as having sufficient intensity)."

        self.SNR = self.mean_intensity/bgstd
        if self.SNR > SNR:
            self.isvalid = True
        else:
            self.isvalid = False

        
    def cos_fit(self):
        self.residual = 0

        # we assume that the number of portraits and lines is the same 
        # for all spots (can't think of a reason why that shouldn't be the case).
        Nportraits = len(self.portraits)
        Nlines     = len(self.portraits[0].lines )

        # for each portrait ---- outermost loop, we do portraits in series
        for pi in range(Nportraits):
            # part I, 'horizontal fitting' of the lines of constant emission angles

            # for each line ---- we do lines in series, __but all spots in parallel__:
            for li in range(Nlines):

                # get excitation angle array (same for all spots!)
                exangles    = self.portraits[pi].lines[li].exangles

                # create list of intensity arrays (one array for each spot)
                intensities = self.portraits[pi].lines[li].intensities

                # turn into numpy array and make sure it's a column
                intensities = np.array( intensities ).reshape((exangles.size,1))

                exa = exangles.copy()

                phase, I0, M, resi, fit, rawfitpars, mm = self.parent.cos_fitter( exa, intensities, \
                                                                               self.parent.Nphases_for_cos_fitter ) 
                # write cosine parameters into line object
                self.portraits[pi].lines[li].set_fit_params( phase, I0, M, resi )
                # and add residual
                self.residual += resi


            # part II, 'vertical fitting' --- we do each spot by itself, but 
            # fit all verticals in parallel

            # collect list of unique emission angles            
            emangles = [l.emangle for l in self.portraits[pi].lines]
            # turn into array, transpose and squeeze
            emangles = np.squeeze(np.array( emangles ).T)

            # evaluate cosine-fit at these em_angles, on a grid of ex_angles:
            fitintensities = np.array( [l.cosValue( self.parent.excitation_angles_grid ) \
                                            for l in self.portraits[pi].lines] )

            phase, I0, M, resi, fit, rawfitpars, mm = self.parent.cos_fitter( emangles, fitintensities, \
                                                                           self.parent.Nphases_for_cos_fitter )
            # store vertical fit params
            self.portraits[pi].vertical_fit_params = [ phase, I0, M, resi, mm ]

    def findModDepths(self):
        sam = self.recover_average_portrait_matrix()
        self.proj_ex = np.mean( sam, axis=0 ).reshape((sam.shape[1],1))
        self.proj_em = np.mean( sam, axis=1 ).reshape((sam.shape[0],1))

        # fitting
        ph_ex, I_ex, M_ex, r_ex, fit_ex, rawfitpars_ex, mm = self.parent.cos_fitter( self.parent.excitation_angles_grid, self.proj_ex, self.parent.Nphases_for_cos_fitter )
        ph_em, I_em, M_em, r_em, fit_em, rawfitpars_em, mm = self.parent.cos_fitter( self.parent.emission_angles_grid, self.proj_em, self.parent.Nphases_for_cos_fitter )

        # assignment
        LS = ph_ex - ph_em
        if LS > np.pi/2:  LS -= np.pi
        if LS < -np.pi/2: LS += np.pi
        self.phase_ex = ph_ex[0]
        self.M_ex     = M_ex[0]
        self.phase_em = ph_em[0]
        self.M_em     = M_em[0]
        self.LS       = LS[0]
        # store in coverage maps
        self.parent.M_ex_image[ self.coords[1]:self.coords[3]+1, self.coords[0]:self.coords[2]+1 ] = self.M_ex
        self.parent.M_em_image[ self.coords[1]:self.coords[3]+1, self.coords[0]:self.coords[2]+1 ] = self.M_em
        self.parent.phase_ex_image[ self.coords[1]:self.coords[3]+1, self.coords[0]:self.coords[2]+1 ] = self.phase_ex
        self.parent.phase_em_image[ self.coords[1]:self.coords[3]+1, self.coords[0]:self.coords[2]+1 ] = self.phase_em
        self.parent.LS_image[ self.coords[1]:self.coords[3]+1, self.coords[0]:self.coords[2]+1 ] = self.LS
