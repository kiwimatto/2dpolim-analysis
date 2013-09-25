from cameradata import CameraData
from fitting import CosineFitter_new
from spot import Spot
from portrait import Portrait
from motors import *
import os
import scipy.optimize as so

class Movie:
    def __init__( self, datadir, basename, phase_offset_in_deg=np.nan ):

        self.cos_fitter = CosineFitter_new
        self.phase_offset_in_deg = phase_offset_in_deg

        self.data_directory = datadir
        self.data_basename  = basename

        self.read_in_EVERYTHING()

        self.spots = []
        self.initContrastImages()

        # fix the excitation and emission angle grids for now
        self.excitation_angles_grid = np.linspace(0,np.pi,91)
        self.emission_angles_grid = np.linspace(0,np.pi,91)
        self.Nphases_for_cos_fitter = 91
        self.precomputed_cosines = [ np.cos(2*(self.emission_angles_grid-ph)) \
                                         for ph in np.linspace(0,np.pi/2,self.Nphases_for_cos_fitter) ]

    def read_in_EVERYTHING(self):
        # change to the data directory
        os.chdir( self.data_directory )
        
        ###### look if we can find the data file, and import it ######
        print 'Looking for file: %s.[spe|SPE]' % (self.data_directory+self.data_basename)

        if os.path.exists( self.data_basename+'.SPE' ):
            self.spefilename = self.data_basename+'.SPE'
        elif os.path.exists( self.data_basename+'.spe' ):
            self.spefilename = self.data_basename+'.spe'
        else:
            print "Couldn't find data SPE file! Bombing out..."
            raise SystemExit
        
        self.camera_data    = CameraData( self.spefilename, compute_frame_average=True )
        self.timeaxis       = self.camera_data.timestamps
        print 'Imported data file %s' % self.spefilename


        ###### look for motor files ######
        print 'Looking for motor file(s)...'
        
        got_motors = False
        for file in os.listdir("."):
            if file=="MS-"+self.data_basename+'.txt':
                print '\t found motor file %s' % file
                self.motorfile = file
                got_motors = True
                break

        # now import it depending on which setup it is
        if got_motors:
            ### Supplying the phase offset here means it overrides the value 
            ### which is read from the header, unless it is NaN (the default). 
            ### This is only useful for AM measurements.
            self.motors = BothMotorsWithHeader( self.motorfile, self.phase_offset_in_deg )
            self.exangles = self.motors.excitation_angles
            self.emangles = self.motors.emission_angles
        else:
            raise IOError("Couldn't find motor files.")
       
        ###### look for blank sample ######
        print 'Looking for blank...',
        self.blanks = []
        for file in os.listdir("."):
            if file.startswith("blank-") and (file.endswith(".spe") or file.endswith(".SPE")):
                print '\t found file %s' % file
                b = CameraData( file, compute_frame_average=True )
                self.blanks.append( b )


    def initContrastImages(self):
        self.spot_coverage_image  = np.ones( (self.camera_data.datasize[1],self.camera_data.datasize[2]) )*np.nan
        self.mean_intensity_image = self.spot_coverage_image.copy()
        self.SNR_image            = self.spot_coverage_image.copy()
        self.M_ex_image           = self.spot_coverage_image.copy()
        self.M_em_image           = self.spot_coverage_image.copy()
        self.phase_ex_image       = self.spot_coverage_image.copy()
        self.phase_em_image       = self.spot_coverage_image.copy()
        self.LS_image             = self.spot_coverage_image.copy()
        self.r_image              = self.spot_coverage_image.copy()
        self.ET_ruler_image       = self.spot_coverage_image.copy()
        self.ET_model_md_fu_image = self.spot_coverage_image.copy()
        self.ET_model_th_fu_image = self.spot_coverage_image.copy()
        self.ET_model_gr_image    = self.spot_coverage_image.copy()
        self.ET_model_et_image    = self.spot_coverage_image.copy()


    def define_background_spot( self, coords, intensity_type='mean' ):
        # create new spot object
        s = Spot( self.camera_data.rawdata, coords, bg=0, int_type=intensity_type, \
                      label='background area', is_bg_spot=True, parent=self )
        self.bg_spot = s        

        # if blank data is present, automatically work out its bg, correct for it,
        # and compute average blank image
        if len(self.blanks)>0:
            # init a blank image
            self.blank_image = np.zeros_like(self.blanks[0])
            # go through all blanks
            for b in self.blanks:
                # get that background spot
                s = Spot( b.resize((1,b.shape[0],b.shape[1])), coords, bg=0, int_type=intensity_type, \
                              label='background area', is_bg_spot=True, parent=self )
                # add bg-corrected blank to blank_image
                self.blank_image += b-s.I
            # divide blank image by number of blanks
            self.blank_image /= len(self.blanks)

        # record background spot in spot coverage image
        self.spot_coverage_image[ s.coords[1]:s.coords[3]+1, s.coords[0]:s.coords[2]+1 ] = -1
        

    def define_spot( self, coords, intensity_type='mean', label=None ):
        """Defines a new spot object and adds it to the list.
        FIXME: make sure coordinate definitions do not exceed frame size
        TO-DO: how about being able to specify center+radius ?
        """

        if hasattr( self, 'bg_spot' ):
            bgself = self.bg_spot.intensity
        else:
            bgself = 0
                
        if hasattr( self, 'blank_image' ):
            bgblank=True
        else:
            bgblank=False

        # create new spot object
        s = Spot( self.camera_data.rawdata, coords, bg=bgself, int_type='mean', \
                      label=label, parent=self, blankdata=bgblank )
        # append spot object to spots list
        self.spots.append( s )

        self.spot_coverage_image[ s.coords[1]:s.coords[3]+1, s.coords[0]:s.coords[2]+1 ] = 1
        self.mean_intensity_image[ s.coords[1]:s.coords[3]+1, s.coords[0]:s.coords[2]+1 ] = s.mean_intensity


    def correct_excitation_intensities( self ):
        for i,s in enumerate(self.spots):
            s.intensity /= self.motors.sample_plane_intensities


    def correct_emission_intensities( self, corrM, corrphase ):
        # correction function
        corrfun = lambda angle: (1+corrM*np.cos(2*(angle + corrphase) ))/(1+corrM)

        # assemble all corrections in a vector
        corrections = np.nan*np.ones((self.emangles.size,))
        for i,emangle in enumerate(self.emangles):
            corrections[i] = corrfun( emangle )

        # now go through all spots and apply that vector
        for i,s in enumerate(self.spots):
            s.intensity /= corrections

    def startstop( self ):
        """
        This function determines the indices at which portraits start and end.
        There's a fair bit of hackery in here, and to understand what is 
        happening one really needs to look at the raw data, frame indices, valid
        frames (no shutter), etc...  This will be need to be documented much more
        thoroughly in the future, but for now: Handle with care.
        """

        emangles_rounded = np.round( self.motors.emission_angles, decimals=2 )

        number_of_lines = np.unique( emangles_rounded ).size

        # edge 'detection' via diff
        d = np.diff( emangles_rounded )
        d[0] = 1
        d[d!=0] = 1
        d = np.concatenate( (d, np.array([1])) )
        edges    = d.nonzero()[0]
        
        # integer division to find out how many complete portraits we have
        Nportraits = np.diff(edges).size / number_of_lines 

        indices = np.zeros( (Nportraits,2), dtype=np.int )
        for i in range(Nportraits):
            indices[i,0] = edges[i*number_of_lines]+1
            indices[i,1] = edges[(i+1)*number_of_lines]

        # These indices will discard the first frame, despite it being a valid frame, so:
        indices[0,0] = 0

        self.portrait_indices = indices
        # this completes the determination of portrait indices

        self.line_edges = edges[:number_of_lines+1]
        self.line_edges[1:] += 1

        self.Nportraits = self.portrait_indices.shape[0]
        self.Nlines     = len( self.line_edges )-1
        

    def write_data( filename, header=False ):
        """Helper-function which takes the output of collect_data() and
        writes it (possibly with header) into a file.
        """
        fhandle = open( filename, 'wt' )
        if self.datamode=='truedata':
            if header:
                fhandle.write("FrameNumber\t Validity\t excitation\t emission\t Intensities....\n")
            np.savetxt( fhandle, self.data, '%f\t' )
        elif self.datamode=='validdata':
            if header:
                fhandle.write("FrameNumber\t excitation\t emission\t Intensities....\n")
            np.savetxt( fhandle, self.data, '%f\t' )
        else:
            fhandle.close()
            raise hell

        fhandle.close()


    def compute_modulation_in_emission( self,portrait ):
        # part II, 'vertical fitting' --- we do each spot by itself, but 
        # fit all verticals in parallel

        # collect list of unique emission angles (same for all spots!)                    
        emangles = [l.emangle for l in portrait.lines]
        # turn into array, transpose and squeeze
        emangles = np.squeeze(np.array( emangles ).T)

        # evaluate cosine-fit at these em_angles, on a grid of ex_angles:
        fitintensities = np.array([l.cosValue( self.excitation_angles_grid ) for l in portrait.lines])

        phase, I0, M, resi, fit, rawfitpars, mm = self.cos_fitter( emangles, fitintensities, \
                                                                       self.Nphases_for_cos_fitter )         
        # phasor addition!
        proj_em = np.real( rawfitpars[0,:] * np.exp( 1j*rawfitpars[1,:] ) \
                               * np.exp( 1j*self.emission_angles_grid ) )

        phase, I0, M, resi, fit, rawfitpars, mm = self.cos_fitter( self.emission_angles_grid, proj_em, \
                                                                       self.Nphases_for_cos_fitter ) 
        portrait.phase_em = phase[0]
        portrait.M_em = M
        

    def fit_all_portraits_spot_parallel_selective( self, myspots=None ):

        if myspots==None:
            myspots = range(len(self.validspots))

        # init average portrait matrices, so that we can write to them without
        # having to store a matrix for each portrait
        for si in myspots:
            self.validspots[si].residual = 0

        # for each portrait ---- outermost loop, we do portraits in series
        for pi in range(self.Nportraits):

            # start and stop indices for this portrait
            pstart = self.portrait_indices[ pi, 0 ]
            pstop  = self.portrait_indices[ pi, 1 ]+1
            
            # part I, 'horizontal fitting' of the lines of constant emission angles

            # for each line ---- we do lines in series, __but all spots in parallel__:
            for li in range(self.Nlines):

                # start and stop indices for this line
                lstart = self.line_edges[ li ]
                lstop  = self.line_edges[ li+1 ]

                # get excitation angle array (same for all spots!)
                exangles = self.exangles[ pstart:pstop ][ lstart:lstop ]

                # create list of intensity arrays (one array for each spot)
                intensities = [self.validspots[si].retrieve_intensity(iportrait=pi, iline=li) for si in myspots]

                # turn into numpy array and transpose
                intensities = np.array( intensities ).T

                exa = exangles.copy()

                phase, I0, M, resi, fit, rawfitpars, mm = self.cos_fitter( exa, intensities, \
                                                                               self.Nphases_for_cos_fitter ) 

                # write cosine parameters into line object
                for sii,si in enumerate(myspots):
                    self.validspots[si].linefitparams[pi,li,0] = phase[sii]
                    self.validspots[si].linefitparams[pi,li,1] = I0[sii]
                    self.validspots[si].linefitparams[pi,li,2] = M[sii]
                    self.validspots[si].linefitparams[pi,li,3] = resi[sii]

            # gather residuals for this protrait
            for si in myspots:
                self.validspots[si].residual = np.sum( self.validspots[si].linefitparams[:,:,3] )

            # part II, 'vertical fitting' --- we do each spot by itself, but 
            # fit all verticals in parallel

            # collect list of unique emission angles (same for all spots!)                    
            emangles = self.emangles[pstart:pstop][self.line_edges[:-1]]
#            emangles = [l.emangle for l in self.validspots[0].portraits[pi].lines]
            # turn into array, transpose and squeeze
            emangles = np.squeeze(np.array( emangles ).T)

            # evaluate cosine-fit at these em_angles, on a grid of ex_angles:
            fitintensities = np.hstack( [ np.array( [ \
                            self.validspots[si].retrieve_line_fit( pi, li, self.excitation_angles_grid ) \
                                for li in range(self.Nlines) ] ) \
                                              for si in myspots ] )
            # fitintensities = [ np.array( \
            #         [ l.cosValue( self.excitation_angles_grid ) \
            #               for l in self.validspots[si].portraits[pi].lines ]) \
            #                        for si in myspots ]
            # fitintensities = np.hstack( fitintensities )

            phase, I0, M, resi, fit, rawfitpars, mm = self.cos_fitter( emangles, fitintensities, \
                                                                           self.Nphases_for_cos_fitter ) 
                
            # store vertical fit params
            phase = np.hsplit(phase, len(myspots))
            I0    = np.hsplit(I0, len(myspots))
            M     = np.hsplit(M, len(myspots))
            resi  = np.hsplit(resi, len(myspots))            
            for sii,si in enumerate(myspots):
                self.validspots[si].verticalfitparams[pi,:,0] = phase[sii]
                self.validspots[si].verticalfitparams[pi,:,1] = I0[sii]
                self.validspots[si].verticalfitparams[pi,:,2] = M[sii]
                self.validspots[si].verticalfitparams[pi,:,3] = resi[sii]


    def find_modulation_depths_and_phases_selective( self, myspots=None ):

        if myspots==None:
            myspots = range(len(self.validspots))

        # projection onto the excitation axis (ie over all emission angles),
        # for all spots, one per column
        proj_ex = []
        proj_em = []
        for si in myspots:
            sam  = self.validspots[si].recover_average_portrait_matrix()
            # sam /= np.max(sam)*255
            # s.sam = sam.astype( np.uint8 )
            self.validspots[si].sam = sam
            self.validspots[si].proj_ex = np.mean( sam, axis=0 )
            self.validspots[si].proj_em = np.mean( sam, axis=1 )
            proj_ex.append( self.validspots[si].proj_ex )
            proj_em.append( self.validspots[si].proj_em )
        proj_ex = np.array(proj_ex).T
        proj_em = np.array(proj_em).T

        # fitting
        ph_ex, I_ex, M_ex, r_ex, fit_ex, rawfitpars_ex, mm = \
            self.cos_fitter( self.excitation_angles_grid, proj_ex, self.Nphases_for_cos_fitter )
        ph_em, I_em, M_em, r_em, fit_em, rawfitpars_em, mm = \
            self.cos_fitter( self.emission_angles_grid, proj_em, self.Nphases_for_cos_fitter )

        # assignment
        LS = ph_ex - ph_em
        LS[LS >  np.pi/2] -= np.pi
        LS[LS < -np.pi/2] += np.pi
        for sii,si in enumerate(myspots):
            self.validspots[si].phase_ex = ph_ex[sii]
            self.validspots[si].M_ex     = M_ex[sii]
            self.validspots[si].phase_em = ph_em[sii]
            self.validspots[si].M_em     = M_em[sii]
            self.validspots[si].LS       = LS[sii]
            # store in coverage maps
            a = self.validspots[si].coords[1]
            b = self.validspots[si].coords[3]+1
            c = self.validspots[si].coords[0]
            d = self.validspots[si].coords[2]+1
            self.M_ex_image[ a:b, c:d ]     = self.validspots[si].M_ex
            self.M_em_image[ a:b, c:d ]     = self.validspots[si].M_em
            self.phase_ex_image[ a:b, c:d ] = self.validspots[si].phase_ex
            self.phase_em_image[ a:b, c:d ] = self.validspots[si].phase_em
            self.LS_image[ a:b, c:d ]       = self.validspots[si].LS        

            #### anisotropy ####
        
            if (self.validspots[si].M_ex < .15):
                iex = np.argmin( np.abs( self.excitation_angles_grid - np.pi/2 ) )
                iem = np.argmin( np.abs( self.emission_angles_grid - np.pi/2 ) )
                Ipara = self.validspots[si].sam[ iex, iem ]
                Iperp = self.validspots[si].sam[ iex, 0 ] ## assumes that the first index is close to 0deg in emission angle grid
            else:
                # find where the found phase matches the phase in the portrait matrix
                # portrait matrices are constructed from the angle grid, so look up which 
                # index corresponds to the found phase
                if self.validspots[si].phase_ex < 0:
                    iphex = np.argmin( np.abs(self.excitation_angles_grid - (np.pi+self.validspots[si].phase_ex)) )
                    iphempara = np.argmin( np.abs(self.emission_angles_grid - (np.pi+self.validspots[si].phase_ex)) )
                else:
                    iphex = np.argmin( np.abs(self.excitation_angles_grid - self.validspots[si].phase_ex) )
                    # for parallel detection, the phase is the same
                    iphempara = np.argmin( np.abs(self.emission_angles_grid - self.validspots[si].phase_ex) )

                # for perpendicular detection, we need to find the phase+90deg
                iphemperp = np.argmin( np.abs(self.emission_angles_grid - (self.validspots[si].phase_ex+np.pi/2)) )
                Ipara = self.validspots[si].sam[ iphex, iphempara ]
                Iperp = self.validspots[si].sam[ iphex, iphemperp ]
            
            if not float(Ipara+2*Iperp)==0:
                self.validspots[si].r = float(Ipara-Iperp)/float(Ipara+2*Iperp)
            else:
                self.validspots[si].r = np.nan
            # store in contrast image
            self.r_image[ a:b, c:d ] = self.validspots[si].r
#            del(self.validspots[si].sam)  # don't need that anymore


    def ETrulerFFT_selective( self, myspots, slope=7, newdatalength=2048 ):
        # we re-sample the 2D portrait matrix along a slanted line to
        # get a 1D array which contains information about both angular
        # dimensions

        # first we compute the indices into the 2D portrait matrix,
        # this depends only on the angular grids and is therefore the same for
        # all spots and portraits.
        
        # along one direction we increase in single steps along the angular grid:
        ind_em = 1*      np.linspace(0,newdatalength-1,newdatalength).astype(np.int)
        # along the other we have a slope
        ind_ex = slope * np.linspace(0,newdatalength-1,newdatalength).astype(np.int)

        # both index arrays need to wrap around their angular axes
        ind_em = np.mod( ind_em, self.emission_angles_grid.size-1 )
        ind_ex = np.mod( ind_ex, self.excitation_angles_grid.size-1 )
        
        # Now we use these to get the new data.
        # We could do this for every portrait of every spot, but we'll 
        # restrict ourselves to the average portrait of each spot.

        # the angular grids likely include redundancy at the edges, and
        # we need to exclude those to not oversample
        newdata = []
        for si in myspots:
            sam = self.validspots[si].sam#recover_average_portrait_matrix()
            if self.excitation_angles_grid[-1]==np.pi:
                sam = sam[:,:-1]
            if self.emission_angles_grid[-1]==np.pi:
                sam = sam[:-1,:]
            newdata.append( sam[ ind_em,ind_ex ] )
        newdata = np.array( newdata )

        # voila, we have new 1d data

        # now we work out the position of the peaks in the FFTs of these new data columns

        f = np.fft.fft( newdata, axis=1 )
        powerspectra = np.real( f*f.conj() )/newdatalength
        normpowerspectra = powerspectra[:,1:newdatalength/2] \
            /np.outer( np.sum(powerspectra[:,1:newdatalength/2],axis=1), np.ones( (1,newdatalength/2-1) ) )

        # first peak index, pops out at newdatalength / grid  (we are awesome.)
        i1 = newdatalength/(self.excitation_angles_grid.size-1.0)
        # second
        i2 = i1*(slope-1)
        i3 = i1*slope
        i4 = i1*(slope+1)
        df = i1/3
        
        self.peaks = np.array( [ np.sum(normpowerspectra[:, np.round(ii-df):np.round(ii+df)], axis=1) \
                      for ii in [i1,i2,i3,i4] ] )

        # now go over all spots
        cet=0
        for sii,si in enumerate(myspots):

            # if we deviate from the normalized sum by more than 5%,
            # we shouldn't use this ruler
            if np.abs(np.sum( self.peaks[:,sii] )-1) > .08:
                print 'fuck. Data peaks are weird... %f' % (np.sum(self.peaks[:,sii]))
                self.validspots[si].ET_ruler = np.nan
            
            # now let's rule
            crossdiff = self.peaks[1,sii]-self.peaks[3,sii]
            
            # 3-dipole model (all of same length, no ET) starts here
            kappa     = .5 * np.arccos( .5*(3*self.validspots[si].M_ex-1) ) 
            alpha     = np.array([ -kappa, 0, kappa ])
            
            phix =   np.linspace(0,newdatalength-1,newdatalength)*2*np.pi/180
            phim = 7*np.linspace(0,newdatalength-1,newdatalength)*2*np.pi/180

            ModelNoET = np.zeros_like( phix )
            for n in range(3):
                ModelNoET[:] += np.cos(phix-alpha[n])**2 * np.cos(phim-alpha[n])**2
            ModelNoET /= 3
            MYY = np.fft.fft(ModelNoET)
            MYpower = np.real( (MYY*MYY.conj()) )/MYY.size

            MYpeaks = np.array( [ np.sum( MYpower[np.round(ii-df):np.round(ii+df)] ) \
                                      for ii in [i1,i2,i3,i4] ] )
            MYpeaks /= np.sum(MYpeaks)
            
            # test again if peaks make sense
            if np.abs(np.sum( MYpeaks )-1) > .08:
                print 'fuck. MYpeaks is off... %f' % (np.sum( MYpeaks ))
                self.validspots[si].ET_ruler = np.nan

            MYcrossdiff = MYpeaks[1]-MYpeaks[3]
            # model done
            
            ruler = 1-(crossdiff/MYcrossdiff)

            if (ruler < -.1) or (ruler > 1.1):
                print "Shit, ruler has gone bonkers (ruler=%f). Spot #%d" % (ruler,si)
                print "Will continue anyways and set ruler to zero or one (whichever is closer)."
                print "You can thank me later."
                cet+=1
                print cet

            if ruler < 0:
                ruler = 0
            if ruler > 1:
                ruler = 1

            self.validspots[si].ET_ruler = ruler
            self.ET_ruler_image[ self.validspots[si].coords[1]:self.validspots[si].coords[3]+1, \
                               self.validspots[si].coords[0]:self.validspots[si].coords[2]+1 ] = ruler


    def ETmodel_selective( self, myspots, fac=1e4, pg=1e-9, epsi=1e-11 ):

        from fitting import fit_portrait_single_funnel_symmetric

        for si in myspots:
            print 'ETmodel fitting spot %d' % si

            # we 'correct' the modulation in excitation to be within 
            # limits of reason (and proper arccos functionality)
            mex = np.clip( self.validspots[si].M_ex, .000001, .999999 )

            a0 = [mex, 0, 1]
            EX, EM = np.meshgrid( self.excitation_angles_grid, self.emission_angles_grid )
            funargs = (EX, EM, self.validspots[si].sam, mex, self.validspots[si].phase_ex, 'fitting')

            LB = [0.001,   -np.pi/2, 0]
            UB = [0.999999, np.pi/2, 2*(1+mex)/(1-mex)*.999]
            print "upper limit: ", 2*(1+self.validspots[si].M_ex)/(1-self.validspots[si].M_ex)
            print "upper limit (fixed): ", 2*(1+mex)/(1-mex)

            fitresult = so.fmin_l_bfgs_b( func=fit_portrait_single_funnel_symmetric, \
                                      x0=a0, \
                                      fprime=None, \
                                      args=funargs, \
                                      approx_grad=True, \
                                      epsilon=epsi, \
                                      bounds=zip(LB,UB), \
                                      factr=fac, \
                                      pgtol=pg )

            et,A = fit_portrait_single_funnel_symmetric( fitresult[0], EX, EM, self.validspots[si].sam, mex, self.validspots[si].phase_ex, \
                                                             mode='show_et_and_A', use_least_sq=True)
            self.validspots[si].ETmodel_md_fu = fitresult[0][0]
            self.validspots[si].ETmodel_th_fu = fitresult[0][1]
            self.validspots[si].ETmodel_gr    = fitresult[0][2]
            self.validspots[si].ETmodel_et    = et            
            a = self.validspots[si].coords[1]
            b = self.validspots[si].coords[3]+1
            c = self.validspots[si].coords[0]
            d = self.validspots[si].coords[2]+1
            self.ET_model_md_fu_image[ a:b, c:d ] = fitresult[0][0]
            self.ET_model_th_fu_image[ a:b, c:d ] = fitresult[0][1]
            self.ET_model_gr_image[ a:b, c:d ]    = fitresult[0][2]
            self.ET_model_et_image[ a:b, c:d ]    = et

            print 'fit done\t',fitresult[0],
            print ' et=',et,
            print ' A=',A


    def chew_AM( self, quiet=False, loud=False, SNR=10 ):
        self.startstop()
        self.are_spots_valid( SNR=SNR )
        if len(self.validspots)<1:
            raise ValueError("No valid spots found! Reduce SNR demands or re-measure...")

        self.fit_all_portraits_spot_parallel()
        self.find_modulation_depths_and_phases()

        for s in self.validspots:
#            print s
            print "M_ex=%3.2f\tM_em=%3.2f\tphase_ex=%3.2fdeg\tphase_em=%3.2fdeg\tLS=%3.2fdeg" % \
                ( s.M_ex,s.M_em, s.phase_ex*180/np.pi, s.phase_em*180/np.pi, s.LS*180/np.pi )


    def are_spots_valid(self, SNR=10, quiet=False):
        # do we actually have the background std
        bgstd = 0
        if hasattr( self, 'bg_spot' ):
            bgstd = self.bg_spot.std
        else:
            print "Dude --- no background spot defined, therefore no standard deviation. Will treat all spots as valid (i.e. as having sufficient intensity)."

        # then we create a new list containing valid spots only
        validspots = []   
        validspotindices = []
        for si,s in enumerate(self.spots):
            s.SNR = s.mean_intensity/bgstd
            if s.SNR > SNR:
                validspots.append(s)
                validspotindices.append(si)
            # store SNR in SNR_image
            s.store_property_in_image( self.SNR_image, 'SNR' )    # let's see if this works...


        # and store in movie object
        self.validspots = validspots
        self.validspotindices = validspotindices        

        if not quiet: print "Got %d valid spots" % len(self.validspots)

