from cameradata import CameraData
from fitting import CosineFitter_new
from spot import Spot
from portrait import Portrait
from motors import *
import os

class Movie:
    def __init__( self, \
                      datadir, basename, \
                      which_setup='cool new setup', \
                      phase_offset_in_deg=np.nan, \
                      excitation_optical_element='l/2 plate' ):

        self.cos_fitter = CosineFitter_new
        self.which_setup = which_setup
        self.phase_offset_in_deg = phase_offset_in_deg

        self.data_directory = datadir
        self.data_basename  = basename

        self.read_in_EVERYTHING()

        self.spots = []
        self.initContrastImages()
        self.datamode = 'validdata'

        # fix the excitation and emission angle grids for now
        self.excitation_angles_grid = np.linspace(0,np.pi,91)
        self.emission_angles_grid = np.linspace(0,np.pi,91)
        self.Nphases_for_cos_fitter = 91
        self.precomputed_cosines = [ np.cos(2*(self.emission_angles_grid-ph)) \
                                         for ph in np.linspace(0,np.pi/2,self.Nphases_for_cos_fitter) ]

    def __init__2( self, \
                      spe_filename, \
                      excitation_motor_filename, \
                      emission_motor_filename=None, \
                      blank_sample_filename=None, \
                      phase_offset_excitation=0, \
                      datamode='validdata', \
                      which_setup='new setup', \
                      use_new_fitter=True, \
                      excitation_optical_element='l/2 plate'):

        self.camera_data    = CameraData( spe_filename, compute_frame_average=True )

        if use_new_fitter:
            self.cos_fitter = CosineFitter_new
        else:
            self.cos_fitter = CosineFitter

        self.which_setup = which_setup
            
        # set up motors --- phase offset in radians!!!
        if which_setup=='old setup':
            self.excitation_motor = ExcitationMotor( excitation_motor_filename, \
                                                         phase_offset_excitation, \
                                                         rotation_direction=-1, \
                                                         optical_element=excitation_optical_element)
            self.emission_motor   = EmissionMotor( emission_motor_filename )

        elif which_setup=='new setup':
            self.excitation_motor = NewSetupMotor( excitation_motor_filename, \
                                                       which_motor='excitation', \
                                                       phase_offset=phase_offset_excitation, \
                                                       optical_element=excitation_optical_element )
            if emission_motor_filename==None:
                emission_motor_filename = excitation_motor_filename
            self.emission_motor = NewSetupMotor( emission_motor_filename, \
                                                     which_motor='emission', \
                                                     phase_offset=0*np.pi/180.0 )
        elif which_setup=='cool new setup':
            self.motors = BothMotors( excitation_motor_filename, \
                                          phase_offset_excitation, \
                                          optical_element=excitation_optical_element )
        else:
            raise hell

        # where do we get our time axis from?  
        self.timeaxis = self.camera_data.timestamps

        # init list of spots
        self.spots = []

        # init contrast images
        self.initContrastImages()

        # store phase offset
        self.phase_offset_excitation = phase_offset_excitation

        # set data mode, used by collect_data() and startstop()
        self.datamode = datamode

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
        
        ### old setup ###
        # the old setup is the only with two separate files, 
        # so if this is the case we'll do it now
        if self.which_setup=='old setup':
            got_motors = False
            got_motor_file_ex = False
            got_motor_file_em = False
            for file in os.listdir("."):
                if (not got_motor_file_ex) and (file=="MSex-"+self.data_basename+'.txt'):
                    print '\t found old-style motor file (for excitation) %s' % file
                    self.motorfile_ex = file
                    got_motor_file_ex = True
                elif (not got_motor_file_em) and (file=="MSem-"+self.data_basename+'.txt'):
                    print '\t found old-style motor file (for emission) %s' % file
                    self.motorfile_em = file
                    got_motor_file_em = True
                # break if both files were found
                if got_motor_file_ex and got_motor_file_em:
                    got_motors = True
                    break
                
            if got_motors:
                # read in motor files
                self.excitation_motor = ExcitationMotor( self.motorfile_ex, \
                                                             self.phase_offset_in_deg, \
                                                             rotation_direction=-1, \
                                                             optical_element=self.excitation_optical_element)
                self.emission_motor   = EmissionMotor( self.motorfile_em )
                # and get the angles
                self.exangles = np.array( [self.excitation_motor.angle(t,exposuretime=self.camera_data.exposuretime) for t in self.timeaxis] )
                self.emangles = np.array( [self.emission_motor.angle(t,exposuretime=self.camera_data.exposuretime) for t in self.timeaxis] )
            else:
                raise IOError("Couldn't find the motor files (old setup was specified).")


        ### newer setups ###
        # the newer versions of the setup all generate just one motor file, let's look for that
        else:
            got_motors = False
            for file in os.listdir("."):
                if file=="MS-"+self.data_basename+'.txt':
                    print '\t found motor file %s' % file
                    self.motorfile = file
                    got_motors = True
                    break
            # now import it depending on which setup it is
            if got_motors:

                if self.which_setup=='new setup':
                    self.motors = NewSetupMotors( self.motorfile, \
                                                      phase_offset=self.phase_offset_in_deg, \
                                                      optical_element=self.excitation_optical_element )
                    self.motors.compute_angles( self.camera_data.timeaxis, self.camera_data.exposuretime )

                elif self.which_setup=='cool new setup':
                    self.motors = BothMotors( self.motorfile, \
                                                  self.phase_offset_in_deg )

                elif self.which_setup=='header':
                    ### Supplying the phase offset here means it overrides the value 
                    ### which is read from the header, unless it is NaN (the default). 
                    ### This is only useful for AM measurements.
                    self.motors = BothMotorsWithHeader( self.motorfile, \
                                                            self.phase_offset_in_deg )

                else:
                    raise ValueError('Dunno what setup="%s" is.' % self.which_setup)
                
                self.exangles = self.motors.excitation_angles
                self.emangles = self.motors.emission_angles                   
 

        if not got_motors:
            raise IOError("Couldn't find motor files.")
       
        ###### look for blank sample ######
        print 'Looking for blank...',
        for file in os.listdir("."):
            if file.startswith("blank-") and (file.endswith(".spe") or file.endswith(".SPE")):
                print '\t found file %s' % file
                b = CameraData( file, compute_frame_average=True )
                self.blank = b.average_image.copy()
                del(b)
                break


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
        if hasattr(self, 'blanks'):
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

    def collect_data( self ):
        """This is a helper-function which collects all the necessary 
        information for further analysis in one array.

        Outputs:
        If mode is 'truedata' then an array with columns:
        [FrameNumber, Validity, excitation angle, emission angle, Intensities (Nspot columns)]
        If mode is 'validdate' then:
        [FrameNumber, excitation angle, emission angle, Intensities (Nspot columns)]    
        """

#        if self.which_setup=='cool new setup':
        if hasattr( self, 'motors'):
            exangles = self.motors.excitation_angles
            emangles = self.motors.emission_angles
        else:
            exangles = np.array( [self.excitation_motor.angle(t,exposuretime=self.camera_data.exposuretime) \
                                      for t in self.timeaxis] )
            emangles = np.array( [self.emission_motor.angle(t,exposuretime=self.camera_data.exposuretime) \
                                      for t in self.timeaxis] )
        validframes = emangles != -1
        self.Nvalidframes = np.sum(validframes)

        Intensity = np.zeros( (self.timeaxis.size, len(self.spots)) )
        for i,s in enumerate(self.spots):
            Intensity[:,i] = self.spots[i].intensity
            del( self.spots[i].intensity )

        Nspots = Intensity.shape[1]

        # make sure that shapes of Intensity and timeaxis match
        assert Intensity.shape[0] == self.timeaxis.size    ### FIXME: is this still needed?

        if self.datamode=='truedata':
            exangles = np.round(exangles, decimals=2)
            emangles = np.round(emangles, decimals=2)

            output = np.zeros( (self.timeaxis.size, 4+Nspots) )
            for i in range( self.timeaxis.size ):
                output[i,0] = i
                output[i,1] = np.int( validframes[i] )
                output[i,2] = exangles[i]
                output[i,3] = emangles[i]
                output[i,4:4+Nspots] = Intensity[i,:]

        elif self.datamode=='validdata':
            emangles = np.round(emangles[validframes], decimals=2)
            exangles = np.round(exangles[validframes], decimals=2)

            output = np.zeros( (self.Nvalidframes, 3+Nspots) )
            trueindices = validframes.nonzero()[0]
            for i in range( self.timeaxis[validframes].size ):            
                output[i,0] = trueindices[i]
                output[i,1] = exangles[i]
                output[i,2] = emangles[i]
                output[i,3:3+Nspots] = Intensity[trueindices[i],:]
        else:
            raise ValueError("Don't understand datamode: %s" % (self.datamode))
                
        self.data = output


    def startstop( self ):
        """
        This function determines the indices at which portraits start and end.
        There's a fair bit of hackery in here, and to understand what is 
        happening one really needs to look at the raw data, frame indices, valid
        frames (no shutter), etc...  This will be need to be documented much more
        thoroughly in the future, but for now: Handle with care.
        """

        if hasattr( self, 'motors' ):
            emangles = self.motors.emission_angles
        else:
            emangles = np.array( [self.emission_motor.angle(t, exposuretime=self.camera_data.exposuretime) \
                                      for t in self.timeaxis] )

        # frames are valid where emangles is not equal to -1
        validframes = emangles != -1
    
        emangles_rounded_valid = np.round(emangles[validframes], decimals=2)

        number_of_lines = np.unique( emangles_rounded_valid ).size

        # edge 'detection' via diff
        d = np.diff( emangles_rounded_valid )
        d[0] = 1
        d[d!=0] = 1
        d = np.concatenate( (d, np.array([1])) )
        edges    = d.nonzero()[0]
        
        # average gap size between edges
        avewidth = np.mean( np.diff( edges )[:-1] )
        # if the last gap size differs more than ten percent from the average,
        # remove it
        if not avewidth*.9 <= np.diff( edges )[-1] <= avewidth*1.1:
            edges = edges[:-1]

        # integer division to find out how many complete portraits we have
        Nportraits = np.diff(edges).size / number_of_lines 

        indices = np.zeros( (Nportraits,2), dtype=np.int )
        for i in range(Nportraits):
            indices[i,0] = edges[i*number_of_lines]+1
            indices[i,1] = edges[(i+1)*number_of_lines]

        # These indices will discard the first frame, despite it being in principle
        # a valid frame. You can 'repair' that by doing:
        #    indices[0,0] = 0
        # However, we think that the first frame could be affected by the shutter, and
        # the overall behaviour of the code is a bit more consistent this way.
        if self.which_setup=='cool new setup':
            indices[0,0] = 0

        #print indices

        if self.datamode=="truedata":
            # look up indices into all (not just valid) frames
            indices = validframes.nonzero()[0][indices]
        elif self.datamode=="validdata":
            pass
        else:
            raise ValueError("You screwed up defining the datamode: %s" % (self.datamode))

        self.portrait_indices = indices

#        return emangles, emangles_rounded_valid, d


    def assign_portrait_data( self ):  #startstop, data, mode ):
        """Generates portrait list _for each spot_. 
        Each list element is a full portrait. 
        We basically split the output of collect_data(), using the indices
        provided by startstop(), and the spot index.
        The output is independent of _mode_, it only contains valid data.
        """
        pind = self.portrait_indices
        Nportraits = pind.shape[0]

        for si,spot in enumerate(self.spots):
            # create empty list
            portraitlist = []
            # go through all portraits
            for n in range(Nportraits):
                # grab nth portrait
                d = self.data[pind[n,0]:pind[n,1]+1, :]
                if self.datamode=='truedata':
                    # use only valid rows
                    d = d[d[:,1]==1]
                    # from the remaining (valid) rows, grab
                    # columns ex, em and Int
                    ex = d[:, 2]
                    em = d[:, 3]
                    I  = d[:, 4+si]
                elif self.datamode=='validdata':            
                    # simply grab columns ex, em and Int
                    ex = d[:, 1]
                    em = d[:, 2]
                    I  = d[:, 3+si]
                # store portrait object in list
                portraitlist.append( Portrait( ex, em, I, spot ) )

            # store list in spot object
            spot.portraits = portraitlist
        

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
        

    def fit_all_portraits_spot_parallel( self ):

        # init average portrait matrices, so that we can write to them without
        # having to store a matrix for each portrait
        for s in self.validspots:
            s.residual = 0

        # we assume that the number of portraits and lines is the same 
        # for all spots (can't think of a reason why that shouldn't be the case).
        Nportraits = self.portrait_indices.shape[0]
        Nlines     = len( self.validspots[0].portraits[0].lines )

        # for each portrait ---- outermost loop, we do portraits in series
        for pi in range(Nportraits):

            # part I, 'horizontal fitting' of the lines of constant emission angles

            # for each line ---- we do lines in series, __but all spots in parallel__:
            for li in range(Nlines):

                # get excitation angle array (same for all spots!)
                exangles    = self.validspots[0].portraits[pi].lines[li].exangles

                # create list of intensity arrays (one array for each spot)
                intensities = [spot.portraits[pi].lines[li].intensities for spot in self.validspots]

                # turn into numpy array and transpose
                intensities = np.array( intensities ).T

                exa = exangles.copy()

                phase, I0, M, resi, fit, rawfitpars, mm = self.cos_fitter( exa, intensities, \
                                                                               self.Nphases_for_cos_fitter ) 

                # write cosine parameters into line object
                for si in range(len(self.validspots)):
                    self.validspots[si].portraits[pi].lines[li].set_fit_params( phase[si], I0[si], M[si], resi[si] )

            # gather residuals for this protrait
            for si in range(len(self.validspots)):
                self.validspots[si].residual = np.sum( [ l.resi for l in self.validspots[si].portraits[pi].lines ] )

            # part II, 'vertical fitting' --- we do each spot by itself, but 
            # fit all verticals in parallel

            # collect list of unique emission angles (same for all spots!)                    
            emangles = [l.emangle for l in self.validspots[0].portraits[pi].lines]
            # turn into array, transpose and squeeze
            emangles = np.squeeze(np.array( emangles ).T)

            # evaluate cosine-fit at these em_angles, on a grid of ex_angles:
            fitintensities = [ np.array( [l.cosValue( self.excitation_angles_grid ) \
                                             for l in s.portraits[pi].lines]) for s in self.validspots ]
            fitintensities = np.hstack( fitintensities )

            phase, I0, M, resi, fit, rawfitpars, mm = self.cos_fitter( emangles, fitintensities, \
                                                                           self.Nphases_for_cos_fitter ) 
                
            # store vertical fit params
            phase = np.hsplit(phase, len(self.validspots))
            I0    = np.hsplit(I0, len(self.validspots))
            M     = np.hsplit(M, len(self.validspots))
            resi  = np.hsplit(resi, len(self.validspots))
            mm    = np.hsplit(mm, len(self.validspots))
            for si,s in enumerate(self.validspots):
                s.portraits[pi].vertical_fit_params = [ phase[si], I0[si], M[si], resi[si], mm[si] ]


    def find_modulation_depths_and_phases( self ):
        # projection onto the excitation axis (ie over all emission angles),
        # for all spots, one per column
        proj_ex = []
        proj_em = []
        for s in self.validspots:
            sam  = s.recover_average_portrait_matrix()
            # sam /= np.max(sam)*255
            # s.sam = sam.astype( np.uint8 )
            s.sam = sam
            s.proj_ex = np.mean( sam, axis=0 )
            s.proj_em = np.mean( sam, axis=1 )
            proj_ex.append( s.proj_ex )
            proj_em.append( s.proj_em )
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
        for si,s in enumerate(self.validspots):
            s.phase_ex = ph_ex[si]
            s.M_ex     = M_ex[si]
            s.phase_em = ph_em[si]
            s.M_em     = M_em[si]
            s.LS       = LS[si]
            # store in coverage maps
            self.M_ex_image[ s.coords[1]:s.coords[3]+1, s.coords[0]:s.coords[2]+1 ] = s.M_ex
            self.M_em_image[ s.coords[1]:s.coords[3]+1, s.coords[0]:s.coords[2]+1 ] = s.M_em
            self.phase_ex_image[ s.coords[1]:s.coords[3]+1, s.coords[0]:s.coords[2]+1 ] = s.phase_ex
            self.phase_em_image[ s.coords[1]:s.coords[3]+1, s.coords[0]:s.coords[2]+1 ] = s.phase_em
            self.LS_image[ s.coords[1]:s.coords[3]+1, s.coords[0]:s.coords[2]+1 ] = s.LS        

            #### anisotropy ####
        
            if (s.M_ex < .15):
                iex = np.argmin( np.abs( self.excitation_angles_grid - 90 ) )
                iem = np.argmin( np.abs( self.emission_angles_grid - 90 ) )
                Ipara = s.sam[ iex, iem ]
                Iperp = s.sam[ iex, 0 ] ## assumes that the first index is close to 0deg in emission angle grid
            else:
                # find where the found phase matches the phase in the portrait matrix
                # portrait matrices are constructed from the angle grid, so look up which 
                # index corresponds to the found phase
                iphex = np.argmin( np.abs(self.excitation_angles_grid - s.phase_ex) )
                # for parallel detection, the phase is the same
                iphempara = np.argmin( np.abs(self.emission_angles_grid - s.phase_ex) )
                # for perpendicular detection, we need to find the phase+90deg
                iphemperp = np.argmin( np.abs(self.emission_angles_grid - np.mod(s.phase_ex+90, 180)) )
                Ipara = s.sam[ iphex, iphempara ]
                Iperp = s.sam[ iphex, iphemperp ]
            
            if not float(Ipara+2*Iperp)==0:
                s.r = float(Ipara-Iperp)/float(Ipara+2*Iperp)
            else:
                s.r = np.nan
            # store in contrast image
            self.r_image[ s.coords[1]:s.coords[3]+1, s.coords[0]:s.coords[2]+1 ] = s.r
            del(s.sam)  # don't need that anymore


    def ETrulerFFT( self, slope=7, newdatalength=2048 ):
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
        for s in self.validspots:
            sam = s.recover_average_portrait_matrix()
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
        for si,s in enumerate(self.validspots):

            # if we deviate from the normalized sum by more than 5%,
            # we shouldn't use this ruler
            if np.abs(np.sum( self.peaks[:,si] )-1) > .08:
                print 'fuck. Data peaks are weird... %f' % (np.sum(self.peaks[:,si]))
                self.validspots[si].ET_ruler = np.nan
            
            # now let's rule
            crossdiff = self.peaks[1,si]-self.peaks[3,si]
            
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


    def ETmodel( self, fac=1e4, pg=1e-9, epsi=1e-11 ):

        from fitting import fit_portrait_single_funnel_symmetric

        for si,s in enumerate(self.validspots):
            print 'ETmodel fitting spot %d' % si

            # we 'correct' the modulation in excitation to be within 
            # limits of reason (and proper arccos functionality)
            mex = np.clip( s.M_ex, .000001, .999999 )

            a0 = [mex, 0, 1]
            EX, EM = np.meshgrid( self.excitation_angles_grid, self.emission_angles_grid )
            funargs = (EX, EM, s.averagematrix, mex, s.phase_ex, 'fitting')

            LB = [0.001,   -np.pi/2, 0]
            UB = [0.999999, np.pi/2, 2*(1+mex)/(1-mex)*.999]
            print "upper limit: ", 2*(1+s.M_ex)/(1-s.M_ex)
            print "upper limit (fixed): ", 2*(1+mex)/(1-mex)

            a = so.fmin_l_bfgs_b( func=fit_portrait_single_funnel_symmetric, \
                                      x0=a0, \
                                      fprime=None, \
                                      args=funargs, \
                                      approx_grad=True, \
                                      epsilon=epsi, \
                                      bounds=zip(LB,UB), \
                                      factr=fac, \
                                      pgtol=pg )

            et,A = fit_portrait_single_funnel_symmetric( a[0], EX, EM, s.averagematrix, mex, s.phase_ex, \
                                                             mode='show_et_and_A', use_least_sq=True)
            s.ETmodel_md_fu = a[0][0]
            s.ETmodel_th_fu = a[0][1]
            s.ETmodel_gr    = a[0][2]
            s.ETmodel_et    = et            

            self.ET_model_md_fu_image[ s.coords[1]:s.coords[3]+1, s.coords[0]:s.coords[2]+1 ] = a[0][0]
            self.ET_model_th_fu_image[ s.coords[1]:s.coords[3]+1, s.coords[0]:s.coords[2]+1 ] = a[0][1]
            self.ET_model_gr_image[ s.coords[1]:s.coords[3]+1, s.coords[0]:s.coords[2]+1 ] = a[0][2]
            self.ET_model_et_image[ s.coords[1]:s.coords[3]+1, s.coords[0]:s.coords[2]+1 ] = et

            print 'fit done\t',a[0],
            print ' et=',et,
            print ' A=',A


    def chew( self, quiet=False, loud=False ):
        print "collecting data..."
        self.collect_data()
        print "startstop..."
        self.startstop()
        print "assigning portrait data..."
        self.assign_portrait_data()        
        self.are_spots_valid()

        # for s in self.validspots:       
        #     print s.intensity
        #     s.intensity = None
#        self.data = None             # this doesn't seem to do anything memory-wise
                                     # (are these all just pointers??)
#        self.camera_data = None      # this helps (as expected)

        import time as stopwatch
        print "fitting all portraits"
        self.fit_all_portraits_spot_parallel()

        print "finding mod depths and phases..."
        self.find_modulation_depths_and_phases()
        print "ETruler..."
        self.ETrulerFFT()
        print "ETmodel..."

        tstart = stopwatch.time()
        self.ETmodel()
        print "time taken: %fs" % (stopwatch.time()-tstart)

        if not quiet:
            print "Quick report:"
            print "-------------"
            print "Number of frames: %d"  % (self.timeaxis.size)
            print "Number of valid frames: %d" % (self.Nvalidframes)
            print "Number of spots: %d" % (len(self.validspots))
            print "Number of portraits: %d" % (len(self.validspots[0].portraits))
            for i,p in enumerate( self.validspots[0].portraits ):
                print "Portrait #%d (%d frames)" % (i,p.exangles.size)
        if loud:
            for si,spot in enumerate( self.validspots ):
                print "Spot #%d:" % (si)
                print "\tModulation depth in excitation: %f" % (spot.M_ex)
                print "\tPhase in excitation: %f" % (spot.phase_ex)
                print "\tModulation depth in emission:   %f" % (spot.M_em)
                print "\tPhase in emission  : %f" % (spot.phase_em)
                print "\tLS: %f" % (spot.LS)


    def chew_a_bit( self, SNR=30, quiet=False, loud=False ):
        self.collect_data()
        self.startstop()
        self.assign_portrait_data()        
        self.are_spots_valid( SNR )
        self.fit_all_portraits_spot_parallel()
        self.find_modulation_depths_and_phases()

    def chew_AM( self, quiet=False, loud=False, SNR=10 ):
        self.collect_data()
        self.startstop()
        self.assign_portrait_data()        
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

