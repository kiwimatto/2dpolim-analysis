import numpy as np
import matplotlib.pyplot as plt
plt.interactive(1)
from pyspec.ccd.files import PrincetonSPEFile
from motors import ExcitationMotor, EmissionMotor
from fitting import CosineFitter


class Movie:
    def __init__( self, \
                      spe_filename, \
                      excitation_motor_filename, \
                      emission_motor_filename, \
                      phase_offset_excitation=0, \
                      datamode='validdata'):
        # get those objects going
        self.camera_data      = CameraData( spe_filename )
        self.excitation_motor = ExcitationMotor( excitation_motor_filename, \
                                                     phase_offset_excitation, \
                                                     rotation_direction=-1)
        self.emission_motor   = EmissionMotor( emission_motor_filename )

        # where do we get our time axis from?  
        self.timeaxis = self.camera_data.timestamps

        # init list of spots
        self.spots = []

        # store phase offset
#        self.phase_offset_excitation = phase_offset_excitation

        # set data mode, used by collect_data() and startstop()
        self.datamode = datamode

        # fix the excitation and emission angle grids for now
        self.excitation_angles_grid = np.linspace(0,180,181)
        self.emission_angles_grid = np.linspace(0,180,181)



    def define_background_spot( self, coords ):
        # create new spot object
        s = Spot( coords, label='background area' )
        # work out frame-dependent intensities for that spot
        # - mean:
        meanI  = np.sum( np.sum( \
                self.camera_data.rawdata[:, coords[0]:coords[2]+1, coords[1]:coords[3]+1], \
                    axis=2), axis=1 )
        meanI /= s.width*s.height
        # - maximum:
        maxI = np.max( np.max( \
                self.camera_data.rawdata[:, coords[0]:coords[2]+1, coords[1]:coords[3]+1], \
                    axis=2), axis=1 )
        
        # - minimum:
        minI = np.min( np.min( \
                self.camera_data.rawdata[:, coords[0]:coords[2]+1, coords[1]:coords[3]+1], \
                    axis=2), axis=1 )
        # store all that
        s.min_I  = minI
        s.max_I  = maxI
        s.mean_I = meanI

        self.bg = s
        


    def define_spot( self, coords, intensity_type='mean', bg_correction=True, label=None ):
        """Defines a new spot object and adds it to the list.
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

        FIXME: make sure coordinate definitions do not exceed frame size

        TO-DO: how about being able to specify center+radius ?

        """

        # create new spot object
        s = Spot( coords, label )
    
        # work out frame-dependent intensities for that spot
        if intensity_type=='mean':
            I  = np.sum( np.sum( \
                    self.camera_data.rawdata[:, coords[0]:coords[2]+1, coords[1]:coords[3]+1], \
                        axis=2), axis=1 )
            I /= s.width*s.height
            if bg_correction:
                I -= self.bg.mean_I
        elif intensity_type=='max':
            I = np.max( np.max( \
                    self.camera_data.rawdata[:, coords[0]:coords[2]+1, coords[1]:coords[3]+1], \
                        axis=2), axis=1 )
            if bg_correction:
                I -= self.bg.max_I
        elif intensity_type=='min':
            I = np.min( np.min( \
                    self.camera_data.rawdata[:, coords[0]:coords[2]+1, coords[1]:coords[3]+1], \
                        axis=2), axis=1 )
            if bg_correction:
                I -= self.bg.min_I
        else:
            raise ValueError("Don't know what to do for intensity_type %s" % (intensity_type))

        # store all that, subtract bg
        s.intensity_type = intensity_type
        s.bg_corrected = bg_correction
        s.intensity = I
        # append spot object to spots list
        self.spots.append( s )


    def collect_data( self ):
        """This is a helper-function which collects all the necessary 
        information for further analysis in one array.

        Outputs:
        If mode is 'truedata' then an array with columns:
        [FrameNumber, Validity, excitation angle, emission angle, Intensities (Nspot columns)]
        If mode is 'validdate' then:
        [FrameNumber, excitation angle, emission angle, Intensities (Nspot columns)]    
        """

        exangles = np.array( [self.excitation_motor.angle(t) for t in self.timeaxis] )
        emangles = np.array( [self.emission_motor.angle(t) for t in self.timeaxis] )
        validframes = emangles != -1
        self.Nvalidframes = np.sum(validframes)

        Intensity = np.zeros( (self.timeaxis.size, len(self.spots)) )
        for i,s in enumerate(self.spots):
            Intensity[:,i] = self.spots[i].intensity
        
        # if Intensity.ndim==1:
        #     Nspots = 1
        #     Intensity = Intensity.reshape( (Intensity.size,1) )
        # elif Intensity.ndim>2:
        #     raise hell
        # else:
        #     Nspots = Intensity.shape[1]

        Nspots = Intensity.shape[1]

        # make sure that shapes of Intensity and timeaxis match
        assert Intensity.shape[0] == self.timeaxis.size    ### FIXME: is this still needed?

        if self.datamode=='truedata':
            exangles = np.round(exangles, decimals=1)
            emangles = np.round(emangles, decimals=1)

            output = np.zeros( (self.timeaxis.size, 4+Nspots) )
            for i in range( self.timeaxis.size ):
                output[i,0] = i
                output[i,1] = np.int( validframes[i] )
                output[i,2] = exangles[i]
                output[i,3] = emangles[i]
                output[i,4:4+Nspots] = Intensity[i,:]

        elif self.datamode=='validdata':
            emangles = np.round(emangles[validframes], decimals=1)
            exangles = np.round(exangles[validframes], decimals=1)

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
        """This function determines the indices at which portraits start and end.
        There's a fair bit of hackery in here, and to understand what is 
        happening one really needs to look at the raw data, frame indices, valid
        frames (no shutter), etc...  This will be need to be documented much more
        thoroughly in the future, but for now: Handle with care.
        """

        emangles = np.array( [self.emission_motor.angle(t) for t in self.timeaxis] )
        validframes = emangles != -1
    
        emangles_rounded_valid = np.round(emangles[validframes], decimals=1)
        number_of_lines = np.unique( emangles_rounded_valid ).size

        # edge 'detection'
        d = np.diff( emangles_rounded_valid )
        d[0] = 1
        d[d!=0] = 1
        d = np.concatenate( (d, np.array([1])) )
        edges    = d.nonzero()[0]
        # average edge width
        avewidth = np.mean( np.diff( edges )[:-1] )

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

        if self.datamode=="truedata":
            # look up indices into all (not just valid) frames
            indices = validframes.nonzero()[0][indices]
        elif self.datamode=="validdata":
            pass
        else:
            raise ValueError("You screwed up defining the datamode: %s" % (self.datamode))

        self.portrait_indices = indices


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
                d = self.data[pind[n,0]:pind[n,1], :]
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
                portraitlist.append( Portrait( ex, em, I ) )

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


    def fit_all_portraits( self, evaluate_pictures=True ):
        # we assume that the number of portraits and lines is the same 
        # for all spots (can't think of a reason why that shouldn't be the case).
        Nportraits = self.portrait_indices.shape[0]
        Nlines     = len( self.spots[0].portraits[0].lines )

        # for each portrait ---- outermost loop, we do portraits in series
        for pi in range(Nportraits):

            # part I, 'horizontal fitting' of the lines of constant emission angles

            # for each line ---- we do lines in series, __but all spots in parallel__:
            for li in range(Nlines):

                # get excitation angle array (same for all spots!)
                exangles    = self.spots[0].portraits[pi].lines[li].exangles

                # create list of intensity arrays (one array for each spot)
                intensities = [spot.portraits[pi].lines[li].intensities for spot in self.spots]
                # turn into numpy array and transpose
                intensities = np.array( intensities ).T

                phase, I0, M, resi = CosineFitter( exangles, intensities ) 

                # write cosine parameters into line object
                for si in range(len(self.spots)):
                    self.spots[si].portraits[pi].lines[li].set_fit_params( phase[si], I0[si], M[si], resi[si] )

            # part II, 'vertical fitting' --- we do each spot by itself, but 
            # fit all verticals in parallel

            for si in range(len(self.spots)):
                # collect list of unique emission angles
                emangles = [l.emangle for l in self.spots[si].portraits[pi].lines]
                # turn into array, transpose and squeeze
                emangles = np.squeeze(np.array( emangles ).T)
                
                # evaluate cosine-fit at these em_angles, on a grid of ex_angles:
                fitintensities = [ l.cosValue( self.excitation_angles_grid ) \
                                     for l in self.spots[si].portraits[pi].lines ]
                fitintensities = np.array( fitintensities )

                phase, I0, M, resi = CosineFitter( emangles, fitintensities ) 
                
                # store vertical fit params
                self.spots[si].portraits[pi].vertical_fit_params = [phase, I0, M, resi]

                if evaluate_pictures:
                    # evaluate picture
                    mycos = lambda a, ph, I, M: I*(1+M*(np.cos(2*(a-ph)*np.pi/180.0)))
                    pic = np.zeros( (self.emission_angles_grid.size, self.excitation_angles_grid.size) )
                    for exi in range( self.excitation_angles_grid.size ):
                        pic[:,exi] = mycos( self.emission_angles_grid, phase[exi], I0[exi], M[exi] )

                    self.spots[si].portraits[pi].picture = pic




    # def perform_fit( self ):
    #     mycos = lambda a, ph, I, M: I*(1+M*(np.cos(2*(a-ph)*np.pi/180.0)))

    #     # fit each line with a cos**2
    #     for l in self.lines:
    #         # fit line
    #         phase, I0, M, resi = CosineFitter( l.exangles, l.intensities ) 
    #         # store results
    #         l.cosParams( phase, I0, M, resi )

    #     # now create grid for excitation angles
    #     exangs = self.excitation_angles_grid

    #     # and get 'data' from the fits of the lines

    #     ### FIXMEFIXMEFIXME:  - here we only support single spots at the moment !!!
    #     ###                   - check Lines.cosValue for correct behaviour !!!
    #     ###                   - data should have shape [Nemangles, Nexangles, Nspots] !!!
        
    #     data   = np.zeros( (len(self.lines), exangs.size+1) )
    #     for i,l in enumerate(self.lines):
    #         data[i,:] = np.hstack( (l.emangle, l.cosValue(exangs)) )
    #     # then fit all those lines in one go, because the em-angles
    #     # at which data is defined, are the same for all 
    #     phases, I0s, M0s, resis = CosineFitter( data[:,0], data[:,1:] ) 
    #     self.phases = phases
    #     self.I0s    = I0s
    #     self.M0s    = M0s
    #     self.resis  = resis

    #     emangs = self.emission_angles_grid
    #     # take a picture:
    #     # we set this up such that each column corresponds to one excitation angle 
    #     # and each row corresponds to one emission angle
    #     picture = np.zeros( (emangs.size, exangs.size) )
    #     for i in range(picture.shape[1]):
    #         picture[:,i] = mycos( emangs, phases[i], I0s[i], M0s[i] )

    #     self.picture = picture



    def find_modulation_depths_and_phases( self ):
        for spot in self.spots:
            averageportrait = np.zeros_like( spot.portraits[0].picture )
            for p in spot.portraits:
                averageportrait += p.picture
            averageportrait /= len(spot.portraits)

            spot.averageportrait = averageportrait
            ap = averageportrait

            # excitation
            a = self.excitation_angles_grid
            d = np.mean(ap, axis=0)
            assert a.size==d.size
            ph, I, M, r   = CosineFitter( a, d ) 
            spot.phase_ex = ph 
            spot.M_ex     = M

            # emission
            a = self.emission_angles_grid
            d = np.mean(ap, axis=1)
            assert a.size==d.size
            ph, I, M, r   = CosineFitter( a, d ) 
            spot.phase_em = ph
            spot.M_em     = M
       
            # luminescence shift
            LS = spot.phase_ex - spot.phase_em
            if LS>90:
                LS -= 180
            elif LS<-90:
                LS += 180
            spot.LS = LS


    def chew( self, quiet=False, loud=False ):
        self.collect_data()
        self.startstop()
        self.assign_portrait_data()
        self.fit_all_portraits()
        self.find_modulation_depths_and_phases()
        

        if not quiet:
            print "Quick report:"
            print "-------------"
            print "Number of frames: %d"  % (self.timeaxis.size)
            print "Number of valid frames: %d" % (self.Nvalidframes)
            print "Number of spots: %d" % (len(self.spots))
            print "Number of portraits: %d" % (len(self.spots[0].portraits))
            for i,p in enumerate( self.spots[0].portraits ):
                print "Portrait #%d (%d frames)" % (i,p.exangles.size)
        if loud:
            for si,spot in enumerate( self.spots ):
                print "Spot #%d:" % (si)
                print "\tModulation depth in excitation: %f" % (spot.M_ex)
                print "\tPhase in excitation: %f" % (spot.phase_ex)
                print "\tModulation depth in emission:   %f" % (spot.M_em)
                print "\tPhase in emission  : %f" % (spot.phase_em)
                print "\tLS: %f" % (spot.LS)


    def show_average_picture( self ):
        plt.matshow( self.averageportrait, origin='bottom')
        plt.plot( [0,180], [0,180], 'k-' )
        plt.xlim( 0, 180 )
        plt.ylim( 0, 180 )      





class CameraData:
    def __init__( self, spe_filename ):
        # load SPE  ---- this will work for SPE format version 2.5 (probably not for 3...)
        self.filename           = spe_filename
        self.rawdata_fileobject = PrincetonSPEFile( self.filename )
        self.rawdata            = self.rawdata_fileobject.getData().astype(np.float64)
        self.exposuretime       = self.rawdata_fileobject.Exposure   # in seconds

        ###  extract or generate time stamps ###
        #  here we do not have timestamps for each frame, so we
        #  generate time axis from SPE exposure data and number of frames
        self.timestamps = np.linspace( 0, self.exposuretime*self.rawdata.shape[0], \
                                           self.rawdata.shape[0], endpoint=False )


class Spot:
    def __init__(self, coords, label):
        self.coords = coords    #[left, bottom, right, top]
        self.label  = label
        self.width  = coords[2] -coords[0] +1
        self.height = coords[3] -coords[1] +1
        self.intensity_type = None
        self.intensity      = None
        self.bg_corrected   = None



class Portrait:
    def __init__( self, exangles, emangles, intensities ):
        self.exangles = exangles
        self.emangles = emangles
        self.intensities = intensities

        self.split_into_emission_lines()

    def split_into_emission_lines( self ):
        """Splits the portrait into the flat emission parts."""
    
        # edge 'detection' to work out where emission angles change
        edges = (np.diff(self.emangles)!=0).nonzero()[0]
        # extend the array: a zero on the left because we know that the
        # first frame is an edge, and a p.shape[0] on the right because
        # we know that the last frame is an edge as well.
        edges = np.concatenate( (np.array( [0] ), edges+1, np.array([self.exangles.shape[0]]) ) )

        # grab the angles and intensities associated with a line, where by 'line'
        # we mean data points for which emission angles do not change (but excitation 
        # does!)
        lines = []
        for i in range(edges.size-1):
            exangles  = self.exangles[edges[i]:edges[i+1]]
            intensity = self.intensities[edges[i]:edges[i+1]]
            emangles  = self.emangles[edges[i]:edges[i+1]]
        
            # test that this is really one emission line (no change in 
            # emission angles)
            assert np.unique( emangles ).size == 1

            l = Line( exangles, intensity, np.unique(emangles) )

            # append line to list of lines
            lines.append( l )
#        lines.append( np.array((exangles, intensity)).T )

        # did we lose any data rows due to our funky indexing?
        assert np.sum( [lines[i].exangles.shape[0] for i in range(len(lines))] ) == self.exangles.shape[0]
    
        self.lines = lines

    def show_picture( self ):
        plt.matshow( self.picture, origin='bottom')
        plt.plot( [0,180], [0,180], 'k-' )
        plt.xlim( 0, 180 )
        plt.ylim( 0, 180 )      



class Line:
    def __init__( self, exangles, intensities, emangle ):
        self.exangles = exangles
        self.intensities = intensities
        self.emangle = emangle

    def set_fit_params( self, phase, I0, M0, residuals ):
        self.phase = phase
        self.I0 = I0
        self.M0 = M0
        self.resi = residuals
        
    def cosValue( self, angle ):
        return self.I0 * (1+self.M0*np.cos( 2*(angle-self.phase)*np.pi/180.0 ))





