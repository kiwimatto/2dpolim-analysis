import numpy as np
from misc import deal_with_date_time_string


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


class BothMotorsWithHeader:
    def __init__( self, filename, phase_offset_in_deg ):
        """Initialize the class: read in the file"""
        self.experiment_start_datetime = None
        self.filename = filename        
        self.phase_offset_in_deg = phase_offset_in_deg

        # read header if possible
        f = open( self.filename, 'r' )
        self.header = {}
        self.header['comments'] = []
        Nheaderlines = 0
        while True:
            line = f.readline().strip().split(':')
#            print line
            Nheaderlines += 1
            if line[0]=='END-OF-HEADER':
                break
            else:
                if len(line)>1:
                    if is_number(line[1]):
                        self.header[line[0]] = float(line[1])
                    else:
                        self.header[line[0]] = line[1].strip().lower()
                else:
                    self.header['comments'].append( line[0] )
            if Nheaderlines >= 1000:
                print 'Tried to read in 1000 header lines --- either this file has no header or the header contains a novel, you tell me...   (bombing out)'
                raise SystemExit

        f.close()

        # grab motor data --- delimiter is _a tab_ !!!
        # uses converter functions to parse date+time and shutter values
        md = np.loadtxt( filename, \
                             delimiter='\t', \
                             skiprows=Nheaderlines+1 )

        # use externally supplied phase offset only if not NaN
        if not np.isnan( self.phase_offset_in_deg ):
            phase_offset_in_deg = self.phase_offset_in_deg
        else:
            if not self.header.has_key('phase offset in deg'):
                raise ValueError("The header is missing the entry for 'phase offset in deg'. Cannot continue.")
            if not np.isnan(self.header['phase offset in deg']):
                phase_offset_in_deg = self.header['phase offset in deg']
                print "Motorfile: got phase offset of %f deg" % (phase_offset_in_deg)
            else:
                print '\n\n\n'
                print '###################################################################\n'
                print "Cannot continue: the phase offset in the header is NaN.\nThis should only be the case for AM measurements.\nMost likely you forgot to specify the phase offset in the header data (fix this now in the motor file, %s).\nIf this is not the case, and this is an AM measurement, then for some reason the presumed phase offset was not passed to the movie class (that is, you're doing something wrong in the AM software).\n" % (filename)
                print '###################################################################\n'
                raw_input('[got that? press enter]')
                raise ValueError("Phase offset in header is NaN (which should only be the case if this is an AM measurement), but no phase offset was specified to the movie class.")

        self.framenumbers = md[:,0]
#        print self.header['optical element in excitation']

        if md.shape[1]==5:  #### NEW STYLE ANGLES!
            self.excitation_angles        = (md[:,2] + phase_offset_in_deg ) * np.pi/180.0
            self.emission_angles          = md[:,4] * np.pi/180.0
            self.sample_plane_intensities = md[:,3]
        elif md.shape[1]==3:
            if not self.header.has_key('optical element in excitation'):
                raise ValueError("The header is missing the entry 'optical element in excitation'. Cannot continue.")
            if self.header['optical element in excitation']=='l/2 plate':
                print 'Header says that l/2 plate was used.'
                self.excitation_angles = (2*md[:,1] + phase_offset_in_deg ) * np.pi/180.0
            else:
                print 'Err... header says that l/2 plate was _not_ used. Correct??'
                self.excitation_angles = (  md[:,1] + phase_offset_in_deg ) * np.pi/180.0
            self.emission_angles   = md[:,2] * np.pi/180.0
        else:
            raise ValueError("Huh? Something's very wrong here... the number of columns in the motor file appears to be neither five (standard) nor three (older style). Cannot continue.")

        self.excitation_angles     = np.mod( self.excitation_angles, np.pi ) 
        self.emission_angles       = np.mod( self.emission_angles, np.pi )
        self.phase_offset_in_deg   = phase_offset_in_deg
#        print self.excitation_angles







###### OLD MOTOR CLASSES #######




class BothMotors:
    def __init__( self, filename, phase_offset_in_deg=0 ):
        """Initialize the class: read in the file"""
        self.experiment_start_datetime = None
        self.filename = filename
        self.phase_offset_in_deg = phase_offset_in_deg

        # deal with motor file
        f = open(filename,'r')        
        optelemstring = f.readline().strip()   # read first line and strip newline character(s)
        f.close()
        if optelemstring.lower()=='l/2 plate':
            self.optical_element = 'l/2 plate'
        elif optelemstring.lower()=='polarizer':
            self.optical_element = 'polarizer'
        else:
            raise ValueError("BothMotors doesn't understand optelemstring %s (should be 'l/2 plate' or 'polarizer'" % optelemstring)

        # grab motor data --- delimiter is _a tab_ !!!
        # uses converter functions to parse date+time and shutter values
        md = np.loadtxt( filename, \
                             delimiter='\t', \
                             skiprows=2 )

        self.framenumbers = md[:,0]
        if self.optical_element.lower()=='l/2 plate':
            self.excitation_angles = (2*md[:,1] + self.phase_offset_in_deg ) * np.pi/180.0 
        else:
            self.excitation_angles = (  md[:,1] + self.phase_offset_in_deg ) * np.pi/180.0
        self.emission_angles   = md[:,2] * np.pi/180.0

        self.excitation_angles     = np.mod( self.excitation_angles, np.pi )
        self.emission_angles       = np.mod( self.emission_angles, np.pi )




class NewSetupMotors:
    def __init__( self, filename, \
                      phase_offset_in_deg=0, \
                      optical_element='polarizer'):

        """Initialize the class: read in the file"""
        self.experiment_start_datetime = None
        self.filename = filename
#        self.timeaxis = timeaxis
        self.phase_offset_in_deg = phase_offset_in_deg
        self.optical_element = optical_element

        self.load_motor_file( filename )


    def load_motor_file( self, filename ):
        """Reads the labview file written for that motor. 
        It assumes that fields (time, motor position(s), shutter) are separated 
        by _two_ spaces (which seems to be the case in all file generated so far).
        This makes it easy to read in (and parse) the date and time as a single 
        string.
        This function currently also parses the shutter value: open/closed is converted 
        into a boolean with an inline function doing string comparison.

        The function returns the time stamps with corresponding motor position.
        """

        # grab motor data --- delimiter is _a tab_ !!!
        # uses converter functions to parse date+time and shutter values
        md = np.loadtxt( filename, \
                             delimiter='\t', \
                             skiprows=1, \
                             converters={0: lambda s: deal_with_date_time_string(self, s), \
                                             3: lambda s: s=='open'} )
                            
        self.timestamps = md[:,0]
        self.emisangles = md[:,1] * np.pi/180.0
        if self.optical_element.lower()=='l/2 plate':
            self.exciangles = (2*md[:,2] + self.phase_offset_in_deg ) * np.pi/180.0
        elif self.optical_element.lower()=='polarizer':
            self.exciangles = (  md[:,2] + self.phase_offset_in_deg ) * np.pi/180.0
        else:
            raise hell
        self.shutter    = md[:,3]


    def compute_angles( self, timeaxis, exposuretime=0.1, respectShutter=True, raw=False ):
        """This function will return the angle of the emission polarizer given
        a time (in seconds, counting from the start of the experiment).
        Since the emission polarizer is turning in a very discrete fashion (0deg, 45deg, etc),
        we can assume that we can return the angle which is closest to the 
        the queried time (nearest interpolation), with the following caveats:
        if the any of the nearest-neighbor timestamps reports the shutter as closed, 
        then we return -1.
        """
        self.excitation_angles = []
        self.emission_angles = []
        for time in timeaxis:
            if respectShutter:
                ### find all timestamps within one exposure time ###
                mintime = time - exposuretime
                maxtime = time + exposuretime
                # left boundary of the window: move time axis by mintime
                tt = self.timestamps - mintime
                # last element which is less than zero must be 
                # the timestamp just before the window starts,
                # so the length of the array tt[tt>0] is the index
                # for that timestamp, and adding one to it we get
                # the first timestamp in the window. Since indices
                # start with zero, we subtract that one again and 
                # we have simply:
                first = tt[tt<0].size
            
                # right boundary,
                tt = self.timestamps - maxtime
                last = tt[tt<0].size -1

                # now test if the shutter is off at any of the timestamps
                # between first and last:
                if np.any( self.shutter[first:last]==0 ):
                    self.excitation_angles.append( -1 )
                    self.emission_angles.append( -1 )
                else:
                    # interpolate angle from known angles
                    self.excitation_angles.append( \
                        np.mod( \
                            np.interp( time, self.timestamps, self.exciangles ) \
                                + self.phase_offset_in_deg*np.pi/180.0, np.pi ) )
                    self.emission_angles.append( \
                        np.mod( \
                            np.interp( time, self.timestamps, self.emisangles ), np.pi ) )
            else:
                self.excitation_angles.append( \
                    np.mod( \
                        np.interp( time, self.timestamps, self.exciangles ) \
                            + self.phase_offset_in_deg*np.pi/180.0, np.pi ) )
                self.emission_angles.append( \
                    np.mod( \
                        np.interp( time, self.timestamps, self.emisangles ), np.pi ) )




class ExcitationMotor:
    """This class will hold the data associated with the excitation polarizer.
    It is used to:
        *   read in the file generated by the labview component
        *   extrapolate the angular function (assumed to be linear)
        *   present the function such that it can be queried
    """

    def __init__(self, filename, \
                     phase_offset_in_deg=0, \
                     rotation_direction=-1, \
                     optical_element='l/2 plate'):
        """Initialize the class: read in the file"""
        self.experiment_start_datetime = None
        self.filename = filename
        self.phase_offset_in_deg = phase_offset_in_deg
        self.rotation_direction = rotation_direction
        self.optical_element = optical_element
        
        self.load_excitation_motor_file( filename )
        self.determine_function()


    def load_excitation_motor_file( self, filename ):
        """Reads the labview file written for that motor. 
        The file currently gives a timestamp, followed by a signal which can take 
        the values START, UP, DOWN and END. 
        START marks the begin of the measurement, and recurring UP (or equivalently DOWN)
        mark complete revolutions of the lambda/2 plate. END denotes the end of the
        measurement.

        The function returns the time stamps with corresponding motor position.
        """

        # grab motor data --- delimiter is _a tab_ !!!
        # uses converter functions to parse date+time and shutter values
        md = np.loadtxt( filename, \
                             delimiter='\t', \
                             skiprows=1, \
                             converters={0: lambda s: deal_with_date_time_string(self, s), \
                                             1: lambda s: s=='UP'} )

        # timestamps in sec
        timestamps = md[:,0]
        # signal is 1 for UP events, 0 for everything else
        signals    = md[:,1]

        self.timestamps = timestamps
        self.signals    = signals
        self.starttime  = timestamps[0]
        self.endtime    = timestamps[-1]


    def angle(self, time, raw_angles=True):
        """This function will return the angle of the excitation polarization given
        a time (in seconds, counting from the start of the experiment).
        Remember, the angle of the polarization is _twice_ the angle of the polarizer.
        """
        if self.starttime <= time <= self.endtime:
            phi = self.anglefun_slope*time + self.anglefun_intercept + self.phase_offset_in_deg*np.pi/180.0
            if raw_angles:
                return phi
            else:
                return phi % np.pi
        else:
            return -1


    def determine_function(self):
        """ This is an interpolating (fitting!) function which is useful
        for the _excitation motor_.
        """
        a = self.signals[self.signals==1]
        t = self.timestamps[self.signals==1]
        if self.optical_element.lower()=='l/2 plate':
            # 720 because the polarization has rotated twice once the 
            # polarizer (the lambda/2) has rotated once
            a *= 4*np.pi*self.rotation_direction
        elif self.optical_element.lower()=='polarizer':
            a *= 2*np.pi*self.rotation_direction
        else:
            raise ValueError("Motor class doesn't know optical_element '%s'" % self.optical_element)
        a = np.cumsum(a)

        # make t into a 2d-array
        t = np.reshape( t, (t.shape[0],1) )

        # prep matrix for fit
        #   has form:  (t_1 1)
        #              (t_2 1)
        #              ( .  .)
        #              (t_n 1)
        M = np.concatenate( (t, np.ones_like(t)), axis=1 )

        s = np.linalg.lstsq( M, a )
        
        self.anglefun_slope     = s[0][0]
        self.anglefun_intercept = s[0][1]

        

class EmissionMotor:
    """This class will hold the data associated with the emission polarizer.
    It is used to:
        *   read in the file generated by the labview component
        *   present the angular function such that it can be queried,
            using something like 'emission_angle(time)'.
        """

    def __init__( self, filename ):
        """Initialize the class: read in the file"""
        self.experiment_start_datetime = None
        self.filename = filename
        self.load_emission_motor_file( filename )


    def load_emission_motor_file( self, filename ):
        """Reads the labview file written for that motor. 
        It assumes that fields (time, motor position(s), shutter) are separated 
        by _two_ spaces (which seems to be the case in all file generated so far).
        This makes it easy to read in (and parse) the date and time as a single 
        string.
        This function currently also parses the shutter value: open/closed is converted 
        into a boolean with an inline function doing string comparison.

        The function returns the time stamps with corresponding motor position.
        """

        # grab motor data --- delimiter is _a tab_ !!!
        # uses converter functions to parse date+time and shutter values
        md = np.loadtxt( filename, \
                             delimiter='\t', \
                             skiprows=1, \
                             converters={0: lambda s: deal_with_date_time_string(self, s), \
                                             3: lambda s: s=='open'} )
                            
        timestamps = md[:,0]
        emisangles = md[:,1]
#        exciangles = md[:,2]
        shutter    = md[:,3]

        self.timestamps = timestamps
        self.angles     = emisangles * np.pi/180.0
        self.shutter    = shutter


    def angle(self, time, exposuretime=.1, respectShutter=True ):
        """This function will return the angle of the emission polarizer given
        a time (in seconds, counting from the start of the experiment).
        Since the emission polarizer is turning in a very discrete fashion (0deg, 45deg, etc),
        we can assume that we can return the angle which is closest to the 
        the queried time (nearest interpolation), with the following caveats:
        if the any of the nearest-neighbor timestamps reports the shutter as closed, 
        then we return -1.
        """

        if respectShutter:
            ### find all timestamps within one exposure time ###
            mintime = time - exposuretime
            maxtime = time + exposuretime
            # left boundary of the window: move time axis by mintime
            tt = self.timestamps - mintime
            # last element which is less than zero must be 
            # the timestamp just before the window starts,
            # so the length of the array tt[tt>0] is the index
            # for that timestamp, and adding one to it we get
            # the first timestamp in the window. Since indices
            # start with zero, we subtract that one again and 
            # we have simply:
            first = tt[tt<0].size
            
            # right boundary,
            tt = self.timestamps - maxtime
            last = tt[tt<0].size -1

            # now test if the shutter is off at any of the timestamps
            # between first and last:
            if np.any( self.shutter[first:last]==0 ):
                return -1
            else:
                # interpolate angle from known angles
                phi = np.interp( time, self.timestamps, self.angles )
                return phi % np.pi
        else:
            phi = np.interp( time, self.timestamps, self.angles )
            return phi % np.pi

