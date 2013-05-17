import numpy as np
from util_misc import deal_with_date_time_string


class NewSetupMotor:
    """This class represents either of the two motors in the new
    setup. For both, emission and excitation motor, this class should
    be used. There are only two reasons we have separate motor classes
    for the old setup: i) because the data generated by labview is a
    separate, differently formatted file for the excitation motor, and
    ii) the polarization phase offset is implemented for the
    excitation motor.

    This new class is _very_ similar to the emission motor class
    below. In fact, it is worth considering to consolidate the three
    different classes into one, or at least to use this generic class
    for the emission motor in the old setup, and have only one extra
    class specialized on the old excitation motor.

    (Alternatively, we could create separate branches for the old and
    new setup, however i think this is unnecessary, since the
    differences in the analysis software are limited to motor data
    import --- so far.)

    Class constructor arguments are the filename to be read in, a flag
    denoting which_motor should be read in (either 'emission' or
    'excitation'), as well as a phase offset, which, by convention,
    will always be zero for the emission motor, but can take non-zero
    values for the excitation motor.
    """

    def __init__( self, filename, which_motor, phase_offset=0 ):
        """Initialize the class: read in the file"""
        self.experiment_start_datetime = None
        self.filename = filename
        self.phase_offset = phase_offset
        self.load_motor_file( filename, which_motor )


    def load_motor_file( self, filename, which_motor ):
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
        exciangles = md[:,2]*2
        shutter    = md[:,3]

        self.timestamps = timestamps
        if which_motor=='excitation':
            self.angles     = exciangles
        elif which_motor=='emission':
            self.angles     = emisangles
        else:
            raise ValueError("Input argument which_motor to class NewSetupMotor must take values 'excitation' or 'emission'. Got: %s" % (which_motor))
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
                phi = np.interp( time, self.timestamps, self.angles ) + self.phase_offset
                return phi % 180.0
        else:
            phi = np.interp( time, self.timestamps, self.angles ) + self.phase_offset
            return phi % 180.0

    



class ExcitationMotor:
    """This class will hold the data associated with the excitation polarizer.
    It is used to:
        *   read in the file generated by the labview component
        *   extrapolate the angular function (assumed to be linear)
        *   present the function such that it can be queried
    """

    def __init__(self, filename, phase_offset_excitation=0, rotation_direction=-1 ):
        """Initialize the class: read in the file"""
        self.experiment_start_datetime = None
        self.filename = filename
        self.phase_offset_excitation = phase_offset_excitation
        self.rotation_direction = rotation_direction

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


    def angle(self, time, do_modulus=True):
        """This function will return the angle of the excitation polarization given
        a time (in seconds, counting from the start of the experiment).
        Remember, the angle of the polarization is _twice_ the angle of the polarizer.
        """
        if self.starttime <= time <= self.endtime:
            phi = self.anglefun_slope*time + self.anglefun_intercept + self.phase_offset_excitation
            if do_modulus:
                return phi % 180.0
            else:
                return phi
        else:
            return -1


    def determine_function(self):
        """ This is an interpolating (fitting!) function which is useful
        for the _excitation motor_.
        """
        a = self.signals[self.signals==1]
        t = self.timestamps[self.signals==1]
        # 720 because the polarization has rotated twice once the 
        # polarizer (the lambda/2) has rotated once
        a *= 720*self.rotation_direction
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
        self.angles     = emisangles
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
                return phi % 180.0
        else:
            phi = np.interp( time, self.timestamps, self.angles )
            return phi % 180.0

