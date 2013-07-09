import numpy as np

class Line:
    def __init__( self, exangles, intensities, emangle ):
        self.exangles = exangles
        self.intensities = intensities
        self.emangle = emangle
        self.was_fitted = False

    def set_fit_params( self, phase, I0, M0, residuals ):
        self.was_fitted = True
        self.phase = phase
        self.I0 = I0
        self.M0 = M0
        self.resi = residuals
        
    def cosValue( self, angle ):
        return self.I0 * ( 1+self.M0*np.cos( 2*(angle-self.phase) ) )


