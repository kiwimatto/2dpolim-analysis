import numpy as np
from line import Line

class Portrait:
    def __init__( self, exangles, emangles, intensities, parent ):
        self.exangles = exangles
        self.emangles = emangles
        self.intensities = intensities
        self.parent = parent
        self.lines = []

        self.split_into_emission_lines()
        del(self.intensities)

    def split_into_emission_lines( self ):
        """Splits the portrait into the flat emission parts."""
    
        # edge 'detection' to work out where emission angles change
        edges = (np.diff(self.emangles)!=0).nonzero()[0]
        # extend the array: a zero on the left because we know that the
        # first frame is an edge, and a p.shape[0] on the right because
        # we know that the last frame is an edge as well.
        edges = np.concatenate( (np.array( [0] ), edges+1, np.array([self.exangles.shape[0]]) ) )

        print '-----------'
        print edges
        raise SystemExit
        
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

        # did we lose any data rows due to our funky indexing?
        assert np.sum( [lines[i].exangles.shape[0] for i in range(len(lines))] ) == self.exangles.shape[0]
    
        self.lines = lines

    def show_portrait_matrix( self ):
        plt.matshow( self.matrix, origin='bottom')
        plt.plot( [0,180], [0,180], 'k-' )
        plt.xlim( 0, 180 )
        plt.ylim( 0, 180 )


    def recover_portrait_matrix(self):
        mycos = lambda a, ph, I, M: I*( 1+M*( np.cos(2*(a-ph)) ) )

        phase = self.vertical_fit_params[0]
        I0    = self.vertical_fit_params[1]
        M     = self.vertical_fit_params[2]

        pic = np.zeros( (self.parent.parent.emission_angles_grid.size, self.parent.parent.excitation_angles_grid.size) )
        for exi in range( self.parent.parent.excitation_angles_grid.size ):
            pic[:,exi] = mycos( self.parent.parent.emission_angles_grid, phase[exi], I0[exi], M[exi] )

        return pic
