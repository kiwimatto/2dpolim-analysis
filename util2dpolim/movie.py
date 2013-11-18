from cameradata import CameraData
from fitting import CosineFitter_new
from spot import Spot
#from portrait import Portrait
#from misc import pixel_list
from motors import *
import os, os.path, time, sys
import scipy.optimize as so
import multiprocessing

class Movie:
    def __init__( self, datadir, basename, phase_offset_in_deg=np.nan ):

        self.cos_fitter = CosineFitter_new
        self.phase_offset_in_deg = phase_offset_in_deg

        # format path correctly for whichever OS we're running under, and add trailing separator
        self.data_directory = os.path.normpath( datadir ) + os.path.sep
        self.data_basename  = basename
 
        self.sample_data = None
        self.blank_data  = None
        self.exspot_data = None

        self.bg_spot_sample    = None
        self.bg_spot_blank     = None
        self.bg_spot_exspot    = None

        self.read_in_EVERYTHING()

        self.spots = []
        self.validspots = None
        self.initContrastImages()

        # fix the excitation and emission angle grids for now
        TargetNAngles = 90
        # ExAngleGridSize = np.round( TargetNAngles/float(len(self.unique_exangles)) ) * len(self.unique_exangles)
        # EmAngleGridSize = np.round( TargetNAngles/float(len(self.unique_emangles)) ) * len(self.unique_emangles)
        ExAngleGridSize = 6
        EmAngleGridSize = 4
        phase_offset_in_rad = self.motors.phase_offset_in_deg * np.pi/180.0
        self.excitation_angles_grid = np.linspace(0,np.pi,ExAngleGridSize, endpoint=False)
        self.emission_angles_grid   = np.linspace(0,np.pi,EmAngleGridSize, endpoint=False)
        self.excitation_angles_grid += phase_offset_in_rad
        self.excitation_angles_grid = np.mod( self.excitation_angles_grid, np.pi )
#        print self.unique_exangles
#        print self.excitation_angles_grid

        self.uexa_portrait_indices = []
        for uexa in self.unique_exangles:
            assert np.any(np.abs(uexa-self.excitation_angles_grid) < 10*np.finfo(np.float).eps)
            self.uexa_portrait_indices.append( np.argmin( np.abs(uexa-self.excitation_angles_grid) ) )

        self.uema_portrait_indices = []
        for uema in self.unique_emangles:
            assert np.any(np.abs(uema-self.emission_angles_grid) < 10*np.finfo(np.float).eps)
            self.uema_portrait_indices.append( np.argmin( np.abs(uema-self.emission_angles_grid) ) )

        self.Nexphases_for_cos_fitter = 181
        self.Nemphases_for_cos_fitter = 181
        # self.precomputed_cosines = [ np.cos(2*(self.emission_angles_grid-ph)) \
        #                                  for ph in np.linspace(0,np.pi/2,self.Nphases_for_cos_fitter) ]


    def run_mp(self, Nprocs, fits=True, mods=True, ETruler=True):

        assert fits >= mods >= ETruler

        def worker( resultqueue, thesespots, whichproc ):
            if fits:
                self.fit_all_portraits_spot_parallel_selective( myspots=thesespots )
            if mods:
                self.find_modulation_depths_and_phases_selective( myspots=thesespots )
            if ETruler:                
                for si in thesespots:
                    self.validspots[si].values_for_ETruler( newdatalength=1024 )
            for si in thesespots:
                a = {'spot':si}
                if fits:
                    a['linefitparams']     = self.validspots[si].linefitparams
                    a['residual']          = self.validspots[si].residual
                    a['verticalfitparams'] = self.validspots[si].verticalfitparams
                if mods:
                    a['M_ex']       = self.validspots[si].M_ex
                    a['M_em']       = self.validspots[si].M_em
                    a['phase_ex']   = self.validspots[si].phase_ex
                    a['phase_em']   = self.validspots[si].phase_em
                    a['LS']         = self.validspots[si].LS
                    a['anisotropy'] = self.validspots[si].anisotropy
                if ETruler:
                    a['ET_ruler'] = self.validspots[si].ET_ruler

                resultqueue.put( a )

            return

        myspots = np.array_split( np.arange(len(self.validspots)), Nprocs )
        resultqueue = multiprocessing.Queue()
        jobs = [multiprocessing.Process(target=worker, args=(resultqueue, myspots[i], i)) for i in range(Nprocs)]
        for job in jobs: job.start()

        myspots = list(np.concatenate(myspots))
        while len(myspots)>0:
            a = resultqueue.get()   # this will block until it gets something
            si = a.pop('spot')
            # print si
            # print myspots
            for key in a.keys():
                setattr( self.validspots[si], key, a[key] )
            i = myspots.index(si)
            myspots.pop(i)

        for job in jobs: job.join()

        assert resultqueue.empty()

        self.update_images()

        # while p.is_alive()==True:
        #     for i in range(Nprocs):
        #         print 'i=',i,' alive --> ',jobs[i].is_alive()
        #     time.sleep(.1)

        # this should block until all processes have finished
        
        # print self.validspots[si].M_ex
        # print 'p:',self.validspots[si].pixel[0][0],',',self.validspots[si].pixel[0][1]
        # print '.',self.M_ex_image[ self.validspots[si].pixel[0][0], self.validspots[si].pixel[0][1] ]
        
        print 'all done i guess'

    def update_images(self):
        allprops = ['M_ex', 'M_em', 'phase_ex', 'phase_em', 'LS', 'anisotropy', 'ET_ruler']
        for s in self.validspots:
            for prop in allprops:
                if hasattr(s,prop):
                    self.store_property_in_image( s, prop+'_image', prop )


    def read_in_EVERYTHING(self):
        # change to the data directory
        os.chdir( self.data_directory )
        
        ###### look if we can find the data file, and import it ######
        print 'Looking for data file: %s.[spe|SPE]... ' % (self.data_directory+self.data_basename)

        if os.path.exists( self.data_basename+'.SPE' ):
            self.spefilename = self.data_basename+'.SPE'
        elif os.path.exists( self.data_basename+'.spe' ):
            self.spefilename = self.data_basename+'.spe'
        elif os.path.exists( self.data_basename+'.npy' ):     # testing data is captured here!
            self.spefilename = self.data_basename+'.npy'
        else:
            raise IOError("Couldn't find data SPE file! Bombing out...")
        
        self.sample_data    = CameraData( self.spefilename, compute_frame_average=True )
        self.timeaxis       = self.sample_data.timestamps
        print 'Imported data file %s' % self.spefilename


        ###### look for motor files ######
        print 'Looking for motor file(s)...',
        
        got_motors = False
        for file in os.listdir("."):
            if file=="MS-"+self.data_basename+'.txt':
                print '\t found motor file %s' % file
                self.motorfile = file
                got_motors = True
                break

        if got_motors:
            ### Supplying the phase offset here means it overrides the value 
            ### which is read from the header, unless it is NaN (the default). 
            ### This is only useful for AM measurements.
            self.motors = BothMotorsWithHeader( self.motorfile, self.phase_offset_in_deg )
            self.exangles = self.motors.excitation_angles
            self.emangles = self.motors.emission_angles
            self.unique_exangles = np.unique( self.exangles )
            self.unique_emangles = np.unique( self.emangles )
        else:
            raise IOError("Couldn't find motor files.")
       
        ###### look for blank sample ######
        print 'Looking for blank...',
        got_blank = False
        for file in os.listdir("."):
            if file.startswith("blank-") and (file.endswith(".spe") or file.endswith(".SPE")):
                print '\t found file %s' % file
                self.blank_data = CameraData( file, compute_frame_average=True )
                got_blank = True
                break
        if not got_blank: print ' none.'

        ###### look for excitation spot sample ######
        print 'Looking for excitation spot data...',
        got_exspot = False
        for file in os.listdir("."):
            if file.startswith("excitation-spot-") and (file.endswith(".spe") or file.endswith(".SPE")):
                print '\t found file %s' % file
                self.exspot_data = CameraData( file, compute_frame_average=True )
                # we also collapse all frames (if more than one) into a single mean frame:
                self.exspot_data.rawdata = np.mean( self.exspot_data.rawdata, axis=0 )
                self.exspot_data.datasize[0] = 1
                self.exspot_data.timeaxis    = self.exspot_data.timeaxis[0]
                got_exspot = True
                break
        if not got_exspot: print ' none.'


    def initContrastImages(self):
        self.spot_coverage_image  = np.ones( (self.sample_data.datasize[1],self.sample_data.datasize[2]) )*np.nan
        self.original_mean_intensity_image = np.mean( self.sample_data.rawdata, axis=0 )
        self.mean_intensity_image = self.spot_coverage_image.copy()
        self.meanSNR_image        = self.spot_coverage_image.copy()
        self.framecountSNR_image  = self.spot_coverage_image.copy()
        self.M_ex_image           = self.spot_coverage_image.copy()
        self.M_em_image           = self.spot_coverage_image.copy()
        self.phase_ex_image       = self.spot_coverage_image.copy()
        self.phase_em_image       = self.spot_coverage_image.copy()
        self.LS_image             = self.spot_coverage_image.copy()
        self.anisotropy_image     = self.spot_coverage_image.copy()
        self.ET_ruler_image       = self.spot_coverage_image.copy()
        self.ET_model_md_fu_image = self.spot_coverage_image.copy()
        self.ET_model_th_fu_image = self.spot_coverage_image.copy()
        self.ET_model_gr_image    = self.spot_coverage_image.copy()
        self.ET_model_et_image    = self.spot_coverage_image.copy()


    def define_background_spot( self, shape, intensity_type='mean' ):
        
        if not type(shape)==dict:
            shape = {'type': 'Rectangle', \
                         'upper': shape[1], \
                         'lower': shape[3], \
                         'left': shape[0], \
                         'right': shape[2] }

        # create new spot object
        s = Spot( shape, int_type=intensity_type, label='sample bg area', parent=self, \
                      is_bg_spot=True, which_bg='sample' )
        self.bg_spot_sample = s

        # if blank data is present, automatically create a bg spot here too
        if not self.blank_data==None:
            s = Spot( shape, int_type=intensity_type, label='blank bg area', parent=self, \
                          is_bg_spot=True, which_bg='blank' )
            self.bg_spot_blank = s

        # if excitation spot data is present, create a bg spot here too
        if not self.exspot_data==None:
            s = Spot( shape, int_type=intensity_type, label='exspot bg area', parent=self, \
                          is_bg_spot=True, which_bg='exspot' )
            self.bg_spot_exspot = s

        # record background spot in spot coverage image
        for p in s.pixel:
            self.spot_coverage_image[ p[0], p[1] ] = -1
        

    def define_spot( self, shape, intensity_type='mean', label=None, \
                         use_blank=True, use_exspot=False, use_borderbg=False ):
        """Defines a new spot object and adds it to the list.
        FIXME: make sure coordinate definitions do not exceed frame size
        """

        if not type(shape)==dict:
            shape = {'type': 'Rectangle', \
                         'upper': shape[1], \
                         'lower': shape[3], \
                         'left': shape[0], \
                         'right': shape[2] }
        # create new spot object
        s = Spot( shape, int_type=intensity_type, label=label, parent=self, \
                      is_bg_spot=False, which_bg='sample', \
                      use_blank=use_blank, use_exspot=use_exspot, use_borderbg=use_borderbg )
        # append spot object to spots list
        self.spots.append( s )

        # store in image
        self.store_property_in_image( s, 'spot_coverage_image', 'value_for_coverage_image' )

        return s

    def fit_blank_image( self, boolimage, verbosity=0 ):
        # The idea here is to fit the blank image to the sample image, but, of course,
        # only in region(s) where there's no interference from the actual sample.
        # For this we'll need a function which can return the list of pixel occupied by 
        # a rectangular shape, or, occupied by the border of a rectangular shape:

        collapsed_blank_data = np.mean( self.blank_data.rawdata, axis=0 )

        # boolimage = np.ones_like( collapsed_blank_data, dtype=np.bool )*True
        # rect = {'left':150, 'right':450, 'upper':200, 'lower':500, 'op':'exclude'}
        # boolimage = pixel_list( rect, boolimage )

        blankint  = collapsed_blank_data[boolimage]
        fitmatrix = np.vstack( [np.ones_like(blankint), blankint] ).T        
        fitblankimg = np.zeros_like( self.sample_data.rawdata )

        if verbosity>0:
            import time, os
            import matplotlib.pyplot as plt
            plt.interactive(True)
            fig = plt.figure( figsize=(16,8) )
            ax1 = fig.add_subplot(1,2,1)
            ax2 = fig.add_subplot(1,2,2)

#        print self.sample_data.datasize
#        print self.blank_data.datasize

        for fi,frame in enumerate(self.sample_data.rawdata):
            sampleint = frame[boolimage]
            res = np.linalg.lstsq( fitmatrix, sampleint.reshape((sampleint.size,1)) )
            fitblankimg[fi,:,:] = res[0][0] + res[0][1]*collapsed_blank_data
            if verbosity>0:
                print res[0]
                ax1.imshow( frame-fitblankimg[fi,:,:], interpolation='nearest', vmin=-20, vmax=100 )
                n,b = np.histogram( (frame-fitblankimg[fi,:,:]).flatten(), bins=np.linspace(-20,100,121) )
                ax2.bar( (b[:-1]+b[1:])/2.0, np.log10(n), width=np.diff(b)[0], alpha=.6 )
                ax2.set_xlim(-20,100)
                plt.draw()
                plt.savefig( 'blaaa_%03d.png' % fi )
                ax1.cla()
                ax2.cla()

        self.blank_data.rawdata = fitblankimg

        if verbosity>0:
            os.system('convert blaaa_*.png blaaa.gif')
            abc = np.sum( np.abs( self.sample_data.rawdata-fitblankimg ), axis=0 )
            plt.figure( figsize=(12,12) )
            plt.imshow( abc, interpolation='nearest' )
            plt.colorbar()

        # plt.figure( figsize=(12,12) )
        # plt.imshow( np.sum( self.sample_data.rawdata, axis=0 ), interpolation='nearest' )
        # plt.colorbar()



#        blankint = [ self.blank_data.rawdata[ p[0], p[1] ] for p in pixel ]
        


    def correct_excitation_intensities( self ):
        Icorr = self.motors.sample_plane_intensities/np.max(self.motors.sample_plane_intensities)
        for i,s in enumerate(self.spots):
            s.intensity /= Icorr


    def correct_emission_intensities( self ): #, corrM, corrphase ):

        corrM     = self.motors.header['em correction modulation depth']
        corrphase = self.motors.header['em correction phase']/180.0*np.pi

        if np.isnan(corrM): corrM=0
        if np.isnan(corrphase): corrM=phase
        
        # correction function
        corrfun = lambda angle: (1+corrM*np.cos(2*(angle + corrphase) ))/(1+corrM)

        # assemble all corrections in a vector
        corrections = np.nan*np.ones((self.emangles.size,))
        for i,emangle in enumerate(self.emangles):
            corrections[i] = corrfun( emangle )

#        print corrections

        # now go through all spots and apply that vector
        for i,s in enumerate(self.spots):
            s.intensity /= corrections


    def find_portraits( self, frameoffset=0 ):
        exangles_rounded = np.round( self.exangles, decimals=2 )
        emangles_rounded = np.round( self.emangles, decimals=2 )
        number_of_verticals = np.unique( exangles_rounded ).size
        number_of_lines     = np.unique( emangles_rounded ).size

        Nframes    = self.motors.excitation_angles.size
        fpp        = number_of_verticals * number_of_lines      # frames per portrait
        Nportraits = Nframes/fpp

        assert np.mod( Nframes, fpp )==0    # sanity check
        
        portrait_indices = []
        for pi in range(Nportraits+1):
            ind   = pi*fpp + frameoffset
            if not ind > Nframes:
                portrait_indices.append( ind )
            else:
                break

        Nportraits = len(portrait_indices)-1
#        print Nportraits
#        print np.array(portrait_indices)

        self.Nportraits = Nportraits
        self.portrait_indices = portrait_indices


    def find_lines( self ):
        all_line_indices = []
        for pi in range(self.Nportraits):
            line_indices = []
            # get angles in this portrait        
            pstart = self.portrait_indices[ pi ]
            pstop  = self.portrait_indices[ pi+1 ]
            emangles_rounded = np.round( self.emangles[ pstart:pstop ], decimals=2 )
            unique_emangles  = np.unique( emangles_rounded )
            number_of_lines  = unique_emangles.size
            for li in range(number_of_lines):
                line_indices.append( unique_emangles[li]==emangles_rounded )
            all_line_indices.append( line_indices )
#        print all_line_indices
        
        # for i in range(self.Nlines):
        #     print self.motors.emission_angles[ pstart:pstop ][ all_line_indices[0][i] ]

        self.line_indices = all_line_indices
        self.Nlines = number_of_lines


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

        # collect list of unique emission angles (same for all spots!)                    
        emangles = [l.emangle for l in portrait.lines]
        # turn into array, transpose and squeeze
        emangles = np.squeeze(np.array( emangles ).T)

        # evaluate cosine-fit at these em_angles, on a grid of ex_angles:
        fitintensities = np.array([l.cosValue( self.excitation_angles_grid ) for l in portrait.lines])

        phase, I0, M, resi, fit, rawfitpars, mm = self.cos_fitter( emangles, fitintensities, \
                                                                       self.Nemphases_for_cos_fitter )         
        # phasor addition!
        proj_em = np.real( rawfitpars[0,:] * np.exp( 1j*rawfitpars[1,:] ) \
                               * np.exp( 1j*self.emission_angles_grid ) )

        phase, I0, M, resi, fit, rawfitpars, mm = self.cos_fitter( self.emission_angles_grid, proj_em, \
                                                                       self.Nphases_for_cos_fitter ) 
        portrait.phase_em = phase[0]
        portrait.M_em = M
        

    def retrieve_angles( self, iportrait, iline ):
        pstart = self.portrait_indices[ iportrait ]
        pstop  = self.portrait_indices[ iportrait+1 ]
        # lstart = self.line_edges[ iline ]
        # lstop  = self.line_edges[ iline+1 ]
        # exangles = self.exangles[ pstart:pstop ][ lstart:lstop ]
        # emangle  = self.emangles[ pstart:pstop ][ lstart:lstop ]
        exangles = self.exangles[ pstart:pstop ][ self.line_indices[iportrait][iline] ]
        emangle  = self.emangles[ pstart:pstop ][ self.line_indices[iportrait][iline] ]
        # print '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
        # print exangles
        # print self.exangles
        # print emangle
        # print self.emangles
        # print '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
        assert np.all(emangle[0]==emangle)  # emangle must be constant with a scan line
        return exangles, emangle[0]


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
            pstart = self.portrait_indices[ pi ]
            pstop  = self.portrait_indices[ pi+1 ]
            
            # part I, 'horizontal fitting' of the lines of constant emission angles

            # for each line ---- we do lines in series, __but all spots in parallel__:
            for li in range(self.Nlines):

                # start and stop indices for this line
#                lstart = self.line_edges[ li ]
#                lstop  = self.line_edges[ li+1 ]

                # get excitation angle array (same for all spots!)
#                exangles = self.exangles[ pstart:pstop ][ lstart:lstop ]
                exangles = self.exangles[ pstart:pstop ][ self.line_indices[pi][li] ]

                # create list of intensity arrays (one array for each spot)
                intensities = [self.validspots[si].retrieve_intensity(iportrait=pi, iline=li) for si in myspots]

                # turn into numpy array and transpose
                intensities = np.array( intensities ).T

                exa = exangles.copy()

                phase, I0, M, resi, fit, rawfitpars, mm = self.cos_fitter( exa, intensities, \
                                                                               self.Nexphases_for_cos_fitter ) 

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
#            emangles = self.emangles[pstart:pstop][self.line_edges[:-1]]
            emangles = []
            for li in range(self.Nlines):
                emangles.append( self.emangles[pstart:pstop][ self.line_indices[pi][li] ][0] )
#            print emangles

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
                                                                           self.Nemphases_for_cos_fitter ) 
                
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

                # print 'vertical fit params:'
                # print self.validspots[si].verticalfitparams[pi,:,:]

                # print 'line fit params:'
                # print self.validspots[si].linefitparams[pi,:,:]

        return self.validspots


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
            self.cos_fitter( self.excitation_angles_grid, proj_ex, self.Nexphases_for_cos_fitter )
        ph_em, I_em, M_em, r_em, fit_em, rawfitpars_em, mm = \
            self.cos_fitter( self.emission_angles_grid, proj_em, self.Nemphases_for_cos_fitter )

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
            self.store_property_in_image( self.validspots[si], 'M_ex_image', 'M_ex' )
            self.store_property_in_image( self.validspots[si], 'M_em_image', 'M_em' )
            self.store_property_in_image( self.validspots[si], 'phase_ex_image', 'phase_ex' )
            self.store_property_in_image( self.validspots[si], 'phase_em_image', 'phase_em' )
            self.store_property_in_image( self.validspots[si], 'LS_image', 'LS' )
            # a = self.validspots[si].coords[1]
            # b = self.validspots[si].coords[3]+1
            # c = self.validspots[si].coords[0]
            # d = self.validspots[si].coords[2]+1
            # self.M_ex_image[ a:b, c:d ]     = self.validspots[si].M_ex
            # self.M_em_image[ a:b, c:d ]     = self.validspots[si].M_em
            # self.phase_ex_image[ a:b, c:d ] = self.validspots[si].phase_ex
            # self.phase_em_image[ a:b, c:d ] = self.validspots[si].phase_em
            # self.LS_image[ a:b, c:d ]       = self.validspots[si].LS        

            #### anisotropy ####
        
            if (self.validspots[si].M_ex < .015):
                iex = np.argmin( np.abs( self.excitation_angles_grid - np.pi/2 ) )
                iem = np.argmin( np.abs( self.emission_angles_grid - np.pi/2 ) )
                # Ipara = self.validspots[si].sam[ iex, iem ]
                # Iperp = self.validspots[si].sam[ iex, 0 ] ## assumes that the first index is close to 0deg in emission angle grid
                Ipara = self.validspots[si].sam[ iem, iex ]
                Iperp = self.validspots[si].sam[ 0, iex ] ## assumes that the first index is close to 0deg in emission angle grid
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
                Ipara = self.validspots[si].sam[ iphempara, iphex ]
                Iperp = self.validspots[si].sam[ iphemperp, iphex ]
                
                # self.validspots[si].sam[ iex, iem ] = 10000
                # self.validspots[si].sam[ iex, 0 ] = 10000
                # self.validspots[si].sam[ iphex, iphempara ] = 7000
                # self.validspots[si].sam[ iphex, iphemperp ] = 4000
                # import matplotlib.pyplot as plt
                # plt.imshow( self.validspots[si].sam, origin='lower', interpolation='nearest' )
                # plt.show()
                # raise SystemExit

            if not float(Ipara+2*Iperp)==0:
                self.validspots[si].anisotropy = float(Ipara-Iperp)/float(Ipara+2*Iperp)
            else:
                self.validspots[si].anisotropy = np.nan
            # store in contrast image
            self.store_property_in_image( self.validspots[si], 'anisotropy_image', 'anisotropy' )
#            self.r_image[ a:b, c:d ] = self.validspots[si].r
#            del(self.validspots[si].sam)  # don't need that anymore
        si=myspots[0]
#        print self.validspots[si].M_ex
#        print 'p:',self.validspots[si].pixel[0][0],',',self.validspots[si].pixel[0][1]
#        print '.',self.M_ex_image[ self.validspots[si].pixel[0][0], self.validspots[si].pixel[0][1] ]

#        return self.validspots


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
            self.store_property_in_image( self.validspots[si], 'ET_ruler_image', 'ET_ruler' )
            

    def ETmodel_selective( self, myspots, fac=1e4, pg=1e-9, epsi=1e-11 ):

        def look_up_cosine( self, spot, exa, emas ):
            mycos = lambda a, ph, I, M: I*( 1+M*( np.cos(2*(a-ph)) ) )

            which_portrait = 0
            iexa = np.argmin( np.abs(exa-self.excitation_angles_grid) )
            assert np.abs(exa-self.excitation_angles_grid[iexa]) < 1e-5

            return mycos( emas, \
                              spot.verticalfitparams[which_portrait,iexa,0], \
                              spot.verticalfitparams[which_portrait,iexa,1], \
                              spot.verticalfitparams[which_portrait,iexa,2] )

        from fitting import fit_portrait_single_funnel_symmetric

        for si in myspots:
            print 'ETmodel fitting spot %d' % si

            # we 'correct' the modulation in excitation to be within 
            # limits of reason (and proper arccos functionality)
            mex = np.clip( self.validspots[si].M_ex, .000001, .999999 )

            a0 = [mex, 0, 1]
            EX, EM = np.meshgrid( self.excitation_angles_grid, self.emission_angles_grid )
            funargs = (EX, EM, self.validspots[si].sam, mex, self.validspots[si].phase_ex, 'fitting')
            # EX = self.excitation_angles_grid
            # EM = np.linspace( 0, 180, 20, endpoint=False )/180*np.pi
            # spotint = np.array([look_up_cosine(self, self.validspots[si], exa, EM) for exa in EX]).T
            # EX, EM = np.meshgrid( EX, EM )
            # funargs = (EX, EM, spotint, mex, self.validspots[si].phase_ex, 'fitting')

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
            # et,A = fit_portrait_single_funnel_symmetric( fitresult[0], \
            #                                                  EX, EM, spotint, mex, \
            #                                                  self.validspots[si].phase_ex, \
            #                                                  mode='show_et_and_A', use_least_sq=True )
            self.validspots[si].ETmodel_md_fu = fitresult[0][0]
            self.validspots[si].ETmodel_th_fu = fitresult[0][1]
            self.validspots[si].ETmodel_gr    = fitresult[0][2]
            self.validspots[si].ETmodel_et    = et            
            # a = self.validspots[si].coords[1]
            # b = self.validspots[si].coords[3]+1
            # c = self.validspots[si].coords[0]
            # d = self.validspots[si].coords[2]+1
            # self.ET_model_md_fu_image[ a:b, c:d ] = fitresult[0][0]
            # self.ET_model_th_fu_image[ a:b, c:d ] = fitresult[0][1]
            # self.ET_model_gr_image[ a:b, c:d ]    = fitresult[0][2]
            # self.ET_model_et_image[ a:b, c:d ]    = et

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


    def are_spots_valid(self, SNR=10, validframesratio=.7, quiet=False, all_are_valid=False):

        # not sure if we have a value for the bg std 
        bgstd=None

        if not self.bg_spot_sample==None:         # yes we do have a bg spot
            bgstd = self.bg_spot_sample.std

        # create a new list containing valid spots only
        validspots = []   
        validspotindices = []
        for si,s in enumerate(self.spots):

            if s.use_borderbg:        # hang on a sec, spot knows its own bg!
                bgstd = s.borderbgstd

            if not bgstd==None:
                s.SNR = np.abs( s.intensity/bgstd )     # --> SNR for each frame
            else:
                s.SNR = np.ones_like(s.intensity)*np.inf    # --> default inf SNR if bg unknown

            s.framecountSNR = np.sum(s.SNR > SNR)
            if s.framecountSNR >= validframesratio*len(s.SNR):
                validspots.append(s)
                validspotindices.append(si)
            s.meanSNR = np.mean(s.SNR)

            # store mean SNR in SNR_image
            self.store_property_in_image( s, 'meanSNR_image', 'meanSNR' )    # let's see if this works...
            self.store_property_in_image( s, 'framecountSNR_image', 'framecountSNR' )    # let's see if this works...

        # and store in movie object
        self.validspots = validspots
        self.validspotindices = validspotindices        

        for s in self.validspots:
            self.store_property_in_image( s, 'spot_coverage_image', 'value_for_coverage_valid_image' )

        if not quiet: print "Got %d valid spots (of %d spots total)" % (len(self.validspots),len(self.spots))

    def store_property_in_image(self, spot, image, prop):
        # record the spot in the mean intensity image
        for p in spot.pixel:
            getattr(self,image)[ p[0], p[1] ] = getattr(spot, prop)
#        image[ self.coords[1]:self.coords[3]+1, self.coords[0]:self.coords[2]+1 ] = getattr(self,prop)

