import numpy as np
from spefiles import MyPrincetonSPEFile

class CameraData:

    def __init__( self, spe_filename, compute_frame_average=False, in_counts_per_sec=True ):
        # load SPE  ---- this will work for SPE format version 2.5 (probably not for 3...)

        self.filename           = spe_filename

        if self.filename.split('.')[-1]=='npy':   # we got test data, presumably
            print "======== TEST DATA IT SEEMS =========="
            self.rawdata      = np.load(self.filename)
            self.datasize     = self.rawdata.shape
            self.exposuretime = .1    # in seconds

        else:                                     # we got real data 
            self.rawdata_fileobject = MyPrincetonSPEFile( self.filename )
            self.rawdata            = self.rawdata_fileobject.return_Array()#.astype(np.float64)
            self.datasize           = self.rawdata_fileobject.getSize()
            self.exposuretime       = self.rawdata_fileobject.Exposure   # in seconds
            # scale signal to counts/second:
            if in_counts_per_sec:
                self.rawdata           /= self.rawdata_fileobject.Exposure
            self.rawdata_fileobject.close_file()
            del(self.rawdata_fileobject)
        if compute_frame_average:
            self.average_image      = np.mean( self.rawdata, axis=0 )

        ###  extract or generate time stamps ###
        #  here we do not have timestamps for each frame, so we
        #  generate time axis from SPE exposure data and number of frames
        self.timestamps = np.linspace( 0, self.exposuretime*self.rawdata.shape[0], \
                                           self.rawdata.shape[0], endpoint=False )

