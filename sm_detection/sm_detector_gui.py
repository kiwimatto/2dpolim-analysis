# coding=utf-8

from PyQt4 import QtCore,QtGui
import sys, os, time
import numpy as np
import matplotlib.cm as cm
import sm_detector_layout
from spefiles import MyPrincetonSPEFile
import smdetect

class sm_detector_gui_logic(QtGui.QMainWindow,sm_detector_layout.Ui_MainWindow):

    def __init__(self,parent=None,app=None):
        """
            Initialization of the class. Call the __init__ for the super classes
        """
        super(sm_detector_gui_logic,self).__init__(parent)
        self.setupUi(self)
        self.connectActions()

        self.pwd = os.path.dirname(os.path.abspath(__file__))
        self.app = app

        self.filedir = self.pwd   # to be changed to the dir containing the last selected file

        self.imageview.clear()
#        self.axes1 = self.imageview.fig.add_subplot(111)
        self.maskview.clear()
#        self.axes2 = self.maskview.fig.add_subplot(111)

        self.defaultfile = None
        self.hi_my_name_is = "the sm detector gui logic class!"

    def main(self):
        self.show()

    def connectActions(self):
        """Connect the user interface controls to the logic """
        self.loadFilePushButton.clicked.connect( self.loadFile )
        self.applyFilterPushButton.clicked.connect( self.applyFourierFilter )
        self.findClustersPushButton.clicked.connect( self.findClusters )
        self.frameSlider.valueChanged.connect( self.drawFrame )
        self.frameAverageToggleButton.clicked.connect( self.frameAverageToggle )
        self.showWhatComboBox.activated.connect( self.drawFrame )
        self.highlightClusterPositionsCheckBox.stateChanged.connect( self.drawFrame )
        self.filterAndClusterAllFramesPushButton.clicked.connect( self.filter_and_cluster_all_frames )
        self.cropPushButton.clicked.connect( self.crop )
        self.resetPushButton.clicked.connect( self.reset )
        self.loadDefaultPushButton.clicked.connect( self.loadDefault )
        self.clusterIntensityThresholdSlider.valueChanged.connect( self.clusterIntensityThresholdChanged )
        self.clusterSizeThresholdSlider.valueChanged.connect( self.clusterSizeThresholdChanged )
        self.createLineAveragePushButton.clicked.connect( self.createLineAverages )
        self.toolBox.currentChanged.connect( self.toolBoxPageChanged )
        self.applyCorrelatorPushButton.clicked.connect( self.correlateFrames )
        self.addFitExclusionZonePushButton.clicked.connect( self.addFitExclusionZone )
        self.resetFitExclusionZonePushButton.clicked.connect( self.resetFitExclusionZone )
        self.blankFitPushButton.clicked.connect( self.blankFit )
        self.showBoolimageCheckBox.stateChanged.connect( self.drawFrame )

    def addFitExclusionZone(self):
        rect = {'left':np.int( np.round(self.imageview.c0) ), \
                    'right':np.int( np.round(self.imageview.c1) ), \
                    'upper':np.int( np.round(self.imageview.r0) ), \
                    'lower':np.int( np.round(self.imageview.r1) ), \
                    'op':'exclude'}
        self.pixel_list( rect )        
        self.drawFrame()

    def resetFitExclusionZone(self):
        self.boolimage = np.ones((self.blankfile.shape[1],self.blankfile.shape[2]), dtype=np.bool)*True        
        self.drawFrame()

    def blankFit(self):
        collapsed_blank_data = np.mean( self.blankfile, axis=0 )
        blankint  = collapsed_blank_data[self.boolimage]
        fitmatrix = np.vstack( [np.ones_like(blankint), blankint] ).T        
        fitblankimg = np.zeros_like( self.spefile )
        for fi,frame in enumerate(self.spefile):
            sampleint = frame[self.boolimage]
            res = np.linalg.lstsq( fitmatrix, sampleint.reshape((sampleint.size,1)) )
            fitblankimg[fi,:,:] = res[0][0] + res[0][1]*collapsed_blank_data

        self.spefile -= fitblankimg
        self.drawFrame()

    def pixel_list(self, shape):
        """ what follows below is mostly copypasta from the function spot.create_pixel_list """
        left  = shape['left']
        right = shape['right']
        lower = shape['lower']
        upper = shape['upper']
        # validate coords:
        assert lower >= upper   # lower in the image means larger row number!
        assert (left >= 0) and (upper >= 0)
        assert (right<self.spefile.shape[2]) and (lower<self.spefile.shape[1])

        for col in range(left,right+1):
            for row in range(upper,lower+1):
                if shape['op']=='include':
                    self.boolimage[ row, col ] = True
                elif shape['op']=='exclude':
                    self.boolimage[ row, col ] = False
                else:
                    raise ValueError("fit_blank_image(): can't comprehend op '%s'" % shape['op'])


    def toolBoxPageChanged(self):
        i = self.toolBox.currentIndex()
        for j in range(self.toolBox.count()):
            if not i==j:
                text = self.toolBox.itemText(j)
                self.toolBox.setItemText(j, '>>'+text[2:])
            else:
                text = self.toolBox.itemText(j)
                self.toolBox.setItemText(j, 'vv'+text[2:])

    def loadDefault(self):
        f = open('lastdir.txt','r')
        self.defaultfile = f.readline()
        f.close()

#        self.defaultfile = "/home/kiwimatto/Desktop/130925 - MEHPPV YUXI/TDM5/TDM5-488-OD1-01.SPE"
#        print "parent: ",self.parent()
        self.loadFile()

    def loadFile(self):
        if self.defaultfile==None:
            fname = QtGui.QFileDialog.getOpenFileName(self, 'Pick SPE file', \
                                                          directory=self.filedir, \
                                                          filter='SPE files (*.SPE *.spe)')
            fname = str(fname)
        else:
            fname = self.defaultfile

        if not fname=='':
            self.SPEFileLabel.setText(fname)
            self.SPEFileName = fname
            self.filedir = os.path.split( os.path.normpath(fname) )[0]
        else:
            return

        data_directory = os.path.normpath( str(fname) ) 

        # save in lastdir
        f = open('lastdir.txt','w')
        f.writelines( [data_directory] )
        f.close()

        self.spefile = MyPrincetonSPEFile(fname).return_Array()
        self.spefile = self.spefile.astype(np.float)

        # is there a blank file around?
        blankname = os.path.split(fname)[0]+os.path.sep+'blank-'+os.path.split(fname)[1]
        if os.path.isfile(blankname):
            self.blankfile = MyPrincetonSPEFile(blankname).return_Array()
            self.blankfile = self.blankfile.astype(np.float)
            self.boolimage = np.ones((self.blankfile.shape[1],self.blankfile.shape[2]), dtype=np.bool)*True

        self.filtered = np.ones_like( self.spefile ) * np.nan
        self.clusters = np.ones_like( self.spefile ) * np.nan
        self.cluster_positions = [None] * self.spefile.shape[0]

        self.frameaverage = np.mean( self.spefile, axis=0 )

        self.frameSlider.setMinimum(0)
        self.frameSlider.setMaximum(self.spefile.shape[0]-1)
        self.frameSlider.setSingleStep(1)
        self.frameSlider.setPageStep(10)
        self.frameSlider.setValue(0)

        self.imageview.myaxes.imshow( self.spefile[0], interpolation='nearest' )
        self.imageview.figure.canvas.draw()

#        self.frameAverageToggleButton.setEnabled(True)
        self.statusbar.showMessage("SPE file loaded. %d frames" % self.spefile.shape[0])

    def reset(self):
        self.spefile = MyPrincetonSPEFile(self.SPEFileName).return_Array()
        self.spefile = self.spefile.astype(np.float)

        self.filtered = np.ones_like( self.spefile ) * np.nan
        self.clusters = np.ones_like( self.spefile ) * np.nan
        self.cluster_positions = [None] * self.spefile.shape[0]

        self.frameaverage = np.mean( self.spefile, axis=0 )

        self.frameSlider.setMinimum(0)
        self.frameSlider.setMaximum(self.spefile.shape[0]-1)
        self.frameSlider.setSingleStep(1)
        self.frameSlider.setPageStep(10)
        self.frameSlider.setValue(0)

        self.imageview.myaxes.imshow( self.spefile[0], interpolation='nearest' )
        self.imageview.figure.canvas.draw()

#        self.frameAverageToggleButton.setEnabled(True)
        self.statusbar.showMessage("Back to square one.")


    def drawFrame(self):
        frame = self.frameSlider.value()
        self.imageview.clear()
        if self.showWhatComboBox.currentIndex()==0:
            self.imageview.myaxes.imshow( self.spefile[frame], interpolation='nearest', cmap=cm.jet )

        elif self.showWhatComboBox.currentIndex()==1:
            self.imageview.myaxes.imshow( self.filtered[frame], interpolation='nearest'  )

        elif self.showWhatComboBox.currentIndex()==2:
            self.imageview.myaxes.imshow( self.clusters[frame], interpolation='nearest' )
        else:
            raise hell

        if self.showBoolimageCheckBox.isChecked():
            if hasattr( self, 'boolimage' ):
                print 'yup'
                self.imageview.myaxes.imshow( self.boolimage.astype(np.float)*np.max(self.spefile[frame]), \
                                                  interpolation='nearest', alpha=.6, cmap=cm.gray, zorder=9 )

        if self.highlightClusterPositionsCheckBox.isChecked():
            if not self.cluster_positions[frame]==None:                
                self.imageview.myaxes.plot( self.cluster_positions[frame][:,1], \
                                     self.cluster_positions[frame][:,0], 'ro', alpha=.5, markersize=5 )
                self.imageview.myaxes.set_xlim( [0, self.spefile.shape[2]] )
                self.imageview.myaxes.set_ylim( [self.spefile.shape[1], 0] )

            if hasattr(self,'clusterindices'):
                print self.uniqueclusters.size
                print self.allclusters.shape[0]
                print np.unique( np.hstack(self.clusterindices) ).size
                # self.imageview.myaxes.plot( self.uniqueclusters[seenatleasttwice,1], \
                #                                 self.uniqueclusters[seenatleasttwice,0], \
                #                                 'k^', alpha=.5, markersize=6, \
                #                                 markerfacecolor='none' )
                for ci in self.clusterindices[frame]:
                    self.imageview.myaxes.text( self.allclusters[ci,1], \
                                                    self.allclusters[ci,0], \
                                                    '%d' % ci )

        self.imageview.figure.canvas.draw()
        self.statusbar.showMessage("frame %d/%d" % (frame+1,self.spefile.shape[0]))


    def frameAverageToggle(self):
        if self.frameAverageToggleButton.isChecked():            
            self.imageview.clear()
            self.imageview.myaxes.imshow( self.frameaverage, interpolation='nearest' )
            self.imageview.figure.canvas.draw()
            self.frameSlider.setEnabled(False)
            self.statusbar.showMessage("Showing frame average.")
        else:
            self.frameSlider.setEnabled(True)
            self.drawFrame()

    def applyWhiteTophatFilter(self):
        if self.frameAverageToggleButton.isChecked():
            image = self.frameaverage
        else:
            frame = self.frameSlider.value()
            if self.showWhatComboBox.currentIndex()==0:
                image = self.spefile[frame]
            else:
                image = self.filtered[frame]

        cutfreq2 = self.highFreqLineEdit.text().toDouble()[0]

        print image.shape
        newimage = smdetect.white_tophat_filter( image, int(cutfreq2) )
        self.filtered[frame] = newimage
        self.showWhatComboBox.setCurrentIndex(1)
        self.drawFrame()


    # def applyGaussianFilter(self):
    #     if self.frameAverageToggleButton.isChecked():
    #         image = self.frameaverage
    #     else:
    #         frame = self.frameSlider.value()
    #         image = self.spefile[frame]

    #     cutfreq2 = self.highFreqLineEdit.text().toDouble()[0]

    #     newimage = smdetect.gaussian_filter( image, cutfreq2 )
    #     self.filtered[frame] = newimage
    #     self.showWhatComboBox.setCurrentIndex(1)
    #     self.drawFrame()

    def applyFourierFilter(self):
        doWhat = self.filterOperationComboBox.currentText()
        #     image = self.frameaverage
        # else:
        frame = self.frameSlider.value()
        image = self.spefile[frame]

        cutfreq1 = self.lowFreqLineEdit.text().toDouble()[0]
        cutfreq2 = self.highFreqLineEdit.text().toDouble()[0]
        # if self.invertedCheckBox.isChecked():
        #     operation = 'inverted'
        # else:
        #     operation = 'normal'

        newimage, self.mask = smdetect.fourier_filter(image, cutfreq1, cutfreq2, operation=doWhat)
        frame = self.frameSlider.value()
        self.filtered[frame] = newimage
        self.showWhatComboBox.setCurrentIndex(1)
        self.drawFrame()

        # draw mask
        self.maskview.myaxes.cla()
        self.maskview.myaxes.imshow( self.mask, interpolation='nearest' )
        self.maskview.figure.canvas.draw()


    def findClusters(self):
        threshold = float( self.clusterIntensityThresholdSlider.value() )

        frame = self.frameSlider.value()

        clusterimage = self.filtered[frame] > np.std(self.filtered[frame].flatten())*threshold
        clusterimage = -clusterimage.astype(np.int)

        cimage = clusterimage.copy()
        # clean up boundaries (got Fourier spills there)
        bskip=10
        cimage[:bskip,:]  = 0   # top edge
        cimage[:,:bskip]  = 0   # left edge
        cimage[-bskip:,:] = 0   # bottom edge
        cimage[:,-bskip:] = 0   # right edge
    
        cimage, nclusters = smdetect.mark_all_clusters(cimage)

        self.clusters[frame] = cimage
        self.statusbar.showMessage("Found %d clusters." % nclusters)

        sizeThreshold = float( self.clusterSizeThresholdSlider.value() )
        pos = smdetect.find_cluster_positions( cimage, size_threshold=sizeThreshold )
        self.cluster_positions[frame] = pos

        self.showWhatComboBox.setCurrentIndex(0)
        self.highlightClusterPositionsCheckBox.setChecked(True)
        self.drawFrame()


    def validateClusters(self):

        distance_threshold = self.clusterDistanceThresholdSlider.value()

        ### WITHIN-FRAME DISTANCE THRESHOLD ###
        # first, we need to clean up clusters _within_ each frame
        for frame,cp in enumerate(self.cluster_positions):
            print '=== %d =====================================' % frame
            print cp.shape

            valid = np.arange( cp.shape[0] )
            for i,c in enumerate(cp):
                for j in range(i+1,cp.shape[0]):
                    d = np.sqrt(np.sum((c-cp[j,:])**2))
                    if d < 2*distance_threshold:
                        valid[j] == i
            print "frame %d, valid: %d/%d" % (frame, np.unique(valid).size, cp.shape[0])
            # retain only valid ones
            self.cluster_positions[frame] = cp[ np.unique(valid) ]

        # next we check for multiple occurences of the same cluster across frames
        
        ### ACROSS-FRAMES DISTANCE THRESHOLD ###
        # first make a list of all clusters
        allclusters = np.vstack( self.cluster_positions )
        # and an indexing array 
        clusterindices = np.arange( allclusters.shape[0] )

        for i in range(clusterindices.size):
            for j in range(i+1, clusterindices.size):
                d = np.sqrt( np.sum((allclusters[ clusterindices[i] ] - allclusters[ clusterindices[j] ])**2) )
                if d < distance_threshold:
                    clusterindices[j] = clusterindices[i]

        Nuniqueclusters = np.unique( clusterindices ).size
        self.uniqueclusters = np.unique( clusterindices )
        self.allclusters = allclusters[ self.uniqueclusters ]
        ucsorted = np.sort( self.uniqueclusters )        
        for i,ci in enumerate(ucsorted):
            np.place( self.uniqueclusters, self.uniqueclusters==ci, i )
            np.place( clusterindices, clusterindices==ci, i )

        # self.uniqueclusters = np.unique( clusterindices )
        # self.allclusters = allclusters

        # split clusterindices into frames
        frameends = np.cumsum( [ cp.shape[0] for cp in self.cluster_positions[:-1] ] )
        cis = np.split(clusterindices, frameends)

        # clean up each frame
        self.clusterindices = []
        for i,f in enumerate(cis):
            self.clusterindices.append( np.unique(f) )

        ### FRAME THRESHOLD ###
        # now go through all unique clusters and count the frames they appear in
        uccount = np.zeros( (self.uniqueclusters.size,) )
        for i,uc in enumerate(self.uniqueclusters):
            for ci in self.clusterindices:
                if np.any( uc==ci ): uccount[i] += 1 

        # clean up again
        frameThreshold = self.clusterFrameThresholdSlider.value()

        #self.allclusters = self.allclusters[ uccount >= frameThreshold ]
        #self.uniqueclusters = self.uniqueclusters[ uccount >= frameThreshold ]

        for i,ci in enumerate(self.clusterindices):
            newci = []
            for j in ci:
                if uccount[j] >= frameThreshold:
                    newci.append( j )
            self.clusterindices[i] = np.array( newci )

        self.validclusters = self.allclusters[ uccount >= frameThreshold ]
        print 'number of validclusters: %d' % self.validclusters.shape[0]

        self.drawFrame()
 

    def writeCoordinates(self):
        fname = 'coords.txt'
        f = open(fname,'w')
        f.writelines( ['%f\t%f\n' % (vc[0],vc[1]) for vc in self.validclusters] )
        f.close()


    def correlateFrames(self):
#        averageFrame = np.sum( np.log(self.spefile), axis=0 )
        
        if self.showWhatComboBox.currentIndex() == 0:
            what = self.spefile
        elif self.showWhatComboBox.currentIndex() == 1:
            what = self.filtered
        else:
            raise hell

        Nrows = self.correlateNRowsSpinBox.value()
        Ncols = self.correlateNColumnsSpinBox.value()

        a = np.zeros( (self.spefile.shape[0], self.spefile.shape[1]-Nrows+1, self.spefile.shape[2]) )
        for i in range(Nrows):
            # print i
            # print self.spefile.shape[1]-(Nrows-i)+1
            a += np.log(self.spefile[ :, i:self.spefile.shape[1]-(Nrows-i)+1, : ])

        b = np.zeros( (self.spefile.shape[0], self.spefile.shape[1], self.spefile.shape[2]-Ncols+1) )
        for i in range(Ncols):
            # print i
            # print self.spefile.shape[2]-(Ncols-i)+1
            b += np.log(self.spefile[ :, :, i:self.spefile.shape[2]-(Ncols-i)+1 ])

        c = a[:,:,:self.spefile.shape[2]-Ncols+1] + b[:,:self.spefile.shape[1]-Nrows+1,:]

        # a = np.log(self.spefile[:,0:-1,:])+np.log(self.spefile[:,1:,:])
        # b = np.log(self.spefile[:,:,0:-1])+np.log(self.spefile[:,:,1:])
        # c = a[:,:,:-1]+b[:,:-1,:]
        self.spefile = np.exp(c)

        self.filtered = np.ones_like( self.spefile ) * np.nan
        self.clusters = np.ones_like( self.spefile ) * np.nan

        # frame = self.frameSlider.value()
        # self.spefile[frame] = (self.spefile[frame] > 6000)


        self.drawFrame()


    def filter_and_cluster_all_frames(self):
        starttime = time.time()

        if self.useCorrelatorCheckBox.isChecked(): 
            self.correlateFrames()

        for f in range(self.spefile.shape[0]):
            self.frameSlider.setValue(f)
            self.applyFourierFilter()
            self.findClusters()

        self.validateClusters()
        self.writeCoordinates()
        
        self.statusbar.showMessage("Done. Took %f seconds." % (time.time()-starttime) )

    def createLineAverages(self):
        lineAverageImages = [np.mean(line,axis=0) for line in np.array_split( self.spefile, 8 )]
        self.spefile = np.array(lineAverageImages)
        self.filtered = np.ones_like( self.spefile ) * np.nan
        self.clusters = np.ones_like( self.spefile ) * np.nan
        self.cluster_positions = [None] * self.spefile.shape[0]

        print self.spefile.shape
        self.frameSlider.setMaximum( self.spefile.shape[0]-1 )
        self.frameSlider.setValue( 0 )
        self.drawFrame()


    def clusterIntensityThresholdChanged(self):
        t = self.clusterIntensityThresholdSlider.value()
        self.clusterIntensityThresholdLabel.setText( "%d std" % t )

        threshold = float( self.clusterIntensityThresholdSlider.value() )
        frame = self.frameSlider.value()
        clusterimage = self.filtered[frame] > np.std(self.filtered[frame].flatten())*threshold
        self.clusters[frame] = clusterimage
        self.showWhatComboBox.setCurrentIndex(2)
        self.drawFrame()

#        self.findClusters()

    def clusterSizeThresholdChanged(self):
        t = self.clusterSizeThresholdSlider.value()
        self.clusterSizeThresholdLabel.setText( "%d std" % t )
#        self.findClusters()

    def crop(self):
        if self.imageview.x0 > self.imageview.x1:
            x0 = self.imageview.x1
            x1 = self.imageview.x0
        else:
            x0 = self.imageview.x0
            x1 = self.imageview.x1
        if self.imageview.y0 > self.imageview.y1:
            y0 = self.imageview.y1
            y1 = self.imageview.y0
        else:
            y0 = self.imageview.y0
            y1 = self.imageview.y1

        self.spefile  = self.spefile[:,y0:y1+1,x0:x1+1]
        self.filtered = self.filtered[:,y0:y1+1,x0:x1+1]
        self.clusters = self.clusters[:,y0:y1+1,x0:x1+1]
        self.drawFrame()


if __name__=='__main__':
    app = QtGui.QApplication(sys.argv)
    instance = sm_detector_gui_logic(app=app)
    instance.main()
    sys.exit(app.exec_())
