# coding=utf-8

from PyQt4 import QtCore,QtGui
import sys, os, time
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
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

        self.header = {}

    def main(self):
        self.show()

    def connectActions(self):
        """Connect the user interface controls to the logic """
        self.loadFilePushButton.clicked.connect( self.loadFile )
        self.applyFilterPushButton.clicked.connect( self.applyFourierFilter )
#        self.findClustersPushButton.clicked.connect( self.findClusters )
        self.findClustersNewPushButton.clicked.connect( self.findClusters2 )
        self.frameSlider.valueChanged.connect( self.drawFrame )
        self.frameAverageToggleButton.clicked.connect( self.frameAverageToggle )
        self.showWhatComboBox.activated.connect( self.drawFrame )
        self.logScaleComboBox.activated.connect( self.drawFrame )
        self.cmapComboBox.activated.connect( self.drawFrame )
        self.highlightClusterPositionsCheckBox.stateChanged.connect( self.drawFrame )
        self.cropPushButton.clicked.connect( self.crop )
        self.resetPushButton.clicked.connect( self.reset )
        self.loadDefaultPushButton.clicked.connect( self.loadDefault )
        self.clusterIntensityThresholdSlider.valueChanged.connect( self.clusterIntensityThresholdChanged )
#        self.clusterBGIntensityThresholdSlider.valueChanged.connect( self.clusterBGIntensityThresholdChanged )
        self.minClusterSizeThresholdSlider.valueChanged.connect( self.minClusterSizeThresholdChanged )
        self.maxClusterSizeThresholdSlider.valueChanged.connect( self.maxClusterSizeThresholdChanged )
        self.createLineAveragePushButton.clicked.connect( self.createLineAverages )
        self.toolBox.currentChanged.connect( self.toolBoxPageChanged )
        self.applyCorrelatorPushButton.clicked.connect( self.correlateFrames )
        self.addFitExclusionZonePushButton.clicked.connect( self.addFitExclusionZone )
        self.resetFitExclusionZonePushButton.clicked.connect( self.resetFitExclusionZone )
        self.blankFitPushButton.clicked.connect( self.blankFit )
        self.showBoolimageCheckBox.stateChanged.connect( self.drawFrame )
        self.histogramToggleButton.clicked.connect( self.drawFrame )
        self.useOriginalDataToggleButton.clicked.connect( self.useOriginalDataForClusters )
        self.useFourierdDataToggleButton.clicked.connect( self.useFourierdDataForClusters )
        self.specialForOldDataPushButton.clicked.connect( self.special_for_old_data )
        self.specialForOldDataExPushButton.clicked.connect( self.special_for_old_data_ex )
        self.specialForOldDataEmPushButton.clicked.connect( self.special_for_old_data_em )
        self.specialForOldDataMapPushButton.clicked.connect( self.special_for_old_data_map )
        self.specialForOldDataAllPushButton.clicked.connect( self.special_for_old_data_all )
        self.clusterDistanceThresholdSlider.valueChanged.connect( self.clusterDistanceThresholdSliderChanged)
#        self.clusterCrossFrameDistanceThresholdSlider.valueChanged.connect( self.clusterValidationSlidersChanged)
        self.clusterFrameThresholdSlider.valueChanged.connect( self.clusterFrameThresholdSliderChanged)


#        self.clusterMergeAcrossFramesCheckBox.stateChanged.connect( self.clusterValidationStateChanged )

        self.writeCoordsPushButton.clicked.connect( self.writeCoordinates )


    def clusterDistanceThresholdSliderChanged(self):
        c = self.clusterDistanceThresholdSlider.value()
        self.clusterDistanceThresholdLabel.setText( str(c) + ' pixel' )

    def clusterFrameThresholdSliderChanged(self):
        c = self.clusterFrameThresholdSlider.value()
        if c==1:
            self.clusterFrameThresholdLabel.setText( str(c) + ' frame' )
        else:
            self.clusterFrameThresholdLabel.setText( str(c) + ' frames' )

    def useOriginalDataForClusters( self ):
        if self.useOriginalDataToggleButton.isChecked():
            self.useFourierdDataToggleButton.setChecked(False)
        else:
            self.useFourierdDataToggleButton.setChecked(True)

    def useFourierdDataForClusters( self ):
        if self.useFourierdDataToggleButton.isChecked():
            self.useOriginalDataToggleButton.setChecked(False)
        else:
            self.useOriginalDataToggleButton.setChecked(True)

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
            
            # print res[0][0]            
            # print res[0][1]            
            # print np.histogram( collapsed_blank_data )

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
        print "current toolbox index: %d" % i
        for j in range(self.toolBox.count()):
            if not i==j:
                text = self.toolBox.itemText(j)
                self.toolBox.setItemText(j, '>>'+text[2:])
            else:
                text = self.toolBox.itemText(j)
                self.toolBox.setItemText(j, 'vv'+text[2:])

        if i==0:          # correlator tab
            self.showWhatComboBox.setCurrentIndex(0)
            self.drawFrame()
        elif i==1:        # fourier tab
            if not np.isnan(self.filtered[0,0,0]):
                self.showWhatComboBox.setCurrentIndex(1)
            else:
                self.showWhatComboBox.setCurrentIndex(0)
            self.drawFrame()
        elif i==2:        # cluster tab
            self.clusterIntensityThresholdChanged()
        else:
            raise thefuckisthis


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
            self.filename = os.path.split( os.path.normpath(fname) )[1]
            self.filedir  = os.path.split( os.path.normpath(fname) )[0] + os.path.sep
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
        print 'Looking for blank file with name: %s ...' % blankname, 
        if os.path.isfile(blankname):
            print 'found it.'
            self.blankfile = MyPrincetonSPEFile(blankname).return_Array()
            self.blankfile = self.blankfile.astype(np.float)
            self.boolimage = np.ones((self.blankfile.shape[1],self.blankfile.shape[2]), dtype=np.bool)*True
            self.blankFitGroupBox.setTitle("blank fit")
            self.blankFitGroupBox.setEnabled(True)
        else:
            print 'nothing with that exact name'
            self.blankFitGroupBox.setTitle("NO BLANK (looked for %s)" % ('blank-'+os.path.split(fname)[1]))
            self.blankFitGroupBox.setEnabled(False)

        self.filtered = np.ones_like( self.spefile ) * np.nan
        self.clusters = np.ones_like( self.spefile ) * np.nan
        self.cluster_positions = [None] * self.spefile.shape[0]

        self.frameaverage = np.mean( self.spefile, axis=0 )

        self.frameSlider.setMinimum(0)
        self.frameSlider.setMaximum(self.spefile.shape[0]-1)
        self.frameSlider.setSingleStep(1)
        self.frameSlider.setPageStep(10)
        self.frameSlider.setValue(0)

        self.drawFrame()
#        self.imageview.myaxes.imshow( self.spefile[0], interpolation='nearest' )
#        self.imageview.figure.canvas.draw()

#        self.frameAverageToggleButton.setEnabled(True)
        self.statusbar.showMessage("SPE file loaded. %d frames" % self.spefile.shape[0])

    def reset(self):
        self.spefile = MyPrincetonSPEFile(self.filedir+self.filename).return_Array()
        self.spefile = self.spefile.astype(np.float)

        if hasattr( self, 'boolimage' ):
            self.boolimage = np.ones((self.blankfile.shape[1],self.blankfile.shape[2]), dtype=np.bool)*True

        self.filtered = np.ones_like( self.spefile ) * np.nan
        self.clusters = np.ones_like( self.spefile ) * np.nan
        self.cluster_positions = [None] * self.spefile.shape[0]

        if hasattr( self, 'clusterindices' ): delattr( self, 'clusterindices' )
        if hasattr( self, 'validclusters' ): delattr( self, 'validclusters' )

        self.frameaverage = np.mean( self.spefile, axis=0 )

        self.frameSlider.setMinimum(0)
        self.frameSlider.setMaximum(self.spefile.shape[0]-1)
        self.frameSlider.setSingleStep(1)
        self.frameSlider.setPageStep(10)
        self.frameSlider.setValue(0)

        # self.imageview.myaxes.imshow( self.spefile[0], interpolation='nearest' )
        # self.imageview.figure.canvas.draw()        
        self.drawFrame()

#        self.frameAverageToggleButton.setEnabled(True)
        self.statusbar.showMessage("Back to square one.")
        

    def drawFrame_old(self):
        frame = self.frameSlider.value()
        cmap = getattr( cm, str(self.cmapComboBox.currentText()) )
        bins = 200
        self.imageview.clear()

        if str(self.logScaleComboBox.currentText())=='normal scale':
            modfun = lambda x: x
        elif str(self.logScaleComboBox.currentText())=='log scale':
            modfun = lambda x: np.log( x+2*np.abs(np.min(x)) )


        if self.showWhatComboBox.currentIndex()==0:
            if self.histogramToggleButton.isChecked():
                self.imageview.myaxes.hist( modfun(self.spefile[frame]).flatten(), bins=bins )
            else:
                self.imageview.myaxes.imshow( modfun( self.spefile[frame] ), interpolation='nearest', cmap=cmap, \
                                                  vmin=modfun(np.min(self.spefile)), vmax=modfun(np.max(self.spefile)) )

        elif self.showWhatComboBox.currentIndex()==1:
            if self.histogramToggleButton.isChecked():
                self.imageview.myaxes.hist( modfun(self.filtered[frame]).flatten(), 500 )
            else:
                self.imageview.myaxes.imshow( modfun( self.filtered[frame] ), interpolation='nearest', cmap=cmap  )

        elif self.showWhatComboBox.currentIndex()==2:
            if self.histogramToggleButton.isChecked():
                self.imageview.myaxes.hist( self.clusters[frame].flatten(), 500 )
            else:
                if np.max(self.clusters[frame])==1:   # most likely got the boolean image
                    self.imageview.myaxes.imshow( self.clusters[frame], interpolation='nearest', vmax=2, cmap=cmap )
                else:
                    self.imageview.myaxes.imshow( self.clusters[frame], interpolation='nearest', cmap=cmap )
        else:
            raise hell

        if self.showBoolimageCheckBox.isChecked():
            if hasattr( self, 'boolimage' ):
                print 'yup'
                self.imageview.myaxes.imshow( self.boolimage.astype(np.float)*np.max(self.spefile[frame]), \
                                                  interpolation='nearest', alpha=.6, cmap=cm.gray, zorder=9 )

        if self.highlightClusterPositionsCheckBox.isChecked():
            if hasattr( self, 'validclusters' ):

                self.imageview.myaxes.plot( self.cluster_positions[frame][:,1], \
                                                self.cluster_positions[frame][:,0], 'ko', \
                                                markersize=8, markerfacecolor='none', zorder=5 )

        # fix zoom (cluster marker plot can screw it up)
        xlim = self.imageview.myaxes.get_xlim()
        ylim = self.imageview.myaxes.get_ylim()
        if (np.min(xlim)<0) or (np.max(xlim)>self.spefile.shape[2]):
            self.imageview.myaxes.set_xlim( [0, self.spefile.shape[2]] )
        if (np.min(ylim)<0) or (np.max(ylim)>self.spefile.shape[1]):
            self.imageview.myaxes.set_ylim( [self.spefile.shape[1], 0] )

        self.imageview.figure.canvas.draw()
        self.statusbar.showMessage("frame %d/%d" % (frame+1,self.spefile.shape[0]))


    def drawFrame(self):
        frame = self.frameSlider.value()
        cmap = getattr( cm, str(self.cmapComboBox.currentText()) )
        bins = 200

        if str(self.logScaleComboBox.currentText())=='normal scale':
            modfun = lambda x: x
        elif str(self.logScaleComboBox.currentText())=='log scale':
            modfun = lambda x: np.log( x+2*np.abs(np.min(x)) )


        if self.histogramToggleButton.isChecked():
            self.imageview.clear()
            if self.showWhatComboBox.currentIndex()==0:
                self.imageview.myaxes.hist( modfun(self.spefile[frame]).flatten(), bins=bins )
            elif self.showWhatComboBox.currentIndex()==1:            
                self.imageview.myaxes.hist( modfun(self.filtered[frame]).flatten(), 500 )
            elif self.showWhatComboBox.currentIndex()==2:
                self.imageview.myaxes.hist( self.clusters[frame].flatten(), 500 )
        else:
            if self.showWhatComboBox.currentIndex()==0:
                if not hasattr(self.imageview, 'current_image'):
                    im = self.imageview.myaxes.imshow( modfun(self.spefile[frame]), \
                                                           interpolation='nearest', cmap=cmap, \
                                                           vmin=modfun(np.min(self.spefile)), \
                                                           vmax=modfun(np.max(self.spefile)) )
                    self.imageview.current_image = im
                else:
                    self.imageview.current_image.set_data( modfun(self.spefile[frame]) )
                    self.imageview.current_image.set_clim( modfun(np.min(self.spefile)), \
                                                               modfun(np.max(self.spefile)) )
                    self.imageview.current_image.set_cmap( cmap )

            elif self.showWhatComboBox.currentIndex()==1:
                if not hasattr(self.imageview, 'current_image'):
                    im = self.imageview.myaxes.imshow( modfun(self.filtered[frame]), \
                                                           interpolation='nearest', cmap=cmap  )
                    self.imageview.current_image = im
                else:
                    self.imageview.current_image.set_data( modfun(self.filtered[frame]) )
                    self.imageview.current_image.set_clim( modfun(np.min(self.filtered)), \
                                                               modfun(np.max(self.filtered)) )
                    self.imageview.current_image.set_cmap( cmap )

            elif self.showWhatComboBox.currentIndex()==2:
                if not hasattr(self.imageview, 'current_image'):
                    im = self.imageview.myaxes.imshow( self.clusters[frame], \
                                                           interpolation='nearest', cmap=cmap )
                    self.imageview.current_image = im
                else:
                    self.imageview.current_image.set_data( self.clusters[frame] )
                    self.imageview.current_image.set_clim( 0, np.max(self.clusters)-np.min(self.clusters) )
                    print ( 0, np.max(self.clusters)-np.min(self.clusters) )
                    self.imageview.current_image.set_cmap( cmap )
            else:
                raise hell

        if self.showBoolimageCheckBox.isChecked():
            if hasattr( self, 'boolimage' ):
                if not hasattr( self.imageview, 'bool_image' ):
                    im = self.imageview.myaxes.imshow( \
                        self.boolimage.astype(np.float)*np.max(self.spefile[frame]), \
                            interpolation='nearest', alpha=.6, cmap=cm.gray, zorder=9 )
                    self.imageview.bool_image = im
                else:
                    self.imageview.bool_image.set_data( \
                        self.boolimage.astype(np.float)*np.max(self.spefile[frame]) )
        else:
            if hasattr( self.imageview, 'bool_image' ):
                self.imageview.bool_image.set_data( None )


        if self.highlightClusterPositionsCheckBox.isChecked():
            if hasattr( self, 'validclusters' ):
                if not hasattr( self.imageview, 'cluster_markers' ):
                    cmarkers = self.imageview.myaxes.plot( self.validclusters[:,1], \
                                                         self.validclusters[:,0], 'ko', \
                                                         markersize=8, markerfacecolor='none', zorder=5 )
                    self.imageview.cluster_markers = cmarkers[0]
                else:
                    self.imageview.cluster_markers.set_data( self.validclusters[:,1], \
                                                                 self.validclusters[:,0] )

        # fix zoom (cluster marker plot can screw it up)
        xlim = self.imageview.myaxes.get_xlim()
        ylim = self.imageview.myaxes.get_ylim()
        if (np.min(xlim)<0) or (np.max(xlim)>self.spefile.shape[2]):
            self.imageview.myaxes.set_xlim( [0, self.spefile.shape[2]] )
        if (np.min(ylim)<0) or (np.max(ylim)>self.spefile.shape[1]):
            self.imageview.myaxes.set_ylim( [self.spefile.shape[1], 0] )

        # and redraw
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
        for iframe in range(self.spefile.shape[0]):
            self.frameSlider.setValue(iframe)
        
            doWhat = self.filterOperationComboBox.currentText()
            #     image = self.frameaverage
            # else:
            image = self.spefile[iframe]

            cutfreq1 = self.lowFreqLineEdit.text().toDouble()[0]
            cutfreq2 = self.highFreqLineEdit.text().toDouble()[0]
            # if self.invertedCheckBox.isChecked():
            #     operation = 'inverted'
            # else:
            #     operation = 'normal'

            newimage, self.mask = smdetect.fourier_filter(image, cutfreq1, cutfreq2, operation=doWhat)
            self.filtered[iframe] = newimage
            self.showWhatComboBox.setCurrentIndex(1)
            self.drawFrame()

            # draw mask
            self.maskview.myaxes.cla()
            self.maskview.myaxes.imshow( self.mask, interpolation='nearest' )
            self.maskview.figure.canvas.draw()

    def findClusters2(self):
        threshold = float( self.clusterIntensityThresholdSlider.value()/10.0 )
        self.clusterIntensityThresholdChanged()
        # we'll work with the average image 
        if self.useOriginalDataToggleButton.isChecked():
            data = self.spefile
        else:
            data = self.filtered
        image = np.mean( data, axis=0 )
        # true/false array showing where pixel are above threshold
        clusterimage = image > np.std(image.flatten())*threshold
        # turn each True element into -1, False into 0 (in preparation for the next line)
        clusterimage = -clusterimage.astype(np.int)
        # mark clusters
        clusterimage, nclusters = smdetect.mark_all_clusters(clusterimage)
#        plt.imshow(clusterimage)
#        plt.savefig(self.filedir+'clusterimage.png')
        print self.filedir+'clusterimage.png'

        Ngoodframesrequired = self.clusterFrameThresholdSlider.value()

        #### now check if clusters are in each image
        fpix = .8     # fraction of pixel which must be above level to declare a frame 'good' 
        # for each cluster
        for ni in range(1,nclusters+1):
            # find out where in the image it is,
            clusterbools = clusterimage==ni
            # and how many pixel it has
            clusterNpix = np.sum(clusterbools)
            # go through all frames and find where it is about level
            goodframes = 0
            for frame in range( data.shape[0] ):
                # we trust that self.clusters was set by clusterIntensityThresholdChanged()
#                if np.sum( self.clusters[frame][clusterbools] ) >= fpix*clusterNpix:
                if np.sum( self.clusters[frame][clusterbools] ) >= self.minClusterSizeThresholdSlider.value():
                    goodframes += 1
            # did we have enough good frames? No? Then remove this cluster from the image:
            if goodframes < Ngoodframesrequired:
                clusterimage[clusterbools] = 0

        self.statusbar.showMessage("Found %d clusters." % nclusters)
        # determine cluster centers and apply size criteria
        minSizeThreshold = float( self.minClusterSizeThresholdSlider.value() )
        maxSizeThreshold = float( self.maxClusterSizeThresholdSlider.value() )
        pos = smdetect.find_cluster_positions( clusterimage, \
                                                   min_size_threshold=minSizeThreshold, \
                                                   max_size_threshold=maxSizeThreshold, \
                                                   originalimage=image )
        for frame in range(self.spefile.shape[0]):
            self.cluster_positions[frame] = pos

        self.validclusters = np.vstack( self.cluster_positions )

        ### clean up potential redundancies ###
        # convert array to list of tuples, so that 'set' operation can work
        vc = [(a[0],a[1]) for a in self.validclusters]
        # create set (thereby removing redundancies), then back to array
        self.validclusters = np.array([ p for p in set(vc) ])

        
        ### min distance check
        distance_threshold = self.clusterDistanceThresholdSlider.value()
        # distance matrix, one if ok, zero if too close
        validity = np.ones( (self.validclusters.shape[0],self.validclusters.shape[0]), dtype=np.int )
        # will check only the upper triangular matrix
        for i,c in enumerate(self.validclusters):
            for j in range(i+1,self.validclusters.shape[0]):
                d = np.sqrt(np.sum((c-self.validclusters[j,:])**2))
                if d < distance_threshold:
                    validity[i,j] = 0
        # the lower triangular matrix and main diagonal can remain all ones: this way we
        # will check for and remove only one of two clusters which are too close to each other
        # Ok, now find the min of each column, and convert to True/False
        validity = np.min( validity, axis=0 ).astype(np.bool)
        # and only keep those clusters marked as valid
        self.validclusters = self.validclusters[validity]

        # now plotting: 
        if self.useOriginalDataToggleButton.isChecked():
            self.showWhatComboBox.setCurrentIndex(0)
        else:
            self.showWhatComboBox.setCurrentIndex(1)
        self.highlightClusterPositionsCheckBox.setChecked(True)
        self.drawFrame()

#        self.validateClusters()

        # if self.clustersMustAppearInAtLeastNFramesCheckBox.isChecked():
        #     bgthreshold = float( self.clusterBGIntensityThresholdSlider.value() )/10.0
        #     self.bgexclusionmap = np.zeros( data[0,:,:].shape, dtype=np.bool )*False
        #     for frame in range(data.shape[0]):
        #         self.bgexclusionmap = np.logical_or( self.bgexclusionmap, data[frame] > np.std(data[frame].flatten())*bgthreshold )

        #     self.imageview.clear()
        #     self.imageview.myaxes.imshow( self.bgexclusionmap )
        #     self.imageview.myaxes.figure.canvas.draw()

#        self.drawFrame()


    def findClusters(self):
        threshold = float( self.clusterIntensityThresholdSlider.value()/10.0 )

        frame = self.frameSlider.value()
        for frame in range( self.spefile.shape[0] ):

            if self.useOriginalDataToggleButton.isChecked():
                image = self.spefile[frame]
            else:
                image = self.filtered[frame]

            clusterimage = image > np.std(image.flatten())*threshold
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

            minSizeThreshold = float( self.minClusterSizeThresholdSlider.value() )
            maxSizeThreshold = float( self.maxClusterSizeThresholdSlider.value() )
            pos = smdetect.find_cluster_positions( cimage, \
                                                       min_size_threshold=minSizeThreshold, \
                                                       max_size_threshold=maxSizeThreshold, \
                                                       originalimage=image )
            self.cluster_positions[frame] = pos

        self.validclusters = np.vstack( self.cluster_positions )

        if self.useOriginalDataToggleButton.isChecked():
            self.showWhatComboBox.setCurrentIndex(0)
        else:
            self.showWhatComboBox.setCurrentIndex(1)

        self.highlightClusterPositionsCheckBox.setChecked(True)

        self.validateClusters()

        self.drawFrame()


    def validateClusters(self):

        ### WITHIN-FRAME DISTANCE THRESHOLD ###
        print "### WITHIN-FRAME DISTANCE THRESHOLD ###"
        if self.clusterWithinFrameMinDistanceCheckBox.isChecked():

            distance_threshold = self.clusterDistanceThresholdSlider.value()

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

            self.validclusters = np.vstack( self.cluster_positions )
                
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print self.validclusters.shape


        # next we check for multiple occurences of the same cluster across frames
        
        ### ACROSS-FRAMES DISTANCE THRESHOLD ###
        print "### ACROSS-FRAMES DISTANCE THRESHOLD ###"
        if self.clusterMergeAcrossFramesCheckBox.isChecked():

            # first make a list of all clusters
            allclusters = np.vstack( self.cluster_positions )
            # and an indexing array 
            clusterindices = np.arange( allclusters.shape[0] )

            for i in range(clusterindices.size):
                for j in range(i+1, clusterindices.size):
                    d = np.sqrt( np.sum((allclusters[ clusterindices[i] ]-allclusters[ clusterindices[j] ])**2) )
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
            self.cluster_positions = []
            for i,f in enumerate(cis):
                self.clusterindices.append( np.unique(f) )
                self.cluster_positions.append( self.allclusters[ np.unique(f) ] )

            self.validclusters = np.vstack( self.cluster_positions )

            print "@@@@@@@@@@@@@@@@@@@@@@@@@@@"
            print self.validclusters.shape
            

        ### FRAME THRESHOLD ###
        print "### FRAME THRESHOLD ###"
        if self.clustersMustAppearInAtLeastNFramesCheckBox.isChecked():
        
            # now go through all unique clusters and count the frames they appear in
            uccount = np.zeros( (self.uniqueclusters.size,) )
            for i,uc in enumerate(self.uniqueclusters):
                for ci in self.clusterindices:
                    if np.any( uc==ci ): uccount[i] += 1 

            # clean up again
            frameThreshold = self.clusterFrameThresholdSlider.value()

            # self.allclusters = self.allclusters[ uccount >= frameThreshold ]
            # self.uniqueclusters = self.uniqueclusters[ uccount >= frameThreshold ]

            for i,ci in enumerate(self.clusterindices):
                newci = []
                for j in ci:
                    if uccount[j] >= frameThreshold:
                        newci.append( j )
                self.clusterindices[i] = np.array( newci )

            self.validclusters = self.allclusters[ uccount >= frameThreshold ]
            print 'number of validclusters: %d' % self.validclusters.shape[0]

            print '###############################'
            print self.validclusters
            print self.validclusters.shape
        self.drawFrame()
 

    def writeCoordinates(self):
#        fname = QtGui.QFileDialog.getSaveFileName( self, 'Specific save file name', directory='' )
        fname = 'spotcoordinates_'+self.filename[:-4]+'.txt'
        fname = self.filedir + fname

        f = open( str(fname),'w')
        # header = {'SPEFileName': self.SPEFileName,
        #           'WaterLevelInSigmas': float(self.bgwaterlevelmap}
        f.writelines( ['%.1f\t%.1f\n' % (vc[0],vc[1]) for vc in self.validclusters] )
        f.close()

        # fname = 'bgexclusionmap_'+self.filename[:-4]+'.txt'
        # fname = self.filedir + fname
        # np.savetxt(fname, self.bgexclusionmap)


    def correlateFrames(self):
#        averageFrame = np.sum( np.log(self.spefile), axis=0 )
        
        # if self.showWhatComboBox.currentIndex() == 0:
        #     what = self.spefile
        # elif self.showWhatComboBox.currentIndex() == 1:
        #     what = self.filtered
        # else:
        #     raise hell

        # Nrows = self.correlateNRowsSpinBox.value()
        # Ncols = self.correlateNColumnsSpinBox.value()

        if np.min(self.spefile) < 0:
            offset = -2*np.min(self.spefile)
        else:
            offset = 0
        self.spefile += offset

#         a = np.ones( (self.spefile.shape[0], self.spefile.shape[1]-Nrows+1, self.spefile.shape[2]) )
#         for i in range(Nrows):
#             # print i
#             # print self.spefile.shape[1]-(Nrows-i)+1
#             a += np.log(self.spefile[ :, i:self.spefile.shape[1]-(Nrows-i)+1, : ])
# #            a *= self.spefile[ :, i:self.spefile.shape[1]-(Nrows-i)+1, : ]

#         b = np.ones( (self.spefile.shape[0], self.spefile.shape[1], self.spefile.shape[2]-Ncols+1) )
#         for i in range(Ncols):
#             # print i
#             # print self.spefile.shape[2]-(Ncols-i)+1
#             b += np.log(self.spefile[ :, :, i:self.spefile.shape[2]-(Ncols-i)+1 ])
# #            b *= self.spefile[ :, :, i:self.spefile.shape[2]-(Ncols-i)+1 ]

#         c = a[:,:,:self.spefile.shape[2]-Ncols+1] + b[:,:self.spefile.shape[1]-Nrows+1,:]

        a = np.log(self.spefile)+np.log( np.roll(self.spefile, 1, axis=1) )
        b = np.log(self.spefile)+np.log( np.roll(self.spefile, 1, axis=2) )
        c = a+b
        if self.fourthRootInCorrelatorCheckBox.isChecked():
            self.spefile = np.power( np.exp(c), .25 ) - offset
        else:
            self.spefile = np.exp(c) - offset
#        self.spefile = c

        print np.sum( np.isnan( self.spefile ) )
        print np.sum( np.isinf( self.spefile ) )
#        print np.histogram(self.spefile.flatten(), 100)# bins=np.linspace(-100,1000))

        self.filtered = np.ones_like( self.spefile ) * np.nan
        self.clusters = np.ones_like( self.spefile ) * np.nan

        # frame = self.frameSlider.value()
        # self.spefile[frame] = (self.spefile[frame] > 6000)

        self.drawFrame()
        print 'done'


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

    def write_motor_file( self, prefix, basename, exangles, emangles ):
        header = """optical element in excitation: L/2 Plate
rotation mode in excitation: stepwise
rotation mode in emission: stepwise
phase offset in deg: 0.000000
full excitation power in mW: 1.0
optical density of filter: 1.000000
excitation power after filter: 1.0
real gain: 1.000000
camera: Cascade
laser wavelength in nm: 488.000000
objective: 60X oil
em correction modulation depth: 0.0
em correction phase: 0.000000
NA: 0.650000
user notes: 
This is a bogus header intended for use with test data.
END-OF-HEADER 
Frame	Excitation Physical Angle	Excitation Polarization Angle	Excitation Power [A.U.]	Emission Angle
"""
        
        fid = open( prefix+'MS-'+basename+'.txt', 'w' )

        fid.write(header)

        frame = 0
        for ema in emangles:
            for exa in exangles:
                fid.write('%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (frame, exa*180/np.pi, exa*180/np.pi, 1.0, ema*180/np.pi))
                frame += 1

        fid.close()


    def special_for_old_data_ex(self):
        self.imageview.c0 = 59
        self.imageview.c1 = 75
        self.imageview.r0 = 240
        self.imageview.r1 = 250
        return self.special_for_old_data()
        
    def special_for_old_data_em(self):
        self.imageview.c0 = 47
        self.imageview.c1 = 68
        self.imageview.r0 = 197
        self.imageview.r1 = 212
        return self.special_for_old_data()

    def special_for_old_data_map(self):
        self.clusterIntensityThresholdSlider.setValue( 40.0 )   # this will also trigger clusterIntensityThresholdChanged()
        boolmap = np.sum( self.clusters, axis=0 ) > 0
        avemap = np.mean( self.spefile, axis=0 )
        avemap[boolmap] = 0                   # set all to zero which are above threshold
#        np.save( self.filedir+self.filename[:-4]+'_avemap.npy', avemap )

        return avemap

    def special_for_old_data_all(self):
        addmap = np.zeros_like( self.spefile[0,:,:] )
        valmap = np.zeros_like( self.spefile[0,:,:] )

        for f in os.listdir(self.filedir):
            if f.endswith(".spe") or f.endswith(".SPE"):
                self.defaultfile = self.filedir+f
                self.loadFile()
                print "------------------------ EX AND EM ANGLES ------------------------"
                exangles_in_rad = self.special_for_old_data_ex()
                emangles_in_rad = self.special_for_old_data_em()
                self.write_motor_file( self.filedir, self.filename[:-4], exangles_in_rad, emangles_in_rad )
#                np.save( self.filedir+self.filename[:-4]+'_exangles_in_rad.npy', exangles_in_rad )
#                np.save( self.filedir+self.filename[:-4]+'_emangles_in_rad.npy', emangles_in_rad )

                self.createLineAverages()
                print "------------------------ MAP -------------------------"
                a = self.special_for_old_data_map()
                addmap[~(a==0)] += 1
                valmap += a
        print set(addmap.flatten())
        print valmap
        blank = valmap/addmap
        np.save( self.filedir+"blank-"+self.filename[:-4]+'.npy', blank )


    def special_for_old_data(self, plots=False):
        if self.imageview.c0 > self.imageview.c1:
            c0 = self.imageview.c1
            c1 = self.imageview.c0
        else:
            c0 = self.imageview.c0
            c1 = self.imageview.c1
        if self.imageview.r0 > self.imageview.r1:
            r0 = self.imageview.r1
            r1 = self.imageview.r0
        else:
            r0 = self.imageview.r0
            r1 = self.imageview.r1

        print [r0,r1,c0,c1]

        abc  = np.sum( np.sum( self.spefile[:,r0:r1+1,c0:c1+1], axis=1 ), axis=1 ) 
        abc -= np.min(abc)
        
        testval = np.max(abc)/10.0     # somewhere between baseline and peaks
        whereabove = np.array([abc > testval]).squeeze()

        framenumbers = np.arange(abc.size)
        if plots:
            import matplotlib.pyplot as plt
            plt.plot( framenumbers, abc, 'b.-' )
            plt.plot( framenumbers, testval*np.ones_like(abc), 'r:' )
            plt.plot( framenumbers[whereabove], abc[whereabove], 'r' )        

        # get the indices of data in each peak separately
        b = []
        c = []
        onpeak = False
        for i in framenumbers[1:-1]:
            if np.logical_and( ~whereabove[i-1], whereabove[i] ): 
                onpeak=True
            if onpeak:
                b.append( i )
            if np.logical_and( whereabove[i], ~whereabove[i+1] ): 
                onpeak=False
                c.append( b )
                b = []

        # peak positions are signal weighted
        peakpositions = []
        for peak in c:
            peakpositions.append( np.sum( abc[peak]*peak ) / np.sum( abc[peak] ) )
    
        # linear fit of peak positions over angles [0, 2*pi, 4*pi, ...]
        angles = np.arange(0,2*np.pi*len(peakpositions),2*np.pi)
        fit = np.linalg.lstsq( np.vstack([angles, np.ones_like(angles)]).T, peakpositions )
        slope,intercept = fit[0]
        print slope, intercept   # units of [frame/angle] and [frame]

        if plots:
            plt.figure()
            plt.plot( angles, peakpositions, 'o' )
            plt.plot( angles, slope*angles+intercept )

        # inverse function 
        frameangles = framenumbers/float(slope) - intercept/float(slope)

        print 1.0/slope, intercept/float(slope)   # units of [angle/frame] and [angle]
        print frameangles
        return frameangles

    def clusterBGIntensityThresholdChanged(self):
        bgthreshold = float( self.clusterBGIntensityThresholdSlider.value() )/10.0
        print bgthreshold
        self.clusterBGIntensityThresholdLabel.setText( "%.1f std" % bgthreshold )

    def clusterIntensityThresholdChanged(self):
        threshold = float( self.clusterIntensityThresholdSlider.value() )/10.0
        print threshold
        self.clusterIntensityThresholdLabel.setText( "%.1f std" % threshold )

        for frame in range( self.spefile.shape[0] ):
#        frame = self.frameSlider.value()
            if self.useOriginalDataToggleButton.isChecked():
                image = self.spefile[frame]
            else:
                image = self.filtered[frame]

            clusterimage = image > np.std(image.flatten())*threshold
            self.clusters[frame] = clusterimage
        self.showWhatComboBox.setCurrentIndex(2)
        self.drawFrame()

#        self.findClusters()

    def minClusterSizeThresholdChanged(self):
        t = self.minClusterSizeThresholdSlider.value()
        self.minClusterSizeThresholdLabel.setText( "%d pixel" % t )

    def maxClusterSizeThresholdChanged(self):
        t = self.maxClusterSizeThresholdSlider.value()
        self.maxClusterSizeThresholdLabel.setText( "%d pixel" % t )
#        self.findClusters()

    def crop(self):
        if self.imageview.c0 > self.imageview.c1:
            c0 = self.imageview.c1
            c1 = self.imageview.c0
        else:
            c0 = self.imageview.c0
            c1 = self.imageview.c1
        if self.imageview.r0 > self.imageview.r1:
            r0 = self.imageview.r1
            r1 = self.imageview.r0
        else:
            r0 = self.imageview.r0
            r1 = self.imageview.r1

        self.spefile   = self.spefile[:,r0:r1+1,c0:c1+1]
        self.blankfile = self.blankfile[:,r0:r1+1,c0:c1+1]
        self.filtered  = self.filtered[:,r0:r1+1,c0:c1+1]
        self.clusters  = self.clusters[:,r0:r1+1,c0:c1+1]
        self.drawFrame()


if __name__=='__main__':
    app = QtGui.QApplication(sys.argv)
    instance = sm_detector_gui_logic(app=app)
    instance.main()
    sys.exit(app.exec_())
