from PyQt4 import QtCore,QtGui
import sys, os, time
import util2dpolim.gui.the2dgui as the2dgui
from sm_detection.sm_detector_gui import sm_detector_gui_logic
import numpy as np
from util2dpolim.movie import Movie
from util2dpolim.misc import *
import matplotlib.cm as cm
from matplotlib.patches import Rectangle

class the2dlogic(QtGui.QMainWindow,the2dgui.Ui_MainWindow):

    def __init__(self,parent=None,app=None):
        """
            Initialization of the class. Call the __init__ for the super classes
        """
        super(the2dlogic,self).__init__(parent)
        self.setupUi(self)
        self.connectActions()
        self.set_up_plots()
        
        self.spefiles = []
        self.m = None
        self.pwd = os.path.dirname(os.path.abspath(__file__))
        self.current_spot = None
        self.app = app
        self.window_sm_detector = None
        self.hi_my_name_is = "the 2d gui logic class!"

    def keyPressEvent(self, event):
        if event.key()==QtCore.Qt.Key_Up:
            self.move_crosshairs('up')
        if event.key()==QtCore.Qt.Key_Down:
            self.move_crosshairs('down')
        if event.key()==QtCore.Qt.Key_Left:
            self.move_crosshairs('left')
        if event.key()==QtCore.Qt.Key_Right:
            self.move_crosshairs('right')
            
    def set_up_plots( self ):
        self.dataview.fig.clear()
        self.dataview.myaxes = self.dataview.fig.add_subplot(111)

    def move_crosshairs( self, direction ):
        r = self.imageview.crosshairs_r
        c = self.imageview.crosshairs_c
        if not self.current_spot==None:
            while self.m.spots[self.current_spot].has_pixel( (r,c) ):
                if direction=='up':   r -= 1
                if direction=='down': r += 1
                if direction=='left':  c -= 1
                if direction=='right': c += 1
        else:
            if direction=='up':   r -= 1
            if direction=='down': r += 1
            if direction=='left':  c -= 1
            if direction=='right': c += 1
        r = np.clip(r,0,self.m.sample_data.datasize[1]-1)
        c = np.clip(c,0,self.m.sample_data.datasize[2]-1)
        self.imageview.crosshairs_r = r
        self.imageview.crosshairs_c = c
        self.imageview.hline.set_ydata( np.array([r,r]) )
        self.imageview.vline.set_xdata( np.array([c,c]) )

        self.crosshair_pick()

    def crosshair_pick(self):

        col = self.imageview.crosshairs_c
        row = self.imageview.crosshairs_r
        print 'row=%f\tcol=%f' % (row,col)
        sys.stdout.flush()        
        
        self.current_spot = None
        if (not self.m==None) and hasattr(self.m,'spots'):
            # find the spot which covers these coordinates
            for si,s in enumerate(self.m.spots):
                if s.has_pixel( (row,col) ): 
                    self.current_spot = si
                    row = s.center[0]
                    col = s.center[1]
                    print 'picked spot #%d' % si
                    sys.stdout.flush()
                    break
            if not self.m.bg_spot_sample==None:
                if self.m.bg_spot_sample.has_pixel( (row,col) ):
                    self.current_spot = -1
                    row = self.m.bg_spot_sample.center[0]
                    col = self.m.bg_spot_sample.center[1]

        self.imageview.crosshairs_c = col
        self.imageview.crosshairs_r = row
        self.imageview.hline.set_ydata( np.array([row,row]) )
        self.imageview.vline.set_xdata( np.array([col,col]) )

        # self.imageview.figure.canvas.draw()
#        self.imageview.myaxes.draw_artist(self.imageview.hline)
#        self.imageview.myaxes.draw_artist(self.imageview.vline)
        self.imageview.figure.canvas.draw()

        self.dataview_updater()
        self.spotInfo_updater()


    def main(self):
        self.show()

    def connectActions(self):
        """Connect the user interface controls to the logic """
        self.centralwidget.keyPressEvent = self.keyPressEvent
        self.selectDataDirPushButton.clicked.connect( self.selectDataDir )
        self.selectSPEComboBox.activated.connect( self.selectSPE )
        self.setBGSpotPushButton.clicked.connect( self.setBGSpot )
        self.addSignalSpotPushButton.clicked.connect( self.addSignalSpot )
        self.createSpotArrayPushButton.clicked.connect( self.createSpotArray )
        self.clearAllSpotsPushButton.clicked.connect( self.clearAllSpots )
#        self.initAnalysisPushButton.clicked.connect( self.initAnalysis )
        self.checkSpotValidityPushButton.clicked.connect( self.checkSpotValidity )
        self.cosineFitPushButton.clicked.connect( self.cosineFit )
        self.findModDepthsPushButton.clicked.connect( self.findModDepths )
        self.showStuffComboBox.activated.connect( self.showStuff )
        self.ETrulerPushButton.clicked.connect( self.ETruler )
        self.showWhatComboBox.activated.connect( self.dataview_updater )
        self.saveContrastImagesPushButton.clicked.connect( self.saveContrastImages_hdf5 )
        self.importSpotCoordsPushButton.clicked.connect( self.importSpotCoords )
        self.selectLastDataDirPushButton.clicked.connect( self.selectLastDataDir )
        self.addSMPushButton.clicked.connect( self.addSM )
        self.testPushButton.clicked.connect( self.test )
        self.clearSpotPushButton.clicked.connect( self.clearSpot )
        self.moveSpotUpPushButton.clicked.connect( self.moveSpotUp )
        self.moveSpotDownPushButton.clicked.connect( self.moveSpotDown )
        self.moveSpotLeftPushButton.clicked.connect( self.moveSpotLeft )
        self.moveSpotRightPushButton.clicked.connect( self.moveSpotRight )
        self.frameSlider.valueChanged.connect( self.frameSliderChanged )
        self.spotTypeComboBox.activated.connect( self.spotTypeChanged )

    def spotTypeChanged( self ):
        if self.spotTypeComboBox.currentText()=='square':
            self.spotAttributeLabel.setText('edge length')
        elif self.spotTypeComboBox.currentText()=='circle':
            self.spotAttributeLabel.setText('diameter')
        else:
            raise taxes

    def test(self):
        b = self.m.spots[0].dilate( borderwidth=2 )
        self.dataview.myaxes.cla()
        self.dataview.myaxes.imshow(b, interpolation='nearest')
        self.dataview.figure.canvas.draw()

    def tryToUpdateFile( self, basefilename, what ):
        filename = basefilename+'_'+what+'_image.txt'
        # if the file exists, try loading and updating it
        if os.path.isfile(filename):
            try: 
                sc = np.loadtxt(filename)
            except IOError:
                sc = getattr(self.m, what+'_image')
            # where the current image is not nan, update the old image
            new_value_indices = ~np.isnan( getattr(self.m, what+'_image') ) 
            sc[new_value_indices] = getattr(self.m, what+'_image')[new_value_indices]
            np.savetxt( filename, sc )
        else:
            np.savetxt( filename, getattr(self.m, what+'_image') )
        
    def saveContrastImages( self ):
        basefilename = self.data_directory + self.spefiles[self.selectSPEComboBox.currentIndex()][:-4]
        self.tryToUpdateFile( basefilename, 'spot_coverage' )
        self.tryToUpdateFile( basefilename, 'M_ex' )
        self.tryToUpdateFile( basefilename, 'M_em' )
        self.tryToUpdateFile( basefilename, 'phase_ex' )
        self.tryToUpdateFile( basefilename, 'phase_em' )
        self.tryToUpdateFile( basefilename, 'LS' )
        self.tryToUpdateFile( basefilename, 'ET_ruler' )
        self.tryToUpdateFile( basefilename, 'ET_model_md_fu' )
        self.tryToUpdateFile( basefilename, 'ET_model_th_fu' )
        self.tryToUpdateFile( basefilename, 'ET_model_gr' )
        self.tryToUpdateFile( basefilename, 'ET_model_et' )
        # np.savetxt( basefilename+'_M_em_image.txt', self.m.M_em_image )
        # np.savetxt( basefilename+'_phase_ex_image.txt', self.m.phase_ex_image )
        # np.savetxt( basefilename+'_phase_em_image.txt', self.m.phase_em_image )
        # np.savetxt( basefilename+'_LS_image.txt', self.m.LS_image )
        # np.savetxt( basefilename+'_ET_ruler_image.txt', self.m.ET_ruler_image )
        # np.savetxt( basefilename+'_ET_model_md_fu_image.txt', self.m.ET_model_md_fu_image )
        # np.savetxt( basefilename+'_ET_model_th_fu_image.txt', self.m.ET_model_th_fu_image )
        # np.savetxt( basefilename+'_ET_model_gr_image.txt', self.m.ET_model_gr_image )
        # np.savetxt( basefilename+'_ET_model_et_image.txt', self.m.ET_model_et_image )
        self.saveContrastImagesPushButton.setChecked(False)

    def saveContrastImages_hdf5( self ):
        outputdir = os.path.normpath( self.data_directory )
        print outputdir
        save_hdf5( movie=self.m, myspots=range(len(self.m.validspots)), \
                       fileprefix=outputdir, proc=0, images=True, spots=True )
        combine_outputs( basename=self.m.data_basename, fileprefix=outputdir )

    def ETruler(self):
        oldtext = self.setButtonWait( self.ETrulerPushButton )

        self.imageview.ET_ruler_rects = []
        self.m.ETrulerFFT()
        # ruler = []
        # for s in self.m.validspots:
        #     ruler.append(s.ET_ruler)
        for s in self.m.validspots:
            color = cm.jet(s.ET_ruler)
            r = Rectangle( (s.coords[0],s.coords[1]), s.width, s.height, \
                               facecolor=color, edgecolor=color, alpha=1, zorder=7 )
            self.imageview.ET_ruler_rects.append( r )

        self.setButtonWait( self.ETrulerPushButton, oldtext )


    def showStuff(self):
        what = self.showStuffComboBox.currentIndex()
        if what==0:  # spots
            self.imageview.show_stuff( what='spots' )
        elif what==1:  # M_ex            
            self.imageview.show_stuff( what='M_ex' )
#            self.imageview.cbar.
        elif what==2:  # M_em
            self.imageview.show_stuff( what='M_em' )
        elif what==3:  # phase_ex
            self.imageview.show_stuff( what='phase_ex' )
        elif what==4:  # phase_ex
            self.imageview.show_stuff( what='phase_em' )
        elif what==5:  # ET_ruler
            self.imageview.show_stuff( what='ET_ruler' )


    def findModDepths(self):
        oldtext = self.setButtonWait( self.findModDepthsPushButton )

        self.imageview.M_ex_rects = []
        self.imageview.M_em_rects = []
        self.imageview.phase_ex_rects = []
        self.imageview.phase_em_rects = []
        self.m.find_modulation_depths_and_phases_selective()
        # for s in self.m.validspots:
        #     color = cm.jet(s.M_ex)
        #     r = Rectangle( (s.coords[0],s.coords[1]), s.width, s.height, \
        #                        facecolor=color, edgecolor=color, alpha=1, zorder=7 )
        #     self.imageview.M_ex_rects.append( r )
        #     color = cm.jet(s.M_em)
        #     r = Rectangle( (s.coords[0],s.coords[1]), s.width, s.height, \
        #                        facecolor=color, edgecolor=color, alpha=1, zorder=7 )
        #     self.imageview.M_em_rects.append( r )
        #     color = cm.hsv(s.phase_ex/np.pi+.5)
        #     r = Rectangle( (s.coords[0],s.coords[1]), s.width, s.height, \
        #                        facecolor=color, edgecolor=color, alpha=1, zorder=7 )
        #     self.imageview.phase_ex_rects.append( r )
        #     color = cm.hsv(s.phase_em/np.pi+.5)
        #     r = Rectangle( (s.coords[0],s.coords[1]), s.width, s.height, \
        #                        facecolor=color, edgecolor=color, alpha=1, zorder=7 )
        #     self.imageview.phase_em_rects.append( r )
#            self.imageview.axes.add_patch( self.imageview.spot_rects[-1] )
#        self.imageview.show_stuff( what='M_ex' )
        self.imageview.show_stuff(what='M_ex')
        self.showStuffComboBox.setCurrentIndex(1)

        self.unsetButtonWait( self.findModDepthsPushButton, oldtext)

    def setButtonWait( self, button ):
        oldtext = button.text()
        button.setText('wait')
        button.setStyleSheet('QPushButton {color: red}')
        self.app.processEvents()
        return oldtext
        
    def unsetButtonWait( self, button, oldtext ):
        button.setText( oldtext )
        button.setStyleSheet('')
        self.app.processEvents()
        
    def cosineFit(self):
        oldtext = self.setButtonWait( self.cosineFitPushButton )
        # self.m.excitation_angles_grid = np.linspace( 0, np.pi, self.NanglesSpinBox.value() )
        # self.m.emission_angles_grid = np.linspace( 0, np.pi, self.NanglesSpinBox.value() )
        self.m.fit_all_portraits_spot_parallel_selective()
        for s in self.m.validspots:
            s.portrait_residuals( iportrait=0 )
        self.unsetButtonWait( self.cosineFitPushButton, oldtext )

    def checkSpotValidity(self):
        oldtext = self.setButtonWait( self.checkSpotValidityPushButton )

        self.m.are_spots_valid( SNR=self.SNRSpinBox.value() )
        # reset color back to red 
        for r in self.imageview.spot_gfxrepr:
            r.set_facecolor('red')
        # turn valid spots back on green
        for vsi in self.m.validspotindices:
            self.imageview.spot_gfxrepr[vsi].set_facecolor('green')
            self.imageview.spot_gfxrepr[vsi].set_edgecolor('green')
#        self.imageview.figure.canvas.draw()
        self.imageview.show_stuff(what='spots')

        self.unsetButtonWait( self.checkSpotValidityPushButton, oldtext )


    def initAnalysis(self):
        oldtext = self.setButtonWait( self.initAnalysisPushButton )
        self.m.collect_data()
        self.m.startstop()
        self.m.assign_portrait_data()
        self.initAnalysisPushButton.setChecked(False)
        self.unsetButtonWait( self.initAnalysisPushButton, oldtext )

    def addSM(self):
        nspotsbefore = len(self.m.spots)
        boxedgelength = np.int( self.spotEdgeLengthSpinBox.value() )

        r = np.int( self.imageview.crosshairs_r )
        c = np.int( self.imageview.crosshairs_c )

        spot_type=self.spotTypeComboBox.currentText()
        use_exspot=self.useExSpotCheckBox.isChecked()
        use_borderbg=self.useBorderBGCheckBox.isChecked()

        if spot_type=='square':
            ri  = np.int( np.floor(r) - (boxedgelength-1)/2 )
            ci  = np.int( np.floor(c) - (boxedgelength-1)/2 )
            ri2 = np.int( np.ceil(r) + (boxedgelength-1)/2 )
            ci2 = np.int( np.ceil(c) + (boxedgelength-1)/2 )
            shape = {'type': 'Rectangle', \
                         'left':  ci, \
                         'right': ci2, \
                         'lower': ri2, \
                         'upper': ri }
        elif spot_type=='circle':
            shape = {'type': 'Circle', \
                         'center': (r,c), \
                         'radius': boxedgelength}
        else:
            raise hell

        self.m.define_spot( shape, use_exspot=use_exspot, use_borderbg=use_borderbg )
        print "Defined",self.m.spots[-1]

        # now populate that list
        for s in self.m.spots[nspotsbefore:]:
            r = s.graphical_representation()
            self.imageview.spot_gfxrepr.append( r )

        # reset selection rectangle
        self.imageview.myrect.set_xy((0,0))
        self.imageview.myrect.set_width(0)
        self.imageview.myrect.set_height(0)

        self.imageview.figure.canvas.draw()
        self.imageview.show_stuff()

    def importSpotCoords(self):
        nspotsbefore  = len(self.m.spots)
        print "nspotsbefore: ",nspotsbefore
        import_spot_positions( self.m, self.spefiles[self.selectSPEComboBox.currentIndex()][:-4], \
                                   boxedgelength=self.spotEdgeLengthSpinBox.value(), \
                                   spot_type=self.spotTypeComboBox.currentText(), \
                                   use_exspot=self.useExSpotCheckBox.isChecked(), \
                                   use_borderbg=self.useBorderBGCheckBox.isChecked() )

        # now populate that list
        for s in self.m.spots[nspotsbefore:]:
            r = s.graphical_representation()
            self.imageview.spot_gfxrepr.append( r )

        # reset selection rectangle
        self.imageview.myrect.set_xy((0,0))
        self.imageview.myrect.set_width(0)
        self.imageview.myrect.set_height(0)

        self.imageview.figure.canvas.draw()
        self.imageview.show_stuff()
        self.app.processEvents()

    def createSpotArray(self):
        oldtext = self.setButtonWait( self.createSpotArrayPushButton )

        # how many spots do we have already?
        nspotsbefore = len(self.m.spots)

        resolution = np.int( self.spotEdgeLengthSpinBox.value() )

        fullbounds = [self.imageview.c0, self.imageview.r0, self.imageview.c1, self.imageview.r1]
        topedges = np.arange( np.int(fullbounds[1]), np.int(fullbounds[3]), resolution )  

        for line in topedges:
            for col in range( np.int(fullbounds[0]), np.int(fullbounds[2]), resolution ):
                bounds = [ col, line, col+resolution-1, line+resolution-1 ]
                self.m.define_spot( bounds )

        # now populate that list
        for s in self.m.spots[nspotsbefore:]:
            r = s.graphical_representation()
            self.imageview.spot_gfxrepr.append( r )

        # reset selection rectangle
        self.imageview.myrect.set_xy((0,0))
        self.imageview.myrect.set_width(0)
        self.imageview.myrect.set_height(0)

        self.imageview.figure.canvas.draw()
        self.imageview.show_stuff()

        self.unsetButtonWait( self.createSpotArrayPushButton, oldtext )

    def setBGSpot(self):
        shape = {'type':'Rectangle', \
                     'left':  np.int(self.imageview.c0), \
                     'right': np.int(self.imageview.c1), \
                     'upper': np.int(self.imageview.r0), \
                     'lower': np.int(self.imageview.r1) }

#        self.bg_region_coords = coords
        self.imageview.bg_rect.set_xy((self.imageview.c0,self.imageview.r0))
        self.imageview.bg_rect.set_width(self.imageview.c1-self.imageview.c0)
        self.imageview.bg_rect.set_height(self.imageview.r1-self.imageview.r0)
        # reset selection rectangle
        self.imageview.myrect.set_xy((0,0))
        self.imageview.myrect.set_width(0)
        self.imageview.myrect.set_height(0)
        # update canvas
        self.imageview.figure.canvas.draw()                       
        self.imageview.show_stuff()
        self.m.define_background_spot( shape, intensity_type='mean' )
        print "BG spot defined."
        print "BG spot: intensity: ",self.m.bg_spot_sample.intensity
        print "BG spot: std      : ",self.m.bg_spot_sample.std

    def addSignalSpot(self):
        nspotsbefore  = len(self.m.spots)

        shape = {'type':'Rectangle', \
                     'left':  np.int(self.imageview.c0), \
                     'right': np.int(self.imageview.c1), \
                     'upper': np.int(self.imageview.r0), \
                     'lower': np.int(self.imageview.r1) }

        # create spot
        s = self.m.define_spot( shape, intensity_type='mean' )

        # now populate that list
        for s in self.m.spots[nspotsbefore:]:
            r = s.graphical_representation()
            self.imageview.spot_gfxrepr.append( r )

        # reset selection rectangle
        self.imageview.myrect.set_xy((0,0))
        self.imageview.myrect.set_width(0)
        self.imageview.myrect.set_height(0)
        # update canvas
        self.imageview.figure.canvas.draw()

    def clearSpot(self):
        if not self.current_spot==None:
            self.imageview.spot_gfxrepr[ self.current_spot ].remove()
            self.imageview.spot_gfxrepr.pop( self.current_spot )
            self.m.spots.pop( self.current_spot )
            self.validspots = None
            self.m.initContrastImages()
            self.imageview.figure.canvas.draw()

    def clearAllSpots(self):
        for r in self.imageview.spot_gfxrepr:
            print r.get_axes()
            r.remove()
        self.imageview.spot_gfxrepr = []
        self.m.spots = []
        self.validspots = None
        self.m.initContrastImages()
        self.imageview.figure.canvas.draw()


    def moveSpot(self, direction, dist):
        if not self.current_spot==None:
            # grab coords
            coords = self.m.spots[self.current_spot].coords
            # clear old spot
            self.clearSpot()
            # move spot
            if direction=='up':
                coords[1] -= dist
                coords[3] -= dist
            elif direction=='down':
                coords[1] += dist
                coords[3] += dist
            elif direction=='left':
                coords[0] -= dist
                coords[2] -= dist
            elif direction=='right':
                coords[0] += dist
                coords[2] += dist
            else:
                raise ValueError('moveSpot: bad direction "%s", should be up|down|left|right.' % direction)

            # check for boundary cases
            if coords[0] < 0:
                coords[0] -= coords[0]
                coords[2] -= coords[0]
            if coords[2] >= self.m.sample_data.datasize[2]:
                coords[0] -= coords[2] -self.m.sample_data.datasize[2] +1
                coords[2] -= coords[2] -self.m.sample_data.datasize[2] +1
            if coords[1] < 0:
                coords[1] -= coords[1]
                coords[3] -= coords[1]
            if coords[3] >= self.m.sample_data.datasize[1]:
                coords[0] -= coords[3] -self.m.sample_data.datasize[3] +1
                coords[3] -= coords[3] -self.m.sample_data.datasize[3] +1


            # re-define spot
            s = self.m.define_spot( coords, intensity_type='mean' )
            # the previous 'current spot' doesn't exist anymore (the index will point to
            # some other spot, or may be invalid) -- the moved spot is at the end
            # of the spot list:
            self.current_spot = len(self.m.spots)-1
            # keep the crosshairs on said spot                        
            x = s.coords[0]+.5*s.width
            y = s.coords[1]+.5*s.height
            self.imageview.hline.set_ydata( np.array([y,y]) )
            self.imageview.vline.set_xdata( np.array([x,x]) )

            # recreate rectangle patch and update canvas
            r = Rectangle( (s.coords[0],s.coords[1]), s.width, s.height, \
                               facecolor='red', edgecolor='red', alpha=.3, zorder=8 )
            self.imageview.spot_gfxrepr.append( r )
            self.imageview.myaxes.add_collection( r )
            self.imageview.figure.canvas.draw()

            # for some reason this needs to happen after re-drawing the canvas
            self.imageview.myaxes.draw_artist(self.imageview.hline)
            self.imageview.myaxes.draw_artist(self.imageview.vline)

            self.dataview_updater()
            self.spotInfo_updater()


    def moveSpotUp(self):
        self.moveSpot( direction='up', dist=1 )
            
    def moveSpotDown(self):
        self.moveSpot( direction='down', dist=1 )

    def moveSpotLeft(self):
        self.moveSpot( direction='left', dist=1 )

    def moveSpotRight(self):
        self.moveSpot( direction='right', dist=1 )

    def dataview_updater(self):
        self.dataview.figure.clear()
        self.dataview.myaxes = self.dataview.figure.add_subplot(111)

        if not self.current_spot==None:
            if self.current_spot==-1:
                spot = self.m.bg_spot_sample
            else:
                spot = self.m.spots[self.current_spot]

            showWhat = self.showWhatComboBox.currentIndex()
            if showWhat==0:    # intensity trace
                if hasattr(spot, 'intensity'):
                    # self.dataview.myaxes.cla()
                    # self.dataview.figure.canvas.draw()
                    self.dataview.myaxes.plot( spot.intensity, 'bx-' )
                    if hasattr(spot,'borderbg'):
                        self.dataview.myaxes.plot( spot.borderbg, 'gx-' )
                    self.dataview.myaxes.set_aspect('normal')

            elif showWhat==1:    # portrait data
                # self.dataview.myaxes.cla()
                # self.dataview.myaxes.set_ylim(0,1)
                # self.dataview.myaxes.set_ylim(auto=True)

                Nemangles = np.unique( self.m.emangles ).size
                scaler = 180.0/Nemangles/2/np.max( spot.intensity )

                for p in range(self.m.Nportraits):
                    for l in range(self.m.Nlines):
                        exangles, emangle = self.m.retrieve_angles( p, l )
#                        print "emangle: ",emangle
                        intensity = spot.retrieve_intensity( p, l )
#                        print intensity
#                        print exangles
                        sortind = np.argsort( exangles )
                        self.dataview.myaxes.plot( exangles[sortind], \
                                                       180/np.pi*emangle + scaler*intensity[sortind], \
                                                       'bx:' )
                        self.dataview.myaxes.fill_between( exangles[sortind], \
                                                               180/np.pi*emangle + scaler*intensity[sortind], \
                                                               180/np.pi*emangle*np.ones_like(intensity[sortind]), \
                                                               facecolor='b', alpha=.4 )
                        
                        if hasattr( spot, 'linefitparams' ):
                            fit = spot.retrieve_line_fit( p, l, exangles )
                            self.dataview.myaxes.plot( exangles[sortind], \
                                                         180/np.pi*emangle + scaler*fit[sortind], \
                                                         'r-' )            

            elif showWhat==2:    # portrait fit
                i=self.dataview.myaxes.imshow( spot.recover_average_portrait_matrix(),\
                                                 origin='lower', interpolation='nearest' )
                self.dataview.myaxes.set_position( [.25,.1,.62,.62] )
                axis_proj_em = self.dataview.fig.add_axes([.05,.1,.2,.62])
                axis_proj_ex = self.dataview.fig.add_axes([.25,.1+.62,.62,.2])

                axis_proj_em.plot( spot.proj_em, self.m.emission_angles_grid )
                axis_proj_em.set_xlim( xmin=np.max(axis_proj_em.get_xlim()), xmax=0 )
                [ label.set_rotation(90) for label in axis_proj_em.get_xticklabels() ]
                axis_proj_ex.plot( self.m.excitation_angles_grid, spot.proj_ex )
                axis_proj_ex.set_ylim( ymin=0, ymax=np.max(axis_proj_ex.get_ylim()) )

        else:
            r = self.imageview.crosshairs_r
            c = self.imageview.crosshairs_c
            self.dataview.myaxes.plot( self.m.sample_data.rawdata[ :, r, c ] )
            frame = self.frameSlider.value()
            self.dataview.myaxes.axvline( frame, color='r', ls=':', zorder=9 )

        self.dataview.figure.canvas.draw()


    def spotInfo_updater(self):
        if not self.current_spot==None:

            if self.current_spot==-1:
                s = self.m.bg_spot_sample
            else:
                s = self.m.spots[self.current_spot]

            # is this a valid spot?
            isvalid = False
            if hasattr(self.m, 'validspotindices'):
                if self.m.validspotindices.count( self.current_spot )>0:
                    isvalid = True

            # prepare info
            # infostring  = "area=( %d--%d, %d--%d )<br>" % (s.coords[0],s.coords[0]+s.width, \
            #                                                s.coords[1],s.coords[1]+s.height)
            infostring  = "shape=%s<br>" % (str(s.shape))
            if hasattr(s,'SNR'):
                if isvalid:
                    infostring += 'mean SNR   = <font color="#00df00">%f</font>    <br>' % s.meanSNR
                else:
                    infostring += 'mean SNR   = <font color="#ef0000">%f</font>    <br>' % s.meanSNR

#            print s.verticalfitparams
#            print s.linefitparams

            if hasattr(s,'M_ex'):
                infostring += 'M<sub>ex</sub>  = %f<br>' % s.M_ex
            if hasattr(s,'M_em'):
                infostring += 'M<sub>em</sub>  = %f<br>' % s.M_em
            if hasattr(s,'phase_ex'):
                infostring += '&phi;<sub>ex</sub>  = %f deg<br>' % (s.phase_ex*180/np.pi)
            if hasattr(s,'phase_em'):
                infostring += '&phi;<sub>em</sub>  = %f deg<br>' % (s.phase_em*180/np.pi)
            if hasattr(s,'residual_full'):
                infostring += 'resi<sub>full</sub>  = %f <br>' % (s.residual_full)
            if hasattr(s,'ET_ruler'):
                infostring += 'ET<sub>ruler</sub>  = %f<br>' % s.ET_ruler
            if hasattr(s,'ET_model_md_fu'):
                infostring += 'ET model: M<sub>funnel</sub>  = %f<br>' % s.ET_model_md_fu
                infostring += 'ET model: &theta;<sub>funnel</sub>  = %f<br>' % (s.ET_model_th_fu*180/np.pi)
                infostring += 'ET model: geom. ratio  = %f<br>' % s.ET_model_gr
                infostring += 'ET model: ET param.  = %f<br>' % s.ET_model_et

            if hasattr(s,'intensity'):
                infostring += 'mean(intensity)= %f<br>' % np.mean(s.intensity)
                infostring += 'min(intensity) = %f<br>' % np.min(s.intensity)
                infostring += 'max(intensity) = %f<br>' % np.max(s.intensity)

            self.spotInfoTextBrowser.setHtml(infostring)


    def selectSPE(self):
        print self.selectSPEComboBox.currentIndex()
        self.load_and_display_spe_file(fileindex=self.selectSPEComboBox.currentIndex())


    def setPhaseOffset(self):
        self.phaseOffset = self.phaseOffsetLineEdit.text().toDouble()[0]
        print "Phase offset set to %f [deg]" % self.phaseOffset

    def selectLastDataDir(self):
        # get from lastdir
        f = open('lastdir.txt','r')
        fname = f.readline()
        f.close()

        self.data_directory = os.path.normpath( str(fname) ) + os.path.sep
        self.dataDirLabel.setText(self.data_directory)

        # go to dir
        os.chdir( self.data_directory )

        self.handleDataDir()

        
    def selectDataDir(self):
        fname = QtGui.QFileDialog.getExistingDirectory(self, 'Show me where the data is', \
                directory='' )
        print "Selected data directory ", fname

        if fname=='':
            return

        self.data_directory = os.path.normpath( str(fname) ) + os.path.sep
        self.dataDirLabel.setText(self.data_directory)

        # save in lastdir
        f = open('lastdir.txt','w')
        f.writelines( [self.data_directory] )
        f.close()

        # go to dir
        os.chdir( self.data_directory )


        self.handleDataDir()

    def handleDataDir(self):

        print "Looking for SPE data..."
        # get all filenames
        for file in os.listdir("."):
            # get all that end in spe
            if file.endswith(".spe") or file.endswith(".SPE") or file.endswith(".npy"):
                # if we don't have it already
                if self.spefiles.count(file)==0:
                    # add it to the list
                    self.spefiles.append(file)
                    print "File %s added to the list." % file
                else:
                    print "File %s is already in the list." % file

        self.spefiles.sort()

        self.motorfiles = ['']*len(self.spefiles)

        print "Looking for corresponding motor data..."
        for file in os.listdir("."):
            if file.endswith(".txt"):
                # do we know of the corresponding spe file?
                match = [file[:-4]=='MS-'+spe[:-4] for spe in self.spefiles]
                if any( match ):
                    where = np.nonzero(match)[0][0]
                    self.motorfiles[where] = file
                    print "Found motorfile matching spe-file basename: %s" % file

        print "Any loose ends?"
        for i,file in enumerate(self.motorfiles):
            if len(file)==0:                
                print "SPE file %s doesn't have a motor file... removed from list." % self.spefiles[i]
                self.spefiles.pop(i)
                self.motorfiles.pop(i)   # is an empty string, but necessary to prune the list
        print "All good."
        
        self.selectSPEComboBox.clear()
        for s in self.spefiles:
            self.selectSPEComboBox.addItem(s)
        self.selectSPEComboBox.setCurrentIndex(0)

        if not len(self.spefiles)==0:
            self.load_and_display_spe_file(0)

    def load_and_display_spe_file(self,fileindex=0):
        del(self.m)
        # from guppy import hpy; h=hpy()
        # w=h.heap()
        # print w
        print "loading file %s ... " % (self.spefiles[fileindex]),
        sys.stdout.flush()        
        # self.m = Movie( self.data_directory+"/"+self.spefiles[fileindex], \
        #                     self.data_directory+"/"+self.motorfiles[fileindex], \
        #                     phase_offset_excitation=self.phase_offset*np.pi/180.0, \
        #                     which_setup=self.setup_list[self.which_setup], \
        #                     excitation_optical_element=self.optical_element)
        # self.m = Movie( self.data_directory, \
        #                     self.spefiles[fileindex][:-4], \
        #                     phase_offset_in_deg=self.phase_offset, \
        #                     which_setup=self.setup_list[self.which_setup], \
        #                     excitation_optical_element=self.optical_element)
        self.m = Movie( self.data_directory, \
                            self.spefiles[fileindex][:-4] )
#                            phase_offset_in_deg=self.phase_offset )
#        self.m.startstop()
        self.m.find_portraits()
        self.m.find_lines()

        print self.imageview.myaxes.artists
        print self.imageview.myaxes.patch
        print self.imageview.myaxes.collections
        print self.imageview.myaxes.images
        print self.imageview.myaxes.lines
        print self.imageview.myaxes.patches
        print self.imageview.myaxes.texts
        print self.imageview.myaxes.xaxis
        print self.imageview.myaxes.yaxis
        print '-----------------------------------------'

        self.imageview.show_image( self.m.sample_data.rawdata[0,:,:], zorder=1, cmap=cm.gray )
#        self.imageview.myaxes.imshow( 
#        self.imageview.draw()
        print "done"
        print self.imageview.myaxes.artists
        print self.imageview.myaxes.patch
        print self.imageview.myaxes.collections
        print self.imageview.myaxes.images
        print self.imageview.myaxes.lines
        print self.imageview.myaxes.patches
        print self.imageview.myaxes.texts
        print self.imageview.myaxes.xaxis
        print self.imageview.myaxes.yaxis

        self.frameSlider.setMaximum( self.m.sample_data.rawdata.shape[0]-1 )


    def frameSliderChanged( self ):
        frame = self.frameSlider.value()
        print frame
        self.imageview.im.set_data( self.m.sample_data.rawdata[frame,:,:] )
        cmin = np.min( self.m.sample_data.rawdata[frame,:,:] )
        cmax = np.max( self.m.sample_data.rawdata[frame,:,:] )
        self.imageview.im.set_clim(cmin,cmax)
        self.imageview.figure.canvas.draw()
        #self.imageview.show_image( self.m.sample_data.rawdata[frame,:,:], zorder=1, cmap=cm.gray )

        self.dataview_updater()





if __name__=='__main__':
    app = QtGui.QApplication(sys.argv)
    instance = the2dlogic(app=app)
    instance.main()
    sys.exit(app.exec_())
