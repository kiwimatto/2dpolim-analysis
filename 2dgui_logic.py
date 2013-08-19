from PyQt4 import QtCore,QtGui
import sys, os, time
import util2dpolim.gui.the2dgui as the2dgui
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

        self.spefiles = []
        self.m = None
        self.pwd = os.path.dirname(os.path.abspath(__file__))
        self.optical_element = 'Polarizer'
        self.which_setup = self.whichSetupComboBox.currentIndex()
        self.setup_list = ['old setup','new setup','cool new setup','header']
        self.phase_offset = self.phaseOffsetLineEdit.text().toDouble()[0]        
        self.current_spot = None
        self.app = app

    def keyPressEvent(self, event):
        if event.key()==QtCore.Qt.Key_Up:
            self.move_crosshairs('up')
        if event.key()==QtCore.Qt.Key_Down:
            self.move_crosshairs('down')
        if event.key()==QtCore.Qt.Key_Left:
            self.move_crosshairs('left')
        if event.key()==QtCore.Qt.Key_Right:
            self.move_crosshairs('right')

    def move_crosshairs( self, direction ):
        x = self.imageview.crosshairs_x
        y = self.imageview.crosshairs_y
        if not self.current_spot==None:
            if direction=='up' or direction=='down':
                delta = self.m.spots[self.current_spot].height
            elif direction=='left' or direction=='right':
                delta = self.m.spots[self.current_spot].width
            else:
                print 'whoa?'
        else:
            delta = 1

        if direction=='up':
            y -= delta
        elif direction=='down':
            y += delta
        elif direction=='left':
            x -= delta
        elif direction=='right':
            x += delta
        x = np.clip(x,0,self.m.camera_data.datasize[1])
        y = np.clip(y,0,self.m.camera_data.datasize[2])
        self.imageview.crosshairs_x = x
        self.imageview.crosshairs_y = y
        self.crosshair_pick()

    def crosshair_pick(self):
        x = self.imageview.crosshairs_x
        y = self.imageview.crosshairs_y
        print 'x0=%f\ty0=%f' % (x,y)
        sys.stdout.flush()        
        
        self.current_spot = None
        if (not self.m==None) and hasattr(self.m,'spots'):
            # find the spot which covers these coordinates
            for si,s in enumerate(self.m.spots):
                if (s.coords[0] <= x) and (s.coords[0]+s.width > x):
                    if (s.coords[1] <= y) and (s.coords[1]+s.height > y):
                        self.current_spot = si
                        x = s.coords[0]+.5*s.width
                        y = s.coords[1]+.5*s.height
                        print 'picked spot #%d' % si
                        sys.stdout.flush()
                        break
        self.imageview.crosshairs_x = x
        self.imageview.crosshairs_y = y
        self.imageview.hline.set_ydata( np.array([y,y]) )
        self.imageview.vline.set_xdata( np.array([x,x]) )

        self.imageview.figure.canvas.restore_region(self.imageview.blitbackground)
        # self.imageview.figure.canvas.draw()
        self.imageview.axes.draw_artist(self.imageview.hline)
        self.imageview.axes.draw_artist(self.imageview.vline)
        self.imageview.figure.canvas.blit(self.imageview.axes.bbox)

        self.dataview_updater()
        self.spotInfo_updater()

    def main(self):
        self.show()

    def connectActions(self):
        """Connect the user interface controls to the logic """
        self.selectDataDirPushButton.clicked.connect( self.selectDataDir )
        self.phaseOffsetLineEdit.editingFinished.connect( self.setPhaseOffset )
        self.whichSetupComboBox.activated.connect( self.setWhichSetup )
        self.selectSPEComboBox.activated.connect( self.selectSPE )
        self.setBGSpotPushButton.clicked.connect( self.setBGSpot )
        self.addSignalSpotPushButton.clicked.connect( self.addSignalSpot )
        self.createSpotArrayPushButton.clicked.connect( self.createSpotArray )
        self.clearAllSpotsPushButton.clicked.connect( self.clearAllSpots )
        self.initAnalysisPushButton.clicked.connect( self.initAnalysis )
        self.checkSpotValidityPushButton.clicked.connect( self.checkSpotValidity )
        self.centralwidget.keyPressEvent = self.keyPressEvent
        self.cosineFitPushButton.clicked.connect( self.cosineFit )
        self.findModDepthsPushButton.clicked.connect( self.findModDepths )
        self.showStuffComboBox.activated.connect( self.showStuff )
        self.ETrulerPushButton.clicked.connect( self.ETruler )
        self.showWhatComboBox.activated.connect( self.dataview_updater )
        self.saveContrastImagesPushButton.clicked.connect( self.saveContrastImages_hdf5 )
        self.toolButton1.clicked.connect( self.tool1 )
        self.toolButton2.clicked.connect( self.tool2 )
        self.toolButton3.clicked.connect( self.tool3 )

    def tool3(self):
        self.tool2()
        self.findModDepths()
        time.sleep(.1)

    def tool2(self):
        self.tool1()
        self.cosineFit()
        time.sleep(.1)

    def tool1(self):
        self.initAnalysis()
        time.sleep(.1)
        self.checkSpotValidity()
        time.sleep(.1)

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
        basefilename = self.data_directory + '/' + self.spefiles[self.selectSPEComboBox.currentIndex()][:-4]
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
        outputdir = os.path.normpath( self.data_directory+'/' )
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
        self.m.find_modulation_depths_and_phases()
        for s in self.m.validspots:
            color = cm.jet(s.M_ex)
            r = Rectangle( (s.coords[0],s.coords[1]), s.width, s.height, \
                               facecolor=color, edgecolor=color, alpha=1, zorder=7 )
            self.imageview.M_ex_rects.append( r )

            color = cm.jet(s.M_em)
            r = Rectangle( (s.coords[0],s.coords[1]), s.width, s.height, \
                               facecolor=color, edgecolor=color, alpha=1, zorder=7 )
            self.imageview.M_em_rects.append( r )

            print s.phase_ex/np.pi+.5
            color = cm.jet(s.phase_ex/np.pi+.5)
            r = Rectangle( (s.coords[0],s.coords[1]), s.width, s.height, \
                               facecolor=color, edgecolor=color, alpha=1, zorder=7 )
            self.imageview.phase_ex_rects.append( r )

            color = cm.hsv(s.phase_ex/np.pi+.5)
            r = Rectangle( (s.coords[0],s.coords[1]), s.width, s.height, \
                               facecolor=color, edgecolor=color, alpha=1, zorder=7 )
            self.imageview.phase_ex_rects.append( r )
            color = cm.hsv(s.phase_em/np.pi+.5)
            r = Rectangle( (s.coords[0],s.coords[1]), s.width, s.height, \
                               facecolor=color, edgecolor=color, alpha=1, zorder=7 )
            self.imageview.phase_em_rects.append( r )

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
        self.m.excitation_angles_grid = np.linspace( 0, np.pi, self.NanglesSpinBox.value() )
        self.m.emission_angles_grid = np.linspace( 0, np.pi, self.NanglesSpinBox.value() )
        self.m.fit_all_portraits_spot_parallel()
        self.unsetButtonWait( self.cosineFitPushButton, oldtext )

    def checkSpotValidity(self):
        oldtext = self.setButtonWait( self.checkSpotValidityPushButton )

        self.m.are_spots_valid( SNR=self.SNRSpinBox.value() )
        # reset color back to red 
        for r in self.imageview.spot_rects:
            r.set_facecolor('red')
        # turn valid spots back on green
        for vsi in self.m.validspotindices:
            self.imageview.spot_rects[vsi].set_facecolor('green')
            self.imageview.spot_rects[vsi].set_edgecolor('green')
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


    def createSpotArray(self):
        oldtext = self.setButtonWait( self.createSpotArrayPushButton )

        res = self.spotEdgeLengthSpinBox.value()
        coords = np.round( np.array( [self.imageview.x0, self.imageview.y0, \
                                          self.imageview.x1, self.imageview.y1] ) ).astype(np.int)
        if coords[0]>coords[2]:
            coords[0],coords[2]=coords[2],coords[0]
        if coords[1]>coords[3]:
            coords[1],coords[3]=coords[3],coords[1]

        # how many spots do we have already?
        nspotsbefore = len(self.m.spots)
        print nspotsbefore

        # now add spots
        grid_image_section_into_squares_and_define_spots( self.m, res, coords )

        # now populate that list
        for s in self.m.spots[nspotsbefore:]:
            r = Rectangle( (s.coords[0],s.coords[1]), s.width, s.height, \
                               facecolor='red', edgecolor='red', alpha=.3, zorder=7 )
            self.imageview.spot_rects.append( r )
            self.imageview.axes.add_patch( self.imageview.spot_rects[-1] )

        # reset selection rectangle
        self.imageview.rect.set_xy((0,0))
        self.imageview.rect.set_width(0)
        self.imageview.rect.set_height(0)

        self.imageview.figure.canvas.draw()
        self.imageview.show_stuff()

        self.unsetButtonWait( self.createSpotArrayPushButton, oldtext )

    def setBGSpot(self):
        coords = np.round( np.array( [self.imageview.x0, self.imageview.y0, \
                                          self.imageview.x1, self.imageview.y1] ) ).astype(np.int)
        if coords[0]>coords[2]:
            coords[0],coords[2]=coords[2],coords[0]
        if coords[1]>coords[3]:
            coords[1],coords[3]=coords[3],coords[1]

#        self.bg_region_coords = coords
        self.imageview.bg_rect.set_xy((coords[0],coords[1]))
        self.imageview.bg_rect.set_width(coords[2]-coords[0])
        self.imageview.bg_rect.set_height(coords[3]-coords[1])
        # reset selection rectangle
        self.imageview.rect.set_xy((0,0))
        self.imageview.rect.set_width(0)
        self.imageview.rect.set_height(0)
        # update canvas
        self.imageview.figure.canvas.draw()                       
        self.imageview.show_stuff()
        self.m.define_background_spot( coords, intensity_type='mean' )


    def addSignalSpot(self):
        coords = np.round( np.array( [self.imageview.x0, self.imageview.y0, \
                                          self.imageview.x1, self.imageview.y1] ) ).astype(np.int)
        if coords[0]>coords[2]:
            coords[0],coords[2]=coords[2],coords[0]
        if coords[1]>coords[3]:
            coords[1],coords[3]=coords[3],coords[1]

        r = Rectangle( (coords[0],coords[1]), coords[2]-coords[0], \
                           coords[3]-coords[1], facecolor='green', edgecolor='green', alpha=.3, zorder=8 )
        self.imageview.spot_rects.append( r )
        self.imageview.axes.add_patch( self.imageview.spot_rects[-1] )

        # reset selection rectangle
        self.imageview.rect.set_xy((0,0))
        self.imageview.rect.set_width(0)
        self.imageview.rect.set_height(0)
        # update canvas
        self.imageview.figure.canvas.draw()
        self.m.define_spot( coords, intensity_type='mean' )


    def clearAllSpots(self):
        for r in self.imageview.spot_rects:
            r.remove()
        self.imageview.spot_rects = []
        self.m.spots = []
        self.m.initContrastImages()
        self.imageview.figure.canvas.draw()


    def dataview_updater(self):
        if not self.current_spot==None:
            showWhat = self.showWhatComboBox.currentIndex()
            if showWhat==0:    # intensity trace
                if hasattr(self.m.spots[self.current_spot], 'intensity'):
#                    print 'showing intensity trace'
                    self.dataview.clear()
                    self.dataview.figure.canvas.draw()
                    self.dataview.axes.plot( self.m.spots[self.current_spot].intensity, 'bx-' )
                    self.dataview.figure.canvas.draw()
            elif showWhat==1:    # portrait data
                if hasattr(self.m.spots[self.current_spot], 'portraits'):
 #                   print 'showing portrait data'
                    self.dataview.clear()
                    self.dataview.figure.canvas.draw()
                    # collect all intensities

                    Nemangles = np.unique(self.m.spots[self.current_spot].portraits[0].emangles).size
                    print self.m.spots[self.current_spot].portraits[0].emangles
                    allints = np.array([])
                    for l in self.m.spots[self.current_spot].portraits[0].lines:
                        allints = np.hstack( (allints,l.intensities) )
                    maxint = np.max(allints)
                    scaler = 180.0/Nemangles/2/maxint
                    
                    for p in self.m.spots[self.current_spot].portraits:
                        for l in p.lines:
                            self.dataview.axes.plot( l.exangles, 180/np.pi*l.emangle + scaler*l.intensities, 'bx:' )
                            self.dataview.axes.fill_between( l.exangles, \
                                     180/np.pi*l.emangle + scaler*l.intensities, \
                                     180/np.pi*l.emangle*np.ones_like(l.intensities), \
                                     facecolor='b', alpha=.4 )
                            print np.min(l.intensities),' --- ',np.max(l.intensities)

                            if l.was_fitted:
                                self.dataview.axes.plot( l.exangles, 180/np.pi*l.emangle + scaler*l.cosValue(l.exangles), 'r-' )
                                print 'phase=%f\tI0=%f\tM0=%f\tresi=%f' % (l.phase,l.I0,l.M0,l.resi)
                                print 'l.exangles=', str(l.exangles)

                    self.dataview.figure.canvas.draw()

            elif showWhat==2:    # portrait fit
                print 'yup'
                self.dataview.clear()
                self.dataview.figure.canvas.draw()
                i=self.dataview.axes.imshow( self.m.spots[self.current_spot].recover_average_portrait_matrix(),\
                                                 origin='lower', interpolation='nearest' )
                self.dataview.axes.set_position( [.25,.1,.62,.62] )
                axis_proj_em = self.dataview.fig.add_axes([.05,.1,.2,.62])
                axis_proj_ex = self.dataview.fig.add_axes([.25,.1+.62,.62,.2])

                axis_proj_em.plot( self.m.spots[self.current_spot].proj_em, self.m.emission_angles_grid )
                axis_proj_em.set_xlim( xmin=np.max(axis_proj_em.get_xlim()), xmax=0 )
                [ label.set_rotation(90) for label in axis_proj_em.get_xticklabels() ]
                axis_proj_ex.plot( self.m.excitation_angles_grid, self.m.spots[self.current_spot].proj_ex )
                axis_proj_ex.set_ylim( ymin=0, ymax=np.max(axis_proj_ex.get_ylim()) )

                self.dataview.figure.canvas.draw()


#            print dir(self.m.spots[self.current_spot])
        else:
            self.dataview.clear()
            self.dataview.figure.canvas.draw()


    def spotInfo_updater(self):
        # is this a valid spot?
        isvalid = False
        s = self.m.spots[self.current_spot]
        if hasattr(self.m, 'validspotindices'):
            if self.m.validspotindices.count( self.current_spot )>0:
                isvalid = True

        # prepare info
        infostring  = "area=( %d--%d, %d--%d )<br>" % (s.coords[0],s.coords[0]+s.width, \
                                                           s.coords[1],s.coords[1]+s.height)
        if hasattr(s,'SNR'):
            if isvalid:
                infostring += 'SNR   = <font color="#00ff00">%f</font>    <br>' % s.SNR
            else:
                infostring += 'SNR   = <font color="#ff0000">%f</font>    <br>' % s.SNR

        if hasattr(s,'M_ex'):
            infostring += 'M<sub>ex</sub>  = %f<br>' % s.M_ex
        if hasattr(s,'M_em'):
            infostring += 'M<sub>em</sub>  = %f<br>' % s.M_em
        if hasattr(s,'phase_ex'):
            infostring += '&phi;<sub>ex</sub>  = %f deg<br>' % (s.phase_ex*180/np.pi)
        if hasattr(s,'phase_em'):
            infostring += '&phi;<sub>em</sub>  = %f deg<br>' % (s.phase_em*180/np.pi)
        if hasattr(s,'ET_ruler'):
            infostring += 'ET<sub>ruler</sub>  = %f<br>' % s.ET_ruler
        if hasattr(s,'ET_model_md_fu'):
            infostring += 'ET model: M<sub>funnel</sub>  = %f<br>' % s.ET_model_md_fu
            infostring += 'ET model: &theta;<sub>funnel</sub>  = %f<br>' % (s.ET_model_th_fu*180/np.pi)
            infostring += 'ET model: geom. ratio  = %f<br>' % s.ET_model_gr
            infostring += 'ET model: ET param.  = %f<br>' % s.ET_model_et

        if hasattr(s,'intensity'):
            infostring += 'min(intensity) = %f<br>' % np.min(s.intensity)
            infostring += 'max(intensity) = %f<br>' % np.max(s.intensity)

        self.spotInfoTextBrowser.setHtml(infostring)


    def selectSPE(self):
        print self.selectSPEComboBox.currentIndex()
        self.load_and_display_spe_file(fileindex=self.selectSPEComboBox.currentIndex())


    def setWhichSetup(self):
        self.which_setup = self.whichSetupComboBox.currentIndex()        

    def setPhaseOffset(self):
        self.phaseOffset = self.phaseOffsetLineEdit.text().toDouble()[0]
        print "Phase offset set to %f [deg]" % self.phaseOffset

    def selectDataDir(self):
        fname = QtGui.QFileDialog.getExistingDirectory(self, 'Show me where the data is', \
                directory='' )
        print "Selected data directory ", fname

        self.data_directory = str(fname)
        self.selectDataDirPushButton.setText(self.data_directory)

        # go to dir
        os.chdir( self.data_directory )
        print "Looking for SPE data..."
        # get all filenames
        for file in os.listdir("."):
            # get all that end in spe
            if file.endswith(".spe") or file.endswith(".SPE"):
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
        self.m = Movie( self.data_directory, \
                            self.spefiles[fileindex][:-4], \
                            phase_offset_in_deg=self.phase_offset, \
                            which_setup=self.setup_list[self.which_setup], \
                            excitation_optical_element=self.optical_element)

        self.imageview.show_image( self.m.camera_data.rawdata[0,:,:], zorder=1, cmap=cm.gray )
#        self.imageview.axes.imshow( 
#        self.imageview.draw()
        print "done"


if __name__=='__main__':
    app = QtGui.QApplication(sys.argv)
    instance = the2dlogic(app=app)
    instance.main()
    sys.exit(app.exec_())
