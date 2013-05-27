#!/usr/bin/env python

# embedding_in_qt4.py --- Simple Qt4 application embedding matplotlib canvases
#
# Copyright (C) 2005 Florent Rougon
#               2006 Darren Dale
#
# This file is an example program for matplotlib. It may be used and
# modified with no restriction; raw copies as well as modified versions
# may be distributed without limitation.

import sys, os
from PyQt4 import QtGui, QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import matplotlib.cm as cmap

import numpy as np
#from pyspec.ccd.files import PrincetonSPEFile

from util_2d import *
import spot_picker

class MyStaticMplCanvas(FigureCanvas):
    """Simple canvas with a sine plot."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        i = self.axes.imshow(np.outer( np.linspace(0,1,10),np.linspace(0,2,10) ), zorder=1 )
        self.cbar = self.fig.colorbar(i)
        # We want the axes cleared every time plot() is called
#        self.axes.hold(False)

        #
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

        self.crosshairs_x = 0
        self.crosshairs_y = 0

        self.bg_rect = Rectangle((0,0), 0, 0, facecolor='blue', edgecolor='blue', alpha=.3, \
                                     zorder=8 )
        self.axes.add_patch(self.bg_rect)

        self.signal_rect = Rectangle((0,0), 0, 0, facecolor='red', edgecolor='red', alpha=.3, \
                                     zorder=9 )
        self.axes.add_patch(self.signal_rect)

        self.anno = spot_picker.Annotate(self.axes)
        self.draw()

    def clear(self):
        self.fig.clear()
#        self.axes = self.fig.add_subplot(111)

    def show_portrait(self,portrait):
        self.fig.clear()
        self.axes = self.fig.add_subplot(111)
        self.axes.imshow( portrait, origin='bottom' )
        self.draw()
        

    # def update_figure(self, data):
    #     self.data = data
    #     self.fig.clear()
    #     self.axes = self.fig.add_subplot(111)
    #     i = self.axes.imshow(data)
    #     self.cb = self.fig.colorbar(i)
    #     self.cb.set_label('W/cm^2')
    #     self.axes.axhspan(511-BACKGROUNDLINES, 511, facecolor=(1.0,0,0), alpha=0.4)
    #     self.axes.text( 255, 511-BACKGROUNDLINES/2.0, 'BACKGROUND TAKEN FROM HERE', \
    #                         horizontalalignment='center', \
    #                         verticalalignment='center', \
    #                         color='red')        
    #     self.axes.set_xlim(0,511)
    #     self.axes.set_ylim(0,511)        
    #     self.axes.axhline( self.crosshairs_y )
    #     self.axes.axvline( self.crosshairs_x )
    #     self.draw()

    # def onclick(self, event):
    #     print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % \
    #         (event.button, event.x, event.y, event.xdata, event.ydata)
    #     self.crosshairs_x = event.xdata 
    #     self.crosshairs_y = event.ydata
    #     self.update_figure(self.data)


class ApplicationWindow(QtGui.QMainWindow):
    rawdata = None
    def __init__(self):

        self.spefiles = []
        self.m = None
        self.pwd = os.path.dirname(os.path.abspath(__file__))
        self.optical_element = 'Polarizer'

        QtGui.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setGeometry(QtCore.QRect(700,50,700,700))
        self.setWindowTitle("Artificial molecule analysis")
        

        self.main_widget = QtGui.QWidget(self)

        l = QtGui.QGridLayout(self.main_widget)

        self.getDirButton = QtGui.QPushButton('Get path to data', self.main_widget)
        self.getDirButton.clicked.connect(self.getDir)

        self.sc = MyStaticMplCanvas(self.main_widget, width=2, height=2, dpi=100)

        self.globalPhaseOffsetEdit = QtGui.QLineEdit(self.main_widget)
        self.globalPhaseOffsetEdit.setText( '9.0' )
#        self.globalPhaseOffsetEdit.setFixedWidth( 50 )
#        poweredit.setValidator(QtGui.QDoubleValidator)
        self.globalPhaseOffsetEdit.editingFinished.connect( self.update_global_phase )
        self.update_global_phase()

        self.SNREdit = QtGui.QLineEdit(self.main_widget)
        self.SNREdit.setText( '10' )

        self.fileChanger = QtGui.QComboBox(self.main_widget)
        self.fileChanger.activated.connect( self.update_file_change )

        self.lambdaOver2Checkbox = QtGui.QCheckBox('lambda/2 plate?',self.main_widget)
        self.lambdaOver2Checkbox.stateChanged.connect(self.optical_element_change)
        self.lambdaOver2Checkbox.setCheckState(False)

        self.setBGRegionButton     = QtGui.QPushButton('Set BG region', self.main_widget)
        self.setBGRegionButton.clicked.connect(self.set_bg_region)
        self.setSignalRegionButton = QtGui.QPushButton('Set signal region', self.main_widget)
        self.setSignalRegionButton.clicked.connect(self.set_signal_region)

        self.runButton = QtGui.QPushButton('run', self.main_widget)
        self.runButton.clicked.connect(self.runAnalysis)
        self.runAllButton = QtGui.QPushButton('run all', self.main_widget)
        self.runAllButton.clicked.connect(self.runAllAnalysis)

        self.pp1 = MyStaticMplCanvas(self.main_widget, width=1, height=1, dpi=100)
        self.pp1.clear()
        self.pp2 = MyStaticMplCanvas(self.main_widget, width=1, height=1, dpi=100)
        self.pp2.clear()
        self.pp3 = MyStaticMplCanvas(self.main_widget, width=1, height=1, dpi=100)
        self.pp3.clear()
        self.pp4 = MyStaticMplCanvas(self.main_widget, width=1, height=1, dpi=100)
        self.pp4.clear()


#        self.pick_reporter  = QtGui.QLabel('0 W/cm^{-2} \t x=0 \t y=0' )

        l.addWidget( QtGui.QLabel('global phase offset in deg:'), 1, 0 )
        l.addWidget( self.globalPhaseOffsetEdit, 1,1 )
        l.addWidget( self.lambdaOver2Checkbox, 1,2)
        l.addWidget( QtGui.QLabel('SNR:'), 1, 3 )
        l.addWidget( self.SNREdit, 1,4 )
        l.addWidget( self.getDirButton,2,0,1,2)
        l.addWidget( self.fileChanger,3,0,1,2)
        l.addWidget( self.setBGRegionButton,2,2 )
        l.addWidget( self.setSignalRegionButton,3,2 )
        l.addWidget( self.sc, 4,0,4,4 )
        l.addWidget( self.pp1, 4,4,1,1 )
        l.addWidget( self.pp2, 5,4,1,1 )
        l.addWidget( self.pp3, 6,4,1,1 )
        l.addWidget( self.pp4, 7,4,1,1 )
        l.addWidget( self.runButton, 8,0 )
        l.addWidget( self.runAllButton, 8,1 )

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

    # def pick_report(self, event):
    #     xi = np.round(event.xdata)
    #     yi = np.round(event.ydata)
    #     self.pick_reporter.setText( "%5.3e W/cm^{-2} \t x=%d \t y=%d" % (self.data[yi,xi], xi, yi) )

    def optical_element_change(self):
        if self.lambdaOver2Checkbox.isChecked():
            self.optical_element = 'l/2'
        else:
            self.optical_element = 'polarizer'

    def runAnalysis(self,fileindex=None):
        if fileindex==None or fileindex==False:
            fileindex = self.fileChanger.currentIndex()

        callstring  = 'python '+"\""+os.path.normpath(self.pwd+'/am_analyse.py')+"\""+' '
        callstring += "\""+os.path.normpath(self.data_directory+'/'+self.spefiles[fileindex])+"\""+' '
        callstring += "\""+os.path.normpath(self.data_directory+'/'+self.motorfiles[fileindex])+"\""+' '
        callstring += str(self.global_phase)+' '
        for b in self.bg_region_coords:
            callstring += str(b)+' '
        for s in self.signal_region_coords:
            callstring += str(s)+' '
        callstring += str(self.SNREdit.text())
        os.system(callstring)
        
        portrait = np.load('spotmatrix.npy')
        if fileindex<4:
            getattr(self, 'pp'+str(fileindex+1)).show_portrait(portrait)

        # self.m.define_background_spot( self.bg_region_coords )
        # self.m.define_spot( self.signal_region_coords )
        # self.m.chew_a_bit()

    def runAllAnalysis(self):
        for i in range(len(self.spefiles)):
#            self.load_and_display_spe_file(fileindex=i)
            self.runAnalysis(fileindex=i)

    def set_bg_region(self):
        coords = np.round( np.array( [self.sc.anno.x0, self.sc.anno.y0, \
                                          self.sc.anno.x1, self.sc.anno.y1] ) ).astype(np.int)
        if coords[0]>coords[2]:
            coords[0],coords[2]=coords[2],coords[0]
        if coords[1]>coords[3]:
            coords[1],coords[3]=coords[3],coords[1]

        self.bg_region_coords = coords
        self.sc.bg_rect.set_xy((coords[0],coords[1]))
        self.sc.bg_rect.set_width(coords[2]-coords[0])
        self.sc.bg_rect.set_height(coords[3]-coords[1])
        # reset selection rectangle
        self.sc.anno.rect.set_xy((0,0))
        self.sc.anno.rect.set_width(0)
        self.sc.anno.rect.set_height(0)
        # update canvas
        self.sc.figure.canvas.draw()                       

    def set_signal_region(self):
        coords = np.round( np.array( [self.sc.anno.x0, self.sc.anno.y0, \
                                          self.sc.anno.x1, self.sc.anno.y1] ) ).astype(np.int)
        if coords[0]>coords[2]:
            coords[0],coords[2]=coords[2],coords[0]
        if coords[1]>coords[3]:
            coords[1],coords[3]=coords[3],coords[1]
        self.signal_region_coords = coords
        self.sc.signal_rect.set_xy((coords[0],coords[1]))
        self.sc.signal_rect.set_width(coords[2]-coords[0])
        self.sc.signal_rect.set_height(coords[3]-coords[1])
        # reset selection rectangle
        self.sc.anno.rect.set_xy((0,0))
        self.sc.anno.rect.set_width(0)
        self.sc.anno.rect.set_height(0)
        # update canvas
        self.sc.figure.canvas.draw()                       

    def update_global_phase(self):
        self.global_phase = self.globalPhaseOffsetEdit.text().toDouble()[0]
        print "Global phase offset set to %f [deg]" % self.globalPhaseOffsetEdit.text().toDouble()[0]
#            self.sc.update_figure(self.data)

    def update_file_change(self):
        print self.fileChanger.currentIndex()
        self.load_and_display_spe_file(fileindex=self.fileChanger.currentIndex())

    def getDir(self):
        fname = QtGui.QFileDialog.getExistingDirectory(self, 'Show me where the data is', \
                directory='' )
        print "Selected data directory ",fname

        self.data_directory = str(fname)

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
        
        self.fileChanger.clear()
        for s in self.spefiles:
            self.fileChanger.addItem(s)
        self.fileChanger.setCurrentIndex(0)

        if not len(self.spefiles)==0:
            self.load_and_display_spe_file(0)

    def load_and_display_spe_file(self,fileindex=0):
        del(self.m)
        # from guppy import hpy; h=hpy()
        # w=h.heap()
        # print w
        print "loading file %s ... " % (self.spefiles[0]),
        sys.stdout.flush()
        self.m = Movie( self.data_directory+"/"+self.spefiles[fileindex], \
                            self.data_directory+"/"+self.motorfiles[fileindex], \
                            phase_offset_excitation=self.global_phase*np.pi/180.0, \
                            use_new_fitter=True, \
                            which_setup='cool new setup', \
                            excitation_optical_element=self.optical_element)
        self.sc.axes.imshow( self.m.camera_data.rawdata[0,:,:], zorder=1, cmap=cmap.gray )
        self.sc.draw()
        print "done"


    def fileQuit(self):
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()


qApp = QtGui.QApplication(sys.argv)

aw = ApplicationWindow()
#aw.setWindowTitle("")
aw.show()
sys.exit(qApp.exec_())
#qApp.exec_()



