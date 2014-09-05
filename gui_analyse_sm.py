# coding=utf-8

from PyQt4 import QtCore,QtGui
import sys, traceback, os, time
import numpy as np
import layout_gui_analyse_sm
from util2dpolim.spefiles import MyPrincetonSPEFile
from util2dpolim.movie import Movie
from util2dpolim.misc import *

class gui_analyse_sm(QtGui.QMainWindow, layout_gui_analyse_sm.Ui_MainWindow):

    def __init__(self, parent=None, app=None):
        super(gui_analyse_sm, self).__init__(parent)
        self.setupUi(self)
        self.connectActions()

        self.pwd = os.path.dirname(os.path.abspath(__file__))
        self.app = app
        self.data_directory = self.pwd   # to be changed to the dir containing the last selected file
        self.spefile = ''

        self.gotmotor = False
        self.gotblank = False
        self.gotspotcoords = False
        self.run_check()
        self.spotTypeChanged()
        
    def main(self):
        self.show()

    def connectActions(self):
        self.selectDirectoryPushButton.clicked.connect( self.selectDirectory )
        self.selectSPEComboBox.currentIndexChanged.connect( self.selectSPE )
        self.runAnalysisPushButton.clicked.connect( self.run_analysis )
        self.spotTypeComboBox.currentIndexChanged.connect( self.spotTypeChanged )

    def selectDirectory(self):
        # try to start from previous dir
        if os.path.isfile( self.pwd+os.path.sep+'lastdir.txt' ):
            f = open( self.pwd+os.path.sep+'lastdir.txt', 'r')
            lastdir = os.path.normpath( f.readlines()[0] )
            f.close()
            if not len(lastdir)==0 and os.path.isdir(lastdir):
                self.data_directory = lastdir

        dirname = QtGui.QFileDialog.getExistingDirectory(self, 'select data directory', \
                                                          directory=self.data_directory )
        dirname = str(dirname)
        self.data_directory = os.path.normpath(dirname)

        if not dirname=='':
            self.statusbar.showMessage("Looking for SPE data...")
            # get all filenames
            self.spefiles = []
            for f in os.listdir( self.data_directory ):
                # get all that end in spe
                if f.endswith(".spe") or f.endswith(".SPE"):
                    if not f.startswith("blank-"):
                        self.spefiles.append(f)
            self.spefiles = list(set(self.spefiles))
            self.spefiles.sort()

            if len(self.spefiles)>0:
                # register in combobox
                for s in self.spefiles:
                    self.selectSPEComboBox.addItem(s)
                self.selectSPEComboBox.setCurrentIndex(0)
                self.statusbar.showMessage("%d spe files here." % len(self.spefiles))
            else:
                self.statusbar.showMessage("could not find any spe files here.")

            self.selectDirectoryPushButton.setText( self.data_directory )

            # save in lastdir
            f = open( self.pwd+os.path.sep+'lastdir.txt', 'w')
            f.writelines( [self.data_directory] )
            f.close()
            

    def selectSPE(self):
        self.spefile = str(self.selectSPEComboBox.currentText())
        # check if we have a motor file
        if os.path.isfile( self.data_directory+os.path.sep+'MS-'+self.spefile[:-4]+'.txt' ):
            self.gotmotor = True
            self.motorFileLabel.setText(" motor file found ")
            self.motorFileLabel.setStyleSheet("QLabel { color : black; font-weight: normal }")
        else:
            self.gotmotor = False
            self.motorFileLabel.setText(" -- NO MOTOR FILE ! -- ")
            self.motorFileLabel.setStyleSheet("QLabel { color : red; font-weight: bold }")
        # check if we have a blank file
        if os.path.isfile( self.data_directory+os.path.sep+'blank-'+self.spefile ):
            self.gotblank = True
            self.blankFileLabel.setText(" blank file found ")
            self.blankFileLabel.setStyleSheet("QLabel { color : black; font-weight: normal }")
        else:
            self.gotblank = False
            self.blankFileLabel.setText(" -- NO BLANK FILE -- ")
            self.blankFileLabel.setStyleSheet("QLabel { color : orange; font-weight: bold }")
        # check if we have spot coordinates file
        if os.path.isfile( self.data_directory+os.path.sep+'spotcoordinates_'+self.spefile[:-4]+'.txt' ):
            self.gotspotcoords = True
            self.spotCoordsFileLabel.setText(" spot coordinates file found ")
            self.spotCoordsFileLabel.setStyleSheet("QLabel { color : black; font-weight: normal }")
        else:
            self.gotspotcoords = False
            self.spotCoordsFileLabel.setText(" -- NO SPOT COORDINATES FILE -- ")
            self.spotCoordsFileLabel.setStyleSheet("QLabel { color : red; font-weight: bold }")
        
        # check if we can do the analysis with what we got
        self.run_check()

    def spotTypeChanged(self):
        if str(self.spotTypeComboBox.currentText())=='circle':
            self.spotSizeLabel.setText('radius (pixel)')
        elif str(self.spotTypeComboBox.currentText())=='square':
            self.spotSizeLabel.setText('edge length (pixel)')
        else:
            raise dough

    def run_check(self):
        if self.gotmotor and self.gotspotcoords:
            self.SNRSpinBox.setEnabled(True)
            self.VFRSpinBox.setEnabled(True)
            self.spotTypeComboBox.setEnabled(True)
            self.spotSizeSpinBox.setEnabled(True)
            self.runAnalysisPushButton.setEnabled(True)
        else:
            self.SNRSpinBox.setEnabled(False)
            self.VFRSpinBox.setEnabled(False)
            self.spotTypeComboBox.setEnabled(False)
            self.spotSizeSpinBox.setEnabled(False)
            self.runAnalysisPushButton.setEnabled(False)

    def run_analysis(self):
        self.statusbar.showMessage("working...")
        self.runAnalysisPushButton.setText("...wait...")
        QtGui.QApplication.processEvents()
        try:
            # definitions
            prefix   = self.data_directory
            basename = self.spefile[:-4]
            SNR = self.SNRSpinBox.value()
            VFR = self.VFRSpinBox.value()
            spotradius = self.spotSizeSpinBox.value()

            # basics
            m = Movie( prefix, basename )
            m.find_portraits( frameoffset=1 )
            m.find_lines()

            # blank
            if self.gotblank:
                boolimage = np.ones( (m.sample_data.rawdata.shape[1], \
                                          m.sample_data.rawdata.shape[2]), \
                                         dtype=np.bool )*True
                m.fit_blank_image( boolimage, verbosity=0 )

            # spot coordinates
            import_spot_positions( m, basename, spotradius, 'circle', use_borderbg=True )

            # the corrections
            m.correct_excitation_intensities()
            m.correct_emission_intensities()

            # spot validity
            m.are_spots_valid( SNR=SNR, validframesratio=VFR )
            if len(m.validspots)==0:
                raise Exception("No valid spots! (reduce SNR, VFR ?)")

            # cosine fits
            m.fit_all_portraits_spot_parallel_selective( range(len(m.validspots)) )

            # mod depths
            m.find_modulation_depths_and_phases_selective( range(len(m.validspots)) )

            # ET ruler
            for s in m.validspots:
                s.values_for_ETruler( newdatalength=1024 )

            # ET model
            m.ETmodel_selective( range(len(m.validspots)) )

            # hdf5 file output
            remove_hdf5_files( m )    
            save_spot_hdf5( m )
            save_hdf5( m )

        except Exception, e:
            self.statusbar.showMessage("Error")
            QtGui.QMessageBox.about( self, "Error", str(e) )
            traceback.print_exc()
        self.statusbar.showMessage("All done")
        self.runAnalysisPushButton.setText("run analysis")

if __name__=='__main__':
    app = QtGui.QApplication(sys.argv)
    instance = gui_analyse_sm(app=app)
    instance.main()
    sys.exit(app.exec_())
