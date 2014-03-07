# coding=utf-8

from PyQt4 import QtCore,QtGui
import sys, traceback, os, time
import numpy as np
import layout_sm_human_validator
from util2dpolim.spot import Spot

class sm_human_validator(QtGui.QMainWindow, layout_sm_human_validator.Ui_MainWindow):

    def __init__(self, parent=None, app=None):
        super(sm_human_validator, self).__init__(parent)
        self.setupUi(self)
        self.connectActions()

        self.pwd = os.path.dirname(os.path.abspath(__file__))
        self.app = app
        self.data_directory = self.pwd   # to be changed to the dir containing the last selected file

        self.grab_directory()


    def main(self):
        self.show()

    def connectActions(self):
        pass

    def grab_directory(self):
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
            self.statusbar.showMessage("Looking for spot output hdf5 files...")
            # get all filenames
            self.hdf5files = []
            for f in os.listdir( self.data_directory ):
                # get all that end in spe
                if f.endswith("spot_output.hdf5"):
                    self.hdf5files.append(f)
            self.hdf5files = list(set(self.hdf5files))
            self.hdf5files.sort()

            # save in lastdir
            f = open( self.pwd+os.path.sep+'lastdir.txt', 'w')
            f.writelines( [self.data_directory] )
            f.close()
            

    def validate_all_files(self):        
        for hdf5file in self.hdf5files:
            self.validate_single_file( hdf5file )

    def validate_single_file( hdf5file ):
        for spotname in hdf5file.keys():
            intensity = hdf5file[ spotname+'/intensity' ]
            self.axes_ints.plot( intensity )


if __name__=='__main__':
    app = QtGui.QApplication(sys.argv)
    instance = sm_human_validator(app=app)
    instance.main()
    sys.exit(app.exec_())
