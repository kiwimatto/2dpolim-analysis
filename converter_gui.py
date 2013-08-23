from PyQt4 import QtCore,QtGui
import sys, os, time
import util2dpolim.gui.converter_gui_layout as converter_gui_layout
import numpy as np
from util2dpolim.movie import Movie
from util2dpolim.misc import *
from util2dpolim.motors import *
from util2dpolim.cameradata import CameraData
import matplotlib.cm as cm
from matplotlib.patches import Rectangle

class converter_gui_logic(QtGui.QMainWindow,converter_gui_layout.Ui_MainWindow):

    def __init__(self,parent=None,app=None):
        """
            Initialization of the class. Call the __init__ for the super classes
        """
        super(converter_gui_logic,self).__init__(parent)
        self.setupUi(self)
        self.connectActions()

        self.pwd = os.path.dirname(os.path.abspath(__file__))
        self.app = app

        self.filedir = self.pwd   # to be changed to the dir containing the last selected file

        self.whichSetupComboBoxChanged()

    def main(self):
        self.show()

    def connectActions(self):
        """Connect the user interface controls to the logic """
        self.whichSetupComboBox.activated.connect( self.whichSetupComboBoxChanged )
        self.pickExMotorFilePushButton.clicked.connect( self.pickExMotorFile )
        self.pickEmMotorFilePushButton.clicked.connect( self.pickEmMotorFile )
        self.pickSPEFilePushButton.clicked.connect( self.pickSPEFile )
#        self.pickOutputFilePushButton.clicked.connect( self.pickOutputFile )
        self.convertPushButton.clicked.connect( self.convert )
        # self.phaseOffsetLineEdit.editingFinished.connect( self.setPhaseOffset )

    def whichSetupComboBoxChanged(self):
        isetup = self.whichSetupComboBox.currentIndex()
        print isetup
        if isetup==0:   # old setup
            self.pickExMotorFilePushButton.setText('Pick excitation motor file')
            self.pickEmMotorFilePushButton.setText('Pick emission motor file')
            self.pickEmMotorFilePushButton.setEnabled(True)
            self.pickSPEFilePushButton.setEnabled(True)
        elif isetup==1:   # new setup old style
            self.pickExMotorFilePushButton.setText('Pick motors file')
            self.pickEmMotorFilePushButton.setEnabled(False)
            self.pickSPEFilePushButton.setEnabled(True)
        elif isetup==2:   # new setup but missing header
            self.pickExMotorFilePushButton.setText('Pick motors file')
            self.pickEmMotorFilePushButton.setEnabled(False)
            self.pickSPEFilePushButton.setEnabled(False)

    def pickExMotorFile(self):
        path = QtGui.QFileDialog.getOpenFileName(self, 'Pick motor file', \
                                                      directory=self.filedir, filter='motor files (*.txt)')
        path = str(path)
        if not path=='':
            self.exMotorFileLabel.setText(path)
            self.exMotorFileName = path
            self.filedir  = os.path.split( os.path.normpath(path) )[0]
            self.filename = os.path.split( os.path.normpath(path) )[1]
            # trigger autocompletion check
            print self.filename
            if self.filename.startswith('MS-'):
                basename = self.filename[3:-4]
                testpath = os.path.normpath( self.filedir+'/'+basename )
                if os.path.isfile( testpath+'.SPE' ):
                    self.setSPEFile( path=testpath+'.SPE' )
                elif os.path.isfile( testpath+'.spe' ):
                    self.setSPEFile( path=testpath+'.spe' )

    def pickEmMotorFile(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Pick motor file (emission)', \
                directory=self.filedir, filter='motor files (*.txt)')
        fname = str(fname)
        if not fname=='':
            self.emMotorFileLabel.setText(fname)
            self.EmMotorFileName = fname
            self.filedir = os.path.split( os.path.normpath(fname) )[0]

    def pickSPEFile(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Pick SPE file (for a time axis)', \
                                                      directory=self.filedir, filter='SPE files (*.spe *.SPE)')
        fname = str(fname)
        if not fname=='':
            self.SPEFileLabel.setText(fname)
            self.SPEFileName = fname
            self.filedir = os.path.split( os.path.normpath(fname) )[0]

    def setSPEFile(self,path):
        self.SPEFileLabel.setText( path )
        self.SPEFileName = path
        self.filedir = os.path.split( os.path.normpath(path) )[0]
    
    # def pickOutputFile(self):
    #     fname = QtGui.QFileDialog.getSaveFileName(self, 'Pick output motor file', \
    #             directory='', filter='output motor files (*.txt)')
    #     self.OutputFileLabel.setText(fname)
    #     self.OutputFileName = fname

    def compile_header(self):
        header = []
        header.append( 'optical element in excitation: %s' % self.opticalElementInExcitationComboBox.currentText() )
        header.append( 'rotation mode in excitation: stepwise' )
        header.append( 'rotation mode in emission: stepwise' )
        header.append( 'phase offset in deg: %f' % self.phaseOffsetLineEdit.text().toFloat()[0] )
        header.append( 'full excitation power in mW: %f' % self.fullExPowerLineEdit.text().toFloat()[0] )
        header.append( 'optical density of filter: %f' % self.ODofFilterLineEdit.text().toFloat()[0] )
        header.append( 'excitation power after filter: %f' % self.exPowerAfterFilterLineEdit.text().toFloat()[0] )
        header.append( 'real gain: %f' % self.realGainLineEdit.text().toFloat()[0] )
        header.append( 'camera: %s' % self.whichCameraComboBox.currentText() )
        header.append( 'laser wavelength in nm: %f' % self.wavelengthLineEdit.text().toFloat()[0] )
        header.append( 'objective: %s' % self.whichObjectiveComboBox.currentText() )
        header.append( 'NA: %f' % self.NALineEdit.text().toFloat()[0] )
        header.append( 'isotropic emission modulation depth: %f' % self.modulationDepthLineEdit.text().toFloat()[0] )
        header.append( 'isotropic emission phase: %f' % self.phaseLineEdit.text().toFloat()[0] )
        header.append( 'user notes:' )
        header.append( str(self.plainTextEdit.toPlainText()) )
        header.append( 'END-OF-HEADER' )
        self.header = header
        
    def convert(self):
        self.compile_header()

        po = self.phaseOffsetLineEdit.text().toFloat()[0]
        if np.isnan(po):    # if phase offset was not set, then this must be a AM measurement, therefore:
            po = 0
        oe = str(self.opticalElementInExcitationComboBox.currentText())

        isetup = self.whichSetupComboBox.currentIndex()
        if isetup==0:    # old setup
            if self.rotationDirectionComboBox.currentIndex==1:
                rotdir=-1
            else:
                rotdir=1

            self.excitation_motor = ExcitationMotor( self.exMotorFileName, \
                                                         phase_offset_in_deg=po, \
                                                         rotation_direction=rotdir, \
                                                         optical_element=oe )
            self.emission_motor   = EmissionMotor( self.emMotorFileName )

            self.camera_data    = CameraData( self.SPEFileName, compute_frame_average=False )
            self.timeaxis       = self.camera_data.timestamps

            self.exangles = np.array( [self.excitation_motor.angle(t,exposuretime=self.camera_data.exposuretime) for t in self.timeaxis] )
            self.emangles = np.array( [self.emission_motor.angle(t,exposuretime=self.camera_data.exposuretime) for t in self.timeaxis] )

            # In principle we got what we need now: header, and which angles belong to which frame;
            # however, the emangles by default respect the shutter, that is, are -1 if their integration
            # time overlaps with the shutter opening/closing. We need to screen for this: 'validdata'.

            # array of valid indices:
            valid = np.nonzero( self.emangles!=-1 )[0]
            self.framenumbers = range(self.camera_data.datasize[0])[valid]
            self.exangles = self.exangles[valid]
            self.emangles = self.emangles[valid]

            # As a consequence, the frame numbers may miss a beat or two, ie [1,2,3,6,7,...]
            # and, in the later analysis, only those frame must be loaded from the SPE file.
            # An alternative would be to re-save the SPE file containing only the valid frames,
            # but I neither have an SPE writer nor the time to implement one.


        elif isetup==1:    # new setup, old style
            self.camera_data    = CameraData( self.SPEFileName, compute_frame_average=False )
            self.timeaxis       = self.camera_data.timestamps
            self.motors = NewSetupMotors( self.exMotorFileName, \
                                              phase_offset_in_deg=po, \
                                              optical_element=oe )
            self.motors.compute_angles( self.camera_data.timeaxis, self.camera_data.exposuretime )
            self.exangles = np.array( self.motors.excitation_angles )
            self.emangles = np.array( self.motors.emission_angles )
            # we have the same issue here as above, the angles may be set to -1 because of shutter
            # overlap, so:
            valid = np.nonzero( self.emangles!=-1 )[0]
            self.framenumbers = np.arange(self.camera_data.datasize[0])[valid]
            self.exangles = self.exangles[valid]
            self.emangles = self.emangles[valid]

        elif isetup==2:    # cool new setup, but without header
            self.motors = BothMotors( self.exMotorFileName, \
                                              phase_offset_in_deg=po )
            self.framenumbers = self.motors.framenumbers
            self.exangles = np.array( self.motors.excitation_angles )
            self.emangles = np.array( self.motors.emission_angles )


        # first move the original exMotorFile
        os.rename(self.exMotorFileName,self.exMotorFileName+'.backup')

        # now open a new file of the same name
        fid = open(self.exMotorFileName,'w')
        # write header
        for h in self.header: 
            fid.write(h+'\n')
            print h
        # write table head
        fid.write('Frame\tExcitation Angle\tEmission Angle\n')
        print 'Frame\tExcitation Angle\tEmission Angle'
        # write table
        for i,f in enumerate(self.framenumbers):
            fid.write( '%d\t%f\t%f\n' % (f,self.exangles[i],self.emangles[i]) )
            print '%d\t%f\t%f' % (f,self.exangles[i],self.emangles[i])
        # and close the file
        fid.close()
        

if __name__=='__main__':
    app = QtGui.QApplication(sys.argv)
    instance = converter_gui_logic(app=app)
    instance.main()
    sys.exit(app.exec_())
