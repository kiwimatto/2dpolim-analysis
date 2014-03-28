# coding=utf-8

from PyQt4 import QtCore,QtGui
import sys, traceback, os, time
import numpy as np
import layout_sm_human_validator
from util2dpolim.spot import Spot
from util2dpolim.fitting import CosineFitter_new
import matplotlib.cm as cm
import h5py

class sm_human_validator(QtGui.QMainWindow, layout_sm_human_validator.Ui_MainWindow):

    def __init__(self, parent=None, app=None):
        super(sm_human_validator, self).__init__(parent)
        self.setupUi(self)
        self.connectActions()

        self.pwd = os.path.dirname(os.path.abspath(__file__))
        self.app = app
        self.data_directory = self.pwd   # to be changed to the dir containing the last selected file

        self.cfi = -1    # current file index
        self.csi = -1    # current spot index

        self.grab_directory()

        self.curmsg = None

        self.next()

    def main(self):
        self.show()

    def connectActions(self):
        self.yarpPushButton.clicked.connect( self.next )
        self.narpPushButton.clicked.connect( self.narpie )
        self.meanintPushButton.clicked.connect( self.draw_mean_int_and_spot_coverage )
        self.coveragePushButton.clicked.connect( self.draw_mean_int_and_spot_coverage )
        self.imlostPushButton.clicked.connect( self.imlost )
        
    def imlost(self):
        if self.imlostPushButton.isChecked():
            self.intPlotWidget.setVisible(False)
            self.portraitPlotWidget.setVisible(False)
            self.projectionsPlotWidget.setVisible(False)
            self.ETmodelPlotWidget.setVisible(False)
            self.spotcoveragePlotWidget.setVisible(False)
            self.intPlotToolsWidget.setVisible(False)
            self.PortraitPlotToolsWidget.setVisible(False)
            self.ProjectionsPlotToolsWidget.setVisible(False)
            self.spotcoveragePlotToolsWidget.setVisible(False)
            self.yarpPushButton.setVisible(False)
            self.narpPushButton.setVisible(False)
            self.centralwidget.setStyleSheet("background-image: url(bg_anno.png);")
            self.imlostPushButton.setStyleSheet("")
        else:
            self.intPlotWidget.setVisible(True)
            self.portraitPlotWidget.setVisible(True)
            self.projectionsPlotWidget.setVisible(True)
            self.ETmodelPlotWidget.setVisible(True)
            self.spotcoveragePlotWidget.setVisible(True)
            self.intPlotToolsWidget.setVisible(True)
            self.PortraitPlotToolsWidget.setVisible(True)
            self.ProjectionsPlotToolsWidget.setVisible(True)
            self.spotcoveragePlotToolsWidget.setVisible(True)
            self.yarpPushButton.setVisible(True)
            self.narpPushButton.setVisible(True)
            self.centralwidget.setStyleSheet("")

    def narpie(self):
        self.curmsg = None    
        self.next()

    def next(self):
        # if there's no file to do anything with, bomb out here
        if self.cfi == -1:
            print 'no files found'

        # save current message (if any)
        if not os.path.exists(self.data_directory + os.path.sep + 'smhv_results.txt'):
            f = open(self.data_directory + os.path.sep + 'smhv_results.txt','w')
            f.writelines(['Int\tM_ex\tM_em\tphase_ex\tphase_em\tET_M\tET_phase\tET_gr\tET_et\n'])
        else:
            f = open(self.data_directory + os.path.sep + 'smhv_results.txt','a')
        if not self.curmsg==None:
            f.writelines([self.curmsg+'\n'])
        f.close()

        # next
        self.intPlotWidget.clear()
        self.portraitPlotWidget.clear()
        self.projectionsPlotWidget.clear()
        self.ETmodelPlotWidget.clear()

        # if we have not run out of spots in this file yet
        if self.csi < self.NspotsList[self.cfi]-1:
            # ask for the next spot
            print 'next spot'
            self.csi += 1        
        else:
            print 'new file?...',
            if self.cfi < len(self.hdf5files)-1:
                print 'yes'
                # ask for next file
                self.cfi += 1
                self.csi = 0
            else:
                print 'nope'
                return

        self.show_spot()


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

        # now make a list of the number of spots in each hdf5
        NspotsList = []
        for h5f in self.hdf5files:
            fid = h5py.File( self.data_directory + os.path.sep + h5f, 'r' )
            Nspots = 0
            for k in fid.keys():
                if k.startswith('spot_'):
                    Nspots += 1
            NspotsList.append( Nspots )
            fid.close()
        self.NspotsList = NspotsList

        # re-build list of hdf5 files to include only those with spots
        flist = []
        for i,N in enumerate(NspotsList):
            if N>0:
                flist.append( self.hdf5files[i] )
        self.hdf5files = flist        

        # if files were found, point the current file index to the first one
        if len(self.hdf5files)>0:
            self.cfi = 0

        print "Files found:"
        print self.hdf5files


    def show_spot(self):
        spotname  = '/spot_00%04d' % self.csi
        filename = self.data_directory + os.path.sep + self.hdf5files[self.cfi]
        labelmsg  = self.hdf5files[self.cfi]
        labelmsg += " ---  spot %d of %d" % (self.csi+1,self.NspotsList[self.cfi])
        self.label.setText( labelmsg )

        print "Looking at %s in file %s" % (spotname, filename)
        f = h5py.File( filename, 'r' )
        
        pis = np.array( f['portrait_indices'] )
        lis = np.array( f['line_indices'] )
        exangles = np.array( f['exangles'] )
        intensity = np.array( f[spotname+'/intensity'] )

        Nportraits = len(lis)
        Nlinesperportrait = len(lis[0])

        # print pis
        # print lis
        # print len(lis)
        # print len(lis[0])
        # print exangles
        # print intensity
        # print '----------------------------'

        # collect intensity data into lines, where exangles are sorted (ascending)
        linedata = []
        for pi in range(len(lis)):               # for each portrait
            portint = intensity[pis[pi]:pis[pi+1]]
            for li in range(len(lis[pi])):       # for each line
                exa = exangles[ lis[pi][li] ]
                exasorti = np.argsort(exa)
                exasort  = exa[exasorti]
                sig = portint[ lis[pi][li] ]
                sigsort  = sig[exasorti]
                linedata.append( np.vstack([exasort,sigsort]).T )

        ymax = 1.2* np.max( np.array(linedata)[:,:,1] )
        ymin = 0 # np.min( np.array(linedata)[:,:,1] )

        ### intensity ###
        # subplots for each line
        self.intPlotWidget.fig.clf()
        intAxes = [self.intPlotWidget.fig.add_subplot(1,len(linedata),i+1) for i in range(len(linedata))]
        # plot line
        for i in range(len(linedata)):
            ld = linedata[i]
            intAxes[i].plot( ld[:,0], ld[:,1], 'b-s', alpha=.8, markerfacecolor='none' )
            intAxes[i].set_xlim( np.min(ld[:,0]), np.max(ld[:,0]) )
            intAxes[i].set_ylim( ymin, ymax )
            intAxes[i].set_yticklabels([])
            # if we're past the first portrait, then also plot the previous as a dotted line
            if i>=Nlinesperportrait:           
                ldold = linedata[ i - Nlinesperportrait ]  # the same line in the previous portrait
                intAxes[i].plot( ldold[:,0], ldold[:,1], 'b:', alpha=.8, markerfacecolor='none' )
                # if the x-axes of old and current line are the same (they should be), then
                # we can also plot the fill-between area
                if np.all(ldold[:,0]==ld[:,0]):                        
                    intAxes[i].fill_between( ldold[:,0], ldold[:,1], ld[:,1], facecolor='blue', alpha=.05 )

        ### line fits ###
        lfps = np.array( f[spotname+'/fits/linefitparams'] )
        lfpsr = lfps.copy()
        lfpsr = lfpsr.reshape( (lfpsr.shape[0]*lfpsr.shape[1], lfpsr.shape[2]) )

        # this time we don't collect first, but plot things straight-away
        mycos = lambda a, ph, I, M: I*( 1+M*( np.cos(2*(a-ph)) ))
        for i in range(len(linedata)):
            ld = linedata[i]
            fineexa = np.linspace( np.min(ld[:,0]), np.max(ld[:,0]), 40 )
            fit       = mycos( fineexa, lfpsr[i,0],lfpsr[i,1],lfpsr[i,2] )
            fitcoarse = mycos( ld[:,0], lfpsr[i,0],lfpsr[i,1],lfpsr[i,2] )
            
            intAxes[i].plot( fineexa, fit, 'r-' )
            intAxes[i].plot( ld[:,0], fitcoarse, 'rx' )
            # local phase
            intAxes[i].plot( [lfpsr[i,0],lfpsr[i,0]], [0,np.interp(lfpsr[i,0], fineexa, fit)], 'r--' )

            # previous portrait fits..
            if i>=Nlinesperportrait:
                j = i - Nlinesperportrait
                ldold = linedata[j]
                fineexa = np.linspace( np.min(ldold[:,0]), np.max(ldold[:,0]), 40 )
                fitprev = mycos( fineexa, lfpsr[j,0],lfpsr[j,1],lfpsr[j,2] )
#                fitcoarse = mycos( ld[:,0], lfpsr[j,0],lfpsr[j,1],lfpsr[j,2] )

                intAxes[i].plot( fineexa, fitprev, 'r:' )
                # local phase
                intAxes[i].plot( [lfpsr[j,0],lfpsr[j,0]], [0,np.interp(lfpsr[j,0], fineexa, fit)], 'r:' )
                # fill between
                if np.all(ldold[:,0]==ld[:,0]):
                    intAxes[i].fill_between( fineexa, fitprev, fit, facecolor='red', alpha=.05 )

        # intensity = np.array( f[spotname+'/intensity'] )
        # self.intPlotWidget.axes_ints.plot( range(intensity.size), intensity )
        # for pi in range(1,len(pis)-1):
        #     self.intPlotWidget.axes_ints.plot(range(pis[pi],pis[pi+1]), intensity[pis[pi-1]:pis[pi]],'b:')
        # self.intPlotWidget.axes_ints.set_xlim( -1, intensity.size )
        
        # ### draw portrait boundaries ###
        # self.intPlotWidget.axes_ints.axvline( pis[0], color='g', ls='-' )   # first frame of first portrait
        # for p in pis[1:-1]:
        #     self.intPlotWidget.axes_ints.axvline( p-1, color='g', ls='--' )    # last frame of this portrait
        #     self.intPlotWidget.axes_ints.axvline( p, color='g', ls='-' )       # first frame of next portrait
        # self.intPlotWidget.axes_ints.axvline( pis[-1]-1, color='g', ls='--' )  # last frame of last portrait

        # ### draw line boundaries ###
        # linestartframes = []
        # for pi in range(len(pis)-1):          # go through all portraits
        #     for li in range(len(lis[pi])):    # go through all lines in this portrait
        #         x = np.nonzero(lis[pi][li])[0][0] + pis[pi]
        #         linestartframes.append( x )
        #         self.intPlotWidget.axes_ints.axvline( x, ls=':' )

        # ### show line fits ###
        # lfps = np.array( f[spotname+'/fits/linefitparams'] )
        # exangles = np.array( f['exangles'] )
        # mycos = lambda a, ph, I, M: I*( 1+M*( np.cos(2*(a-ph)) ) )
        # for pi in range(len(pis)-1):          # go through all portraits
        #     for li in range(len(lis[pi])):
        #         exa = exangles[ pis[pi]:pis[pi+1] ][ lis[pi][li] ]
        #         x = range(linestartframes[li],linestartframes[li+1])+pis[pi]
        #         y = mycos(exa, lfps[pi][li][0],lfps[pi][li][1],lfps[pi][li][2])
        #         self.intPlotWidget.axes_ints.plot( x, y, 'r-' )

        ### show numerical portraits ###
        # first two numerical data portraits
        emangles = np.array( f['emangles'] ) 
#        ndp = np.zeros( (len(pis)-1, emangles.size, exangles.size) )
        ndp = np.zeros( (lis.shape[0], lis.shape[1], lis.shape[2]/lis.shape[1]) )
        for pi in range(len(pis)-1):
            for li in range(len(lis[pi])):
                ndp[pi,li,:] = intensity[ pis[pi]:pis[pi+1] ][ lis[pi][li] ]
        # first two portraits:
        for pi in range(len(pis)-1):
            self.portraitPlotWidget.axes[0,pi].imshow( ndp[pi], interpolation='nearest', \
                                                           extent=[np.min(exangles),np.max(exangles),\
                                                                       np.min(emangles),np.max(emangles)], \
                                                           origin='bottom' )
        # average numerical portrait
        self.portraitPlotWidget.axes[0,2].imshow( np.mean(ndp,axis=0), interpolation='nearest',\
                                                      extent=[np.min(exangles),np.max(exangles),\
                                                                  np.min(emangles),np.max(emangles)], \
                                                      origin='bottom' )

        ### build and display fitted portraits ###
        # need to first evaluate the line fits for finely spaced excitation angles
        fineexa = np.linspace(0, np.pi, 40, endpoint=False)
        fineema = np.linspace(0, np.pi, 40, endpoint=False)
        finelines = np.zeros( (lis.shape[1], fineexa.size) )
        fineportraits = np.zeros( (len(pis)-1, fineema.size, fineexa.size) )
        for pi in range(len(pis)-1):          # go through all portraits
            for li in range(len(lis[pi])):    # go through all lines in this portrait
                finelines[li,:] = mycos(fineexa, lfps[pi][li][0],lfps[pi][li][1],lfps[pi][li][2])
            # vertical fit for all those fine-ex-angles
            ema = np.unique(emangles)
            phase, I0, M, resi, fit, rawfitpars, mm = CosineFitter_new( ema, finelines, 91 ) 
            # and store in portrait array
            for j in range(len(phase)):
                fineportraits[pi,:,j] = mycos( fineema, phase[j], I0[j], M[j] )
        # first two fitted portraits
        for pi in range(len(pis)-1):
            self.portraitPlotWidget.axes[1,pi].imshow( fineportraits[pi], interpolation='nearest',\
                                                           origin='bottom')
        # mean fitted portrait
        self.portraitPlotWidget.axes[1,2].imshow( np.mean(fineportraits,axis=0), interpolation='nearest',\
                                                      origin='bottom')

        ### modulation depths and their fits ###
        exa = np.unique( np.array( f['excitation_angles_grid'] ) )
        ema = np.unique( np.array( f['emission_angles_grid'] ) )
        exafine = np.linspace( np.min(exa), np.max(exa), 40 )
        emafine = np.linspace( np.min(ema), np.max(ema), 40 )

        modexdata = np.array( f[spotname+'/fits/modulation_ex_data'] )
        modemdata = np.array( f[spotname+'/fits/modulation_em_data'] )
        # reconstruct fits
        I_ex     = np.array( f[spotname+'/contrasts/I_ex'] )
        M_ex     = np.array( f[spotname+'/contrasts/M_ex'] )
        phase_ex = np.array( f[spotname+'/contrasts/phase_ex'] )
        I_em     = np.array( f[spotname+'/contrasts/I_em'] )
        M_em     = np.array( f[spotname+'/contrasts/M_em'] )
        phase_em = np.array( f[spotname+'/contrasts/phase_em'] )
        modexfit = mycos( exafine, phase_ex, I_ex, M_ex )
        modemfit = mycos( emafine, phase_em, I_em, M_em )

        modexfit2 = np.array( f[spotname+'/fits/modulation_ex_fit'] )
        modemfit2 = np.array( f[spotname+'/fits/modulation_em_fit'] )
        self.projectionsPlotWidget.axes_mex.plot( exa, modexdata, 'b-s', markerfacecolor='none' )
        self.projectionsPlotWidget.axes_mex.plot( exa, modexfit2, 'rx' )
        self.projectionsPlotWidget.axes_mex.plot( exa, modexdata-modexfit2, 'r:' )
        self.projectionsPlotWidget.axes_mex.plot( exafine, modexfit, 'r-' )
        self.projectionsPlotWidget.axes_mex.plot( [phase_ex,phase_ex], \
                                                      [0, np.interp(phase_ex,exafine,modexfit)], 'r--' )
        self.projectionsPlotWidget.axes_mex.set_xlim( np.min(exa),np.max(exa) )
        self.projectionsPlotWidget.axes_mex.set_ylim( -.2*np.max(modexfit), 1.2*np.max(modexfit) )
        

        self.projectionsPlotWidget.axes_mem.plot( ema, modemdata, 'b-s', markerfacecolor='none' )
        self.projectionsPlotWidget.axes_mem.plot( ema, modemfit2, 'rx' )
        self.projectionsPlotWidget.axes_mem.plot( ema, modemdata-modemfit2, 'r:' )
        self.projectionsPlotWidget.axes_mem.plot( emafine, modemfit, 'r-' )
        self.projectionsPlotWidget.axes_mem.plot( [phase_em,phase_em], \
                                                      [0, np.interp(phase_em,emafine,modemfit)], 'r--' )
        self.projectionsPlotWidget.axes_mem.set_xlim( np.min(ema),np.max(ema) )
        self.projectionsPlotWidget.axes_mem.set_ylim( -.2*np.max(modemfit), 1.2*np.max(modemfit) )

        ### ET model ###

        md_ex = np.clip( M_ex, .000001, .999999 )

        # this is mostly copypasta from the fitting function
        md_fu = np.array( f[spotname+'/contrasts/ETmodel_md_fu'] )
        th_fu = np.array( f[spotname+'/contrasts/ETmodel_th_fu'] )
        gr = np.array( f[spotname+'/contrasts/ETmodel_gr'] )
        et = np.array( f[spotname+'/contrasts/ETmodel_et'] )

        # flattened meshgrid of the angles of the first portrait
#        EX, EM = exangles[pis[0]:pis[1]].flatten(), emangles[pis[0]:pis[1]].flatten()
        EX, EM = np.meshgrid( exafine, emafine )

        # calculate angle between off-axis dipoles in symmetric model
        if np.isnan(gr): gr=1.0
        alpha = 0.5 * np.arccos( .5*(((gr+2)*md_ex)-gr) )

        ph_ii_minus = phase_ex -alpha
        ph_ii_plus  = phase_ex +alpha

        # print EX
        # print EM
        # print np.cos( EX-ph_ii_minus )**2 * np.cos( EM-ph_ii_minus )**2
        # raise SystemExit
        Fnoet  =    np.cos( EX-ph_ii_minus )**2 * np.cos( EM-ph_ii_minus )**2
        Fnoet += gr*np.cos( EX-phase_ex )**2 * np.cos( EM-phase_ex )**2
        Fnoet +=    np.cos( EX-ph_ii_plus )**2 * np.cos( EM-ph_ii_plus )**2
        Fnoet /= (2.0+gr)
    
        Fet   = .25 * (1+md_ex*np.cos(2*(EX-phase_ex))) * (1+md_fu*np.cos(2*(EM-th_fu-phase_ex)))

        model = et*Fet + (1-et)*Fnoet

        self.portraitPlotWidget.axes[2,0].imshow( Fnoet.reshape( (emafine.size, exafine.size) ), \
                                                      interpolation='nearest', origin='bottom' )
        self.portraitPlotWidget.axes[2,1].imshow( Fet.reshape( (emafine.size, exafine.size) ), \
                                                      interpolation='nearest', origin='bottom' )
        self.portraitPlotWidget.axes[2,2].imshow( model.reshape( (emafine.size, exafine.size) ), \
                                                      interpolation='nearest', origin='bottom' )

        # now draw the model dipoles
        self.ETmodelPlotWidget.axes.plot( [phase_ex,phase_ex+np.pi], [gr,gr], 'k-', lw=2 )
        self.ETmodelPlotWidget.axes.plot( [ph_ii_minus,ph_ii_minus+np.pi], [1,1], 'r-', lw=2 )
        self.ETmodelPlotWidget.axes.plot( [ph_ii_plus,ph_ii_plus+np.pi], [1,1], 'r-', lw=2 )
        self.ETmodelPlotWidget.axes.set_yticklabels([])

        f.close()

        self.intPlotWidget.fig.canvas.draw()
        self.portraitPlotWidget.fig.canvas.draw()
        self.projectionsPlotWidget.fig.canvas.draw()
        self.ETmodelPlotWidget.fig.canvas.draw()

        # message for status bar
        msg  = 'Int=%.2f    ' % np.mean(intensity)
        msg += 'Mex=%.2f    ' % M_ex
        msg += 'Mem=%.2f    ' % M_em
        msg += 'phex=%.1f    ' % (phase_ex * 180.0/np.pi)
        msg += 'phem=%.1f    ' % (phase_em * 180.0/np.pi)
        msg += 'md_fu=%.1f    ' % md_fu
        msg += 'ph_fu=%.1f    ' % (th_fu * 180.0/np.pi)
        msg += 'gr=%.2f    ' % gr
        msg += 'et=%.2f    ' % et
        self.statusbar.showMessage( msg )

        # message for file
        msg  = '%.2f\t' % np.mean(intensity)
        msg += '%.2f\t' % M_ex
        msg += '%.2f\t' % M_em
        msg += '%.1f\t' % (phase_ex * 180.0/np.pi)
        msg += '%.1f\t' % (phase_em * 180.0/np.pi)
        msg += '%.1f\t' % md_fu
        msg += '%.1f\t' % (th_fu * 180.0/np.pi)
        msg += '%.2f\t' % gr
        msg += '%.2f\t' % et
        self.curmsg = msg

        ### spot coverage ###
        self.draw_mean_int_and_spot_coverage()


    def draw_mean_int_and_spot_coverage(self):
        self.spotcoveragePlotWidget.clear()
        spotname  = '/spot_00%04d' % self.csi
        filename = self.data_directory + os.path.sep + self.hdf5files[self.cfi]
        filename2 = self.data_directory + os.path.sep + self.hdf5files[self.cfi][:-17] + '_output.hdf5'
        f  = h5py.File( filename, 'r' )
        f2 = h5py.File( filename2, 'r' )
        meanint = np.array(f2['/images/original_mean_intensity_image'])
        spotcov = np.array(f2['/images/spot_coverage_image'])
        center  = np.array(f[spotname+'/setup/center'])
        margin = 15
        r0 = np.clip( center[0]-margin, 0, meanint.shape[0]-1 )
        r1 = np.clip( center[0]+margin, 0, meanint.shape[0]-1 )
        c0 = np.clip( center[1]-margin, 0, meanint.shape[1]-1 )
        c1 = np.clip( center[1]+margin, 0, meanint.shape[1]-1 )
        if self.meanintPushButton.isChecked():
            self.spotcoveragePlotWidget.axes.imshow( meanint, \
                                                         interpolation='nearest', cmap=cm.jet )
        if self.coveragePushButton.isChecked():
            self.spotcoveragePlotWidget.axes.imshow( spotcov, \
                                                         interpolation='nearest', cmap=cm.spring, alpha=.3 )
        self.spotcoveragePlotWidget.axes.set_xlim( xmin=c0,xmax=c1 )
        self.spotcoveragePlotWidget.axes.set_ylim( ymin=r0,ymax=r1 )
        self.spotcoveragePlotWidget.axes.plot( center[1], center[0], 'wx', \
                                                   markersize=10, markeredgewidth=2, alpha=1 )

        self.spotcoveragePlotWidget.fig.canvas.draw()
        
        f.close()
        f2.close()


if __name__=='__main__':
    app = QtGui.QApplication(sys.argv)
    instance = sm_human_validator(app=app)
    instance.main()
    sys.exit(app.exec_())
