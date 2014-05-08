import sys, os
from PyQt4 import QtGui, QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle, Circle
import matplotlib.cm as cm
import matplotlib

import numpy as np

class MyMplCanvas(FigureCanvas):
    """Simple canvas with a sine plot."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        matplotlib.rcParams.update( {'font size': 9} )
        self.myaxes = self.fig.add_subplot(111)
        self.refimagejet = self.myaxes.imshow( np.random.random(size=(512,512)), cmap=cm.jet )
        self.cbar = self.fig.colorbar(self.refimagejet)
        # push cbar to right edge of canvas
        self.cbaraxes = self.fig.axes[1]
        self.cbaraxes.set_position( [.9,.1,.05,.8] )
        # now size up the main plot a bit
        self.myaxes.set_position( [.05, .1, .8, .8] )
#        self.figure.canvas.draw()

        # We want the axes cleared every time plot() is called
#        self.axes.hold(False)

        #
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

        self.myrect = Rectangle((0,0), 1, 1, facecolor='red', edgecolor='red', alpha=.2, zorder=8 )
        self.myaxes.add_patch(self.myrect)
        self.c0 = 10   #None
        self.r0 = 10   #None
        self.c1 = 60   #None
        self.r1 = 200  #None
        self.is_pressed = False

        self.selectionrectangle = Rectangle((0,0), 1, 1, \
                                        ls='dashed', facecolor='none', edgecolor='w', alpha=.7, zorder=9 )
        self.selectioncircle    = Circle((0,0), 1, \
                                        ls='dashed', facecolor='none', edgecolor='w', alpha=.7, zorder=9 )
        self.myaxes.add_patch(self.selectionrectangle)
        self.selectiontype = None
        self.selectionc0 = None
        self.selectionr0 = None
        self.selectionc1 = None
        self.selectionr1 = None

        self.myaxes.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.myaxes.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.myaxes.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.draw()


    def clear(self):
        self.fig.clf()
        self.myaxes = self.fig.add_subplot(111)
        self.myaxes.cla()
        self.cbaraxes.cla()
#        self.myaxes = self.fig.add_subplot(111)
        self.myrect = Rectangle((0,0), 1, 1, facecolor='red', edgecolor='red', alpha=.2, zorder=8 )
        self.myaxes.add_patch(self.myrect)
        self.selectionrectangle = Rectangle((0,0), 1, 1, \
                                         ls='dashed', facecolor='none', edgecolor='w', alpha=.5, zorder=9 )
        self.selectioncircle    = Circle((0,0), 1, \
                                         ls='dashed', facecolor='none', edgecolor='w', alpha=.5, zorder=9 )
        self.myaxes.add_patch(self.selectionrectangle)

    def set_selection_area(self, shape, c0,r0,c1,r1):
        self.selectionc0 = c0
        self.selectionr0 = r0
        self.selectionc1 = c1
        self.selectionr1 = r1
        if shape=='Rectangle':
            if self.selectionc0 > self.selectionc1:
                self.selectionc0,self.selectionc1=self.selectionc1,self.selectionc0
            if self.selectionr0 > self.selectionr1:
                self.selectionr0,self.selectionr1=self.selectionr1,self.selectionr0
            if not np.any( [c==None for c in [self.selectionc0,self.selectionr0,self.selectionc1,self.selectionr1]] ):
                self.selectionrectangle.set_width(self.selectionc1 - self.selectionc0)
                self.selectionrectangle.set_height(self.selectionr1 - self.selectionr0)
                self.selectionrectangle.set_xy((self.selectionc0, self.selectionr0))
        elif shape=='Circle':
            self.selectioncircle.center = (self.selectionc0, self.selectionr0)
            r = np.sqrt( (self.selectionc1-self.selectionc0)**2 + (self.selectionr1-self.selectionr0)**2 )
            self.selectioncircle.set_radius(r)
        else:
            raise bloodpressure
        self.myaxes.figure.canvas.draw()

    def on_press(self, event):
        self.is_pressed = True
        if self.selectiontype==None: 
            self.c0 = event.xdata
            self.r0 = event.ydata
            self.c1 = event.xdata
            self.r1 = event.ydata
            if not np.any( [c==None for c in [self.c0,self.r0,self.c1,self.r1]] ):
                self.myrect.set_width(self.c1 - self.c0)
                self.myrect.set_height(self.r1 - self.r0)
                self.myrect.set_xy((self.c0, self.r0))
                self.myrect.set_linestyle('dashed')
        elif self.selectiontype=='Rectangle':
            self.selectionc0 = event.xdata
            self.selectionr0 = event.ydata
            self.selectionc1 = event.xdata
            self.selectionr1 = event.ydata
            if not np.any( [c==None for c in [self.selectionc0,self.selectionr0,self.selectionc1,self.selectionr1]] ):
                self.selectionrectangle.set_width(self.selectionc1 - self.selectionc0)
                self.selectionrectangle.set_height(self.selectionr1 - self.selectionr0)
                self.selectionrectangle.set_xy((self.selectionc0, self.selectionr0))
        elif self.selectiontype=='Circle':
            self.selectionc0 = event.xdata
            self.selectionr0 = event.ydata
            self.selectionc1 = event.xdata
            self.selectionr1 = event.ydata
            if not np.any( [c==None for c in [self.selectionc0,self.selectionr0,self.selectionc1,self.selectionr1]] ):
                self.selectioncircle.center = (self.selectionc0, self.selectionr0)
                r = np.sqrt( (self.selectionc1-self.selectionc0)**2 + (self.selectionr1-self.selectionr0)**2 )
                self.selectioncircle.set_radius(r)
                print self.selectioncircle
        else:
            raise tax
        self.myaxes.figure.canvas.draw()

    def on_motion(self,event):
        if self.is_pressed is True:
            if self.selectiontype==None:
                self.c1 = event.xdata
                self.r1 = event.ydata
                if not np.any( [c==None for c in [self.c0,self.r0,self.c1,self.r1]] ):
                    self.myrect.set_width(self.c1 - self.c0)
                    self.myrect.set_height(self.r1 - self.r0)
                    self.myrect.set_xy((self.c0, self.r0))
                    self.myrect.set_linestyle('dashed')
            elif self.selectiontype=='Rectangle':
                self.selectionc1 = event.xdata
                self.selectionr1 = event.ydata
                if not np.any( [c==None for c in [self.selectionc0,self.selectionr0,self.selectionc1,self.selectionr1]] ):
                    self.selectionrectangle.set_width(self.selectionc1 - self.selectionc0)
                    self.selectionrectangle.set_height(self.selectionr1 - self.selectionr0)
                    self.selectionrectangle.set_xy((self.selectionc0, self.selectionr0))
            elif self.selectiontype=='Circle':
                self.selectionc1 = event.xdata
                self.selectionr1 = event.ydata
                if not np.any( [c==None for c in [self.selectionc0,self.selectionr0,self.selectionc1,self.selectionr1]] ):
                    r = np.sqrt( (self.selectionc1-self.selectionc0)**2 \
                                     + (self.selectionr1-self.selectionr0)**2 )
                    self.selectioncircle.set_radius(r)
                    print self.selectioncircle
            else:
                raise wage
            self.myaxes.figure.canvas.draw()

    def on_release(self, event):
        self.is_pressed = False
        if self.selectiontype==None:
            self.c1 = event.xdata
            self.r1 = event.ydata
            if self.c0 > self.c1:
                self.c0,self.c1=self.c1,self.c0
            if self.r0 > self.r1:
                self.r0,self.r1=self.r1,self.r0
            if not np.any( [c==None for c in [self.c0,self.r0,self.c1,self.r1]] ):
                self.myrect.set_width(self.c1 - self.c0)
                self.myrect.set_height(self.r1 - self.r0)
                self.myrect.set_xy((self.c0, self.r0))
                self.myrect.set_linestyle('solid')
        elif self.selectiontype=='Rectangle':
            self.selectionc1 = event.xdata
            self.selectionr1 = event.ydata
            if self.selectionc0 > self.selectionc1:
                self.selectionc0,self.selectionc1=self.selectionc1,self.selectionc0
            if self.selectionr0 > self.selectionr1:
                self.selectionr0,self.selectionr1=self.selectionr1,self.selectionr0
            if not np.any( [c==None for c in [self.selectionc0,self.selectionr0,self.selectionc1,self.selectionr1]] ):
                self.selectionrectangle.set_width(self.selectionc1 - self.selectionc0)
                self.selectionrectangle.set_height(self.selectionr1 - self.selectionr0)
                self.selectionrectangle.set_xy((self.selectionc0, self.selectionr0))
        elif self.selectiontype=='Circle':
            self.selectionc1 = event.xdata
            self.selectionr1 = event.ydata
            if not np.any( [c==None for c in [self.selectionc0,self.selectionr0,self.selectionc1,self.selectionr1]] ):
                r = np.sqrt( (self.selectionc1-self.selectionc0)**2 + (self.selectionr1-self.selectionr0)**2 )
                self.selectioncircle.set_radius(r)
                print self.selectioncircle
        else:
            raise seat

        self.myaxes.figure.canvas.draw()
        #print self.c0,self.r0,self.c1,self.r1

