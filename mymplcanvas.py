import sys, os
from PyQt4 import QtGui, QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import matplotlib.cm as cm

import numpy as np

import spot_picker

class MyMplCanvas(FigureCanvas):
    """Simple canvas with a sine plot."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        i = self.axes.imshow(np.outer( np.linspace(0,1,10),np.linspace(0,1,10) ), zorder=1 )
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
        self.hline = self.axes.axhline(0, color='w', ls=':')
        self.vline = self.axes.axvline(0, color='w', ls=':')

        self.bg_rect = Rectangle((0,0), 0, 0, facecolor='blue', edgecolor='blue', alpha=.3, zorder=8 )
        self.axes.add_patch(self.bg_rect)

        # self.signal_rect = Rectangle((0,0), 0, 0, facecolor='red', edgecolor='red', alpha=.3, zorder=9 )
        # self.axes.add_patch(self.signal_rect)

        self.M_ex_rects = []
        self.M_em_rects = []
        self.phase_ex_rects = []
        self.phase_em_rects = []
        self.spot_rects = []

#        self.anno = spot_picker.Annotate(self.axes)


        self.rect = Rectangle((0,0), 0, 0, facecolor='yellow', edgecolor='yellow', alpha=.5, \
                                  zorder=10, clip_on=False)
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.is_pressed = False
        self.axes.add_patch(self.rect)
        self.axes.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.axes.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.axes.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.axes.figure.canvas.mpl_connect('key_press_event', self.on_key)
        self.draw()

    def on_key(self, event):
        print event
        print event.key
        # if self.crosshairs_active:
        #     if event.key=='up':
        #         self.           
            

    def on_press(self, event):
        self.is_pressed = True
        print 'press'
        self.x0 = event.xdata
        self.y0 = event.ydata    
        self.x1 = event.xdata
        self.y1 = event.ydata
        if not np.any( [c==None for c in [self.x0,self.y0,self.x1,self.y1]] ):
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            self.rect.set_linestyle('dashed')
        self.axes.figure.canvas.draw()
    def on_motion(self,event):
        if self.is_pressed is True:
            self.x1 = event.xdata
            self.y1 = event.ydata
            if not np.any( [c==None for c in [self.x0,self.y0,self.x1,self.y1]] ):
                self.rect.set_width(self.x1 - self.x0)
                self.rect.set_height(self.y1 - self.y0)
                self.rect.set_xy((self.x0, self.y0))
                self.rect.set_linestyle('dashed')
                self.axes.figure.canvas.draw()
    def on_release(self, event):
        self.is_pressed = False
        print 'release'
        self.x1 = event.xdata
        self.y1 = event.ydata
        if not np.any( [c==None for c in [self.x0,self.y0,self.x1,self.y1]] ):
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            self.rect.set_linestyle('solid')
            self.axes.figure.canvas.draw()
        print self.x0,self.y0,self.x1,self.y1
        
        self.crosshairs_active = False
        if self.x0==self.x1 and self.y0==self.y1:
            self.crosshairs_active = True
            self.crosshairs_x = np.round(self.x0)
            self.crosshairs_y = np.round(self.y0)
            self.hline.set_ydata( np.array([self.crosshairs_y,self.crosshairs_y]) )
            self.vline.set_xdata( np.array([self.crosshairs_x,self.crosshairs_x]) )
            self.axes.figure.canvas.draw()
            self.parent().parent().crosshair_pick()


    def clear(self):
        self.fig.clear()
        self.axes = self.fig.add_subplot(111)

    def show_image(self,image, origin='upper', zorder=1, cmap=cm.gray):
        # self.fig.clear()
        # self.axes = self.fig.add_subplot(111)
        self.axes.imshow( image, origin=origin, zorder=zorder, cmap=cmap )
        self.figure.canvas.draw()

    def show_stuff(self, what='spots' ):
        while len(self.axes.patches)>0:
            self.axes.patches[0].remove()

        if what=='spots':
            for r in self.spot_rects:
                self.axes.add_patch( r )
        elif what=='M_ex':
            for r in self.M_ex_rects:
                self.axes.add_patch( r )
        elif what=='M_em':
            for r in self.M_em_rects:
                self.axes.add_patch( r )
        elif what=='phase_ex':
            for r in self.phase_ex_rects:
                self.axes.add_patch( r )
        elif what=='phase_em':
            for r in self.phase_em_rects:
                self.axes.add_patch( r )
        elif what=='ET_ruler':
            for r in self.ET_ruler_rects:
                self.axes.add_patch( r )
            
        self.figure.canvas.draw()
