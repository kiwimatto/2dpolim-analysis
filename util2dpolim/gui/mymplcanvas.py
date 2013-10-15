import sys, os
from PyQt4 import QtGui, QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
import matplotlib

import numpy as np

import spot_picker

class MyMplCanvas(FigureCanvas):
    """Simple canvas with a sine plot."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):

#        self.use_blit = False

        self.fig = Figure(figsize=(width, height), dpi=dpi)
        matplotlib.rcParams.update( {'font size': 9} )
        self.myaxes = self.fig.add_subplot(111)
        self.refimagehsv = self.myaxes.imshow(180*(np.outer( np.linspace(0,1,10),np.linspace(0,1,10) )-.5), cmap=cm.hsv, zorder=1 )
        self.refimagejet = self.myaxes.imshow(np.outer( np.linspace(0,1,10),np.linspace(0,1,10) ), cmap=cm.jet, zorder=1 )
        self.cbar = self.fig.colorbar(self.refimagejet)
        # push cbar to right edge of canvas
        self.cbaraxis = self.fig.axes[1]
        self.cbaraxis.set_position( [.9,.1,.05,.8] )
        # now size up the main plot a bit
        self.myaxes.set_position( [.05, .1, .8, .8] )
#        self.figure.canvas.draw()

        # We want the axes cleared every time plot() is called
#        self.myaxes.hold(False)

        #
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

        self.crosshairs_c = 0
        self.crosshairs_r = 0
        self.hline = self.myaxes.axhline(0, color='w', ls=':', zorder=9)
        self.vline = self.myaxes.axvline(0, color='w', ls=':', zorder=9)
        self.hline.set_gid( 'crosshair_horizontal' )
        self.vline.set_gid( 'crosshair_vertical' )

        self.bg_rect = Rectangle((0,0), 0, 0, facecolor='blue', edgecolor='blue', alpha=.15, zorder=8 )
        self.bg_rect.set_gid( 'background rectangle' )
        self.myaxes.add_artist(self.bg_rect)

        # self.signal_rect = Rectangle((0,0), 0, 0, facecolor='red', edgecolor='red', alpha=.3, zorder=9 )
        # self.myaxes.add_patch(self.signal_rect)

        self.M_ex_rects = []
        self.M_em_rects = []
        self.phase_ex_rects = []
        self.phase_em_rects = []
        self.spot_gfxrepr = []

#        self.anno = spot_picker.Annotate(self.myaxes)


        self.myrect = Rectangle((0,0), 0, 0, facecolor='yellow', edgecolor='yellow', alpha=.15, \
                                  zorder=10, clip_on=False)
        self.myrect.set_gid( 'selection rectangle' )
        self.c0 = None
        self.r0 = None
        self.c1 = None
        self.r1 = None
        self.is_pressed = False
        self.myaxes.add_artist(self.myrect)
        self.myaxes.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.myaxes.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.myaxes.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
#        self.myaxes.figure.canvas.mpl_connect('key_press_event', self.on_key)
        self.draw()

    def on_key(self, event):
        print event
        print event.key
        # if self.crosshairs_active:
        #     if event.key=='up':
        #         self.           
            

    def on_press(self, event):
        if self.parent().parent().imageview_navbar.mode=='':
            self.is_pressed = True
            print 'press'
            self.c0 = np.round( event.xdata )
            self.r0 = np.round( event.ydata )
            self.c1 = np.round( event.xdata )
            self.r1 = np.round( event.ydata )
            if not np.any( [c==None for c in [self.c0,self.r0,self.c1,self.r1]] ):
                self.myrect.set_width(self.c1 - self.c0)
                self.myrect.set_height(self.r1 - self.r0)
                self.myrect.set_xy((self.c0, self.r0))
                self.myrect.set_linestyle('dashed')
                self.myaxes.figure.canvas.draw()
    def on_motion(self,event):
        if self.parent().parent().imageview_navbar.mode=='':
            if self.is_pressed is True:
                self.c1 = np.round( event.xdata )
                self.r1 = np.round( event.ydata )
                if not np.any( [c==None for c in [self.c0,self.r0,self.c1,self.r1]] ):
                    self.myrect.set_width(self.c1 - self.c0)
                    self.myrect.set_height(self.r1 - self.r0)
                    self.myrect.set_xy((self.c0, self.r0))
                    self.myrect.set_linestyle('dashed')
                    self.myaxes.figure.canvas.draw()
    def on_release(self, event):
        if self.parent().parent().imageview_navbar.mode=='':
            self.is_pressed = False
            print 'release'
            self.c1 = np.round( event.xdata )
            self.r1 = np.round( event.ydata )
            if self.c0 > self.c1:
                self.c0,self.c1=self.c1,self.c0
            if self.r0 > self.r1:
                self.r0,self.r1=self.r1,self.r0
            if not np.any( [c==None for c in [self.c0,self.r0,self.c1,self.r1]] ):
                self.myrect.set_width(self.c1 - self.c0)
                self.myrect.set_height(self.r1 - self.r0)
                self.myrect.set_xy((self.c0, self.r0))
                self.myrect.set_linestyle('solid')
                self.myaxes.figure.canvas.draw()
            print self.c0,self.r0,self.c1,self.r1
        
            self.crosshairs_active = False
            if self.c0==self.c1 and self.r0==self.r1:
                self.crosshairs_active = True
                self.crosshairs_c = np.round(self.c0)
                self.crosshairs_r = np.round(self.r0)
                # self.hline.set_ydata( np.array([self.crosshairs_y,self.crosshairs_y]) )
                # self.vline.set_xdata( np.array([self.crosshairs_x,self.crosshairs_x]) )
                self.parent().parent().crosshair_pick()


    def show_image(self,image, origin='upper', zorder=1, cmap=cm.gray):
        self.myaxes.cla()

        im = self.myaxes.imshow( image, origin=origin, \
                                 zorder=zorder, cmap=cmap, \
                                 interpolation='nearest', alpha=1 )

        self.cbar = self.fig.colorbar(im, cax=self.fig.axes[1] )
        self.fig.axes[0].set_position( [.05, .1, .8, .8] )
        self.fig.axes[1].set_position( [.9,.1,.05,.8] )
        self.figure.canvas.draw()
        self.show_stuff(what='')

        # print 'artists in mymplcanvas:'
        # print self.myaxes.artists
        # for a in self.myaxes.artists:
        #     print a," --> ",a.get_gid()

    def show_stuff(self, what='spots' ):        

        for g in self.spot_gfxrepr:
            if not g.get_axes()==None:
                g.remove()

        self.myaxes.add_artist( self.hline )
        self.myaxes.add_artist( self.vline )
        self.myaxes.add_artist( self.myrect )
        self.myaxes.add_artist( self.bg_rect )

        if what=='spots':
            for r in self.spot_gfxrepr:
                self.myaxes.add_artist( r )
                self.myaxes.draw_artist( r )
#            self.cbar = self.fig.colorbar(self.refimagejet, cax=self.cbaraxis )
        elif what=='M_ex':
            mex_im = self.myaxes.imshow( self.parent().parent().m.M_ex_image, alpha=1, zorder=2 )
            self.cbar = self.fig.colorbar(mex_im, cax=self.fig.axes[1] )
#            self.cbar = self.fig.colorbar(self.refimagejet, cax=self.cbaraxis )
        elif what=='M_em':
            for r in self.M_em_rects:
                self.myaxes.add_patch( r )
#            self.cbar = self.fig.colorbar(self.refimagejet, cax=self.cbaraxis )
        elif what=='phase_ex':
            for r in self.phase_ex_rects:
                self.myaxes.add_patch( r )
#            self.cbar = self.fig.colorbar(self.refimagehsv, cax=self.cbaraxis )
        elif what=='phase_em':
            for r in self.phase_em_rects:
                self.myaxes.add_patch( r )
#            self.cbar = self.fig.colorbar(self.refimagehsv, cax=self.cbaraxis )
        elif what=='ET_ruler':
            for r in self.ET_ruler_rects:
                self.myaxes.add_patch( r )
#            self.cbar = self.fig.colorbar(self.refimagejet, cax=self.cbaraxis )

#         self.hline.set_visible(False)
#         self.vline.set_visible(False)
            
        self.figure.canvas.draw()

# #        if self.use_blit:
# #            self.blitbackground = self.figure.canvas.copy_from_bbox(self.myaxes.bbox)

#         self.hline.set_visible(True)
#         self.vline.set_visible(True)

#         self.myaxes.draw_artist(self.hline)
#         self.myaxes.draw_artist(self.vline)
