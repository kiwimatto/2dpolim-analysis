# -*- coding: utf-8 -*-

import sys, os
from PyQt4 import QtGui, QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
import matplotlib

import numpy as np

import spot_picker

# def format_coord(x, y):
#     col = int(x)
#     row = int(y)
#     if not self.image==None:
#         numrows, numcols = self.image.shape
#         if col>=0 and col<numcols and row>=0 and row<numrows:
#             z = self.image[row,col]
#             return 'x=%1.4f, y=%1.4f, z=%1.4f' % (x, y, z)
#         else:
#             return 'x=%1.4f, y=%1.4f' % (x, y)

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

        # self.image = None    # needed in format_coord
        # self.myaxes.format_coord = format_coord

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
        self.myaxes.figure.canvas.mpl_connect('key_press_event', self.on_key)
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


    def show_image(self, image, origin='upper', zorder=1, cmap=cm.gray):
        self.myaxes.cla()
        self.image = image
        # self.im = self.myaxes.imshow( image, origin=origin, \
        #                          zorder=zorder, cmap=cmap, \
        #                          interpolation='nearest', alpha=1 )
        self.im = new_imshow( self.myaxes, range(image.shape[1]), range(image.shape[0]), image, \
                                  origin=origin, zorder=zorder, cmap=cmap, \
                                  interpolation='nearest', alpha=1 )

        self.cbar = self.fig.colorbar(self.im, cax=self.fig.axes[1] )
        self.fig.axes[0].set_position( [.05, .1, .8, .8] )
        self.fig.axes[1].set_position( [.9,.1,.05,.8] )
        self.figure.canvas.draw()
        self.show_stuff(what='')

        # print 'artists in mymplcanvas:'
        # print self.myaxes.artists
        # for a in self.myaxes.artists:
        #     print a," --> ",a.get_gid()

    def show_stuff(self, what='spots' ):        
        
        print self.myaxes.artists

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
            mex_im = self.myaxes.imshow( self.parent().parent().m.M_ex_image, \
                                             alpha=.7, zorder=2, interpolation='nearest' )
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





class imshow_show_z:

    def __init__(self, ax, z, x, y):
        self.ax = ax
        self.x  = x
        self.y  = y
        self.z  = z
        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.numrows, self.numcols = self.z.shape
        self.ax.format_coord = self.format_coord
    def format_coord(self, x, y):
        col = int(x/self.dx+0.5)
        row = int(y/self.dy+0.5)
        #print "Nx, Nf = ", len(self.x), len(self.y), "    x, y =", x, y, "    dx, dy =", self.dx, self.dy, "    col, row =", col, row
        xyz_str = ''
        if ((col>=0) and (col<self.numcols) and (row>=0) and (row<self.numrows)):
            zij = self.z[row,col]
            #print "zij =", zij, '  |zij| =', abs(zij)
            if (np.iscomplex(zij)):
                amp = abs(zij)
                phs = np.angle(zij) / np.pi
                if (zij.imag >= 0.0):
                    signz = '+'
                else:
                    signz = '-'
                xyz_str = 'x=' + str('%.4g' % x) + ', y=' + str('%.4g' % y) + ',' \
                            + ' z=(' + str('%.4g' % zij.real) + signz + str('%.4g' % abs(zij.imag)) + 'j)' \
                            + '=' + str('%.4g' % amp) + r'*exp{' + str('%.4g' % phs) + u' π j})'
            else:
                xyz_str = 'r=' + str('%.4g' % y) + '  c=' + str('%.4g' % x) + ', z=' + str('%.4g' % zij)
        else:
            xyz_str = 'x=%1.4f, y=%1.4f'%(x, y)
#        print xyz_str
        return xyz_str

def new_imshow(ax, x, y, z, *args, **kwargs):
    assert(len(x) == z.shape[1])
    assert(len(y) == z.shape[0])
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    if (np.iscomplex(z).any()):
        zabs = abs(z)
    else:
        zabs = z
    # Use this to center pixel around (x,y) values
#    extent = (x[0]-dx/2.0, x[-1]+dx/2.0, y[-1]-dy/2.0, y[0]+dy/2.0)
    # Use this to let (x,y) be the lower-left pixel location (upper-left when origin = 'lower' is not used)
    #extent = (x[0]-dx/2.0, x[-1]+dx/2.0, y[0]-dy/2.0, y[-1]+dy/2.0)
#    im = ax.imshow(zabs, extent = extent, *args, **kwargs)
    im = ax.imshow(zabs, *args, **kwargs)
    imshow_show_z(ax, z, x, y)
#    ax.set_xlim((x[0], x[-1]))
#    ax.set_ylim((y[0], y[-1]))
    return im
