import sys, os
from PyQt4 import QtGui, QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import numpy as np

class IntPlotMplCanvas(FigureCanvas):
    """Simple canvas with a sine plot."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        matplotlib.rcParams.update( {'font size': 9} )
        self.axes_ints = self.fig.add_subplot(1,1,1)
        self.fig.subplots_adjust(left=0, bottom=0, top=1, right=1, wspace=0, hspace=0)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def clear(self):
        self.axes_ints.clear()
        


class PortraitMplCanvas(FigureCanvas):
    """Simple canvas with a sine plot."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        matplotlib.rcParams.update( {'font size': 9} )
        # self.fig = Figure(figsize=(width, height), dpi=dpi)
        # self.axes = self.fig.add_subplot(1,1,1)
        self.fig, self.axes = plt.subplots( 2,3, figsize=(width, height), dpi=dpi, frameon=False ) 
        self.fig.subplots_adjust(left=0., bottom=0., top=1, right=1, wspace=0.01, hspace=0.01)
        for a in self.axes.flatten():
            a.axis('off')

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def clear(self):
        for a in self.axes.flatten():
            a.clear()
            a.axis('off')


class ProjectionsMplCanvas(FigureCanvas):
    """Simple canvas with a sine plot."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        matplotlib.rcParams.update( {'font size': 9} )
        self.axes_mex = self.fig.add_subplot(1,2,1)
        self.axes_mem = self.fig.add_subplot(1,2,2)
        self.fig.subplots_adjust(left=0, bottom=0, top=1, right=1, wspace=0, hspace=0)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def clear(self):
        self.axes_mex.clear()
        self.axes_mem.clear()


class ETmodelMplCanvas(FigureCanvas):
    """Simple canvas with a sine plot."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        matplotlib.rcParams.update( {'font size': 9} )
        self.axes = self.fig.add_subplot(1,1,1, polar=True)
        self.fig.subplots_adjust(left=0, bottom=0, top=1, right=1, wspace=0, hspace=0)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def clear(self):
        self.axes.clear()


class spotcoverageMplCanvas(FigureCanvas):
    """Simple canvas with a sine plot."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        matplotlib.rcParams.update( {'font size': 9} )
        self.axes = self.fig.add_subplot(1,1,1)
        self.fig.subplots_adjust(left=0, bottom=0, top=1, right=1, wspace=0, hspace=0)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def clear(self):
        self.axes.clear()

class residualsMplCanvas(FigureCanvas):
    """Simple canvas with a sine plot."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        matplotlib.rcParams.update( {'font size': 9} )
        self.axes_portrait = self.fig.add_subplot(2,1,1)
        self.axes_ETportrait = self.fig.add_subplot(2,1,2)
        self.fig.subplots_adjust(left=0, bottom=0, top=1, right=1, wspace=0, hspace=0)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def clear(self):
        self.axes_portrait.clear()
        self.axes_ETportrait.clear()
