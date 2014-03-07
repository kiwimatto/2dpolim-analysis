# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'layout_sm_human_validator.ui'
#
# Created: Mon Mar  3 16:51:33 2014
#      by: PyQt4 UI code generator 4.9.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui
from smHumanValidatorCanvasses import *
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(1000, 800)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.intPlotWidget = IntPlotMplCanvas(self.centralwidget, width=10, height=10, dpi=75)
        self.intPlotWidget.setGeometry(QtCore.QRect(10, 10, 981, 311))
        self.intPlotWidget.setObjectName(_fromUtf8("intPlotWidget"))
        self.intPlotToolsWidget = NavigationToolbar(self.intPlotWidget, self.centralwidget)
        self.intPlotToolsWidget.setGeometry(QtCore.QRect(10, 320, 981, 31))
        self.intPlotToolsWidget.setObjectName(_fromUtf8("intPlotToolsWidget"))
        self.portraitPlotWidget = PortraitMplCanvas(self.centralwidget, width=10, height=10, dpi=75)
        self.portraitPlotWidget.setGeometry(QtCore.QRect(10, 370, 371, 371))
        self.portraitPlotWidget.setObjectName(_fromUtf8("portraitPlotWidget"))
        self.projectionsPlotWidget = ProjectionsMplCanvas(self.centralwidget, width=10, height=10, dpi=75)
        self.projectionsPlotWidget.setGeometry(QtCore.QRect(410, 370, 371, 371))
        self.projectionsPlotWidget.setObjectName(_fromUtf8("projectionsPlotWidget"))
        self.yarpPushButton = QtGui.QPushButton(self.centralwidget)
        self.yarpPushButton.setGeometry(QtCore.QRect(810, 370, 161, 101))
        self.yarpPushButton.setStyleSheet(_fromUtf8("background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 rgba(255, 255, 255, 255), stop:1 rgba(100, 255, 100, 255));\n"
"font-weight: bold;\n"
"font-size: 20pt;\n"
"color: green;\n"
""))
        self.yarpPushButton.setObjectName(_fromUtf8("yarpPushButton"))
        self.narpPushButton = QtGui.QPushButton(self.centralwidget)
        self.narpPushButton.setGeometry(QtCore.QRect(810, 480, 161, 101))
        self.narpPushButton.setStyleSheet(_fromUtf8("background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 rgba(255, 255, 255, 255), stop:1 rgba(255, 100, 100, 255));\n"
"font-weight: bold;\n"
"font-size: 20pt;\n"
"color: red;\n"
""))
        self.narpPushButton.setObjectName(_fromUtf8("narpPushButton"))
        self.PortraitPlotToolsWidget = NavigationToolbar(self.portraitPlotWidget, self.centralwidget)
        self.PortraitPlotToolsWidget.setGeometry(QtCore.QRect(10, 740, 371, 31))
        self.PortraitPlotToolsWidget.setObjectName(_fromUtf8("PortraitPlotToolsWidget"))
        self.ProjectionsPlotToolsWidget = NavigationToolbar(self.projectionsPlotWidget, self.centralwidget)
        self.ProjectionsPlotToolsWidget.setGeometry(QtCore.QRect(410, 740, 371, 31))
        self.ProjectionsPlotToolsWidget.setObjectName(_fromUtf8("ProjectionsPlotToolsWidget"))
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.yarpPushButton.setText(QtGui.QApplication.translate("MainWindow", "YARP", None, QtGui.QApplication.UnicodeUTF8))
        self.narpPushButton.setText(QtGui.QApplication.translate("MainWindow", "NARP", None, QtGui.QApplication.UnicodeUTF8))

