# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'layout_sm_human_validator.ui'
#
# Created: Fri Apr  4 14:27:50 2014
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
        MainWindow.resize(1200, 810)
        MainWindow.setStyleSheet(_fromUtf8(""))
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.intPlotWidget = IntPlotMplCanvas(self.centralwidget, width=10, height=10, dpi=75)
        self.intPlotWidget.setGeometry(QtCore.QRect(10, 30, 1181, 161))
        self.intPlotWidget.setObjectName(_fromUtf8("intPlotWidget"))
        self.intPlotToolsWidget = NavigationToolbar(self.intPlotWidget, self.centralwidget)
        self.intPlotToolsWidget.setGeometry(QtCore.QRect(10, 0, 341, 31))
        self.intPlotToolsWidget.setObjectName(_fromUtf8("intPlotToolsWidget"))
        self.portraitPlotWidget = PortraitMplCanvas(self.centralwidget, width=10, height=10, dpi=75)
        self.portraitPlotWidget.setGeometry(QtCore.QRect(290, 220, 561, 351))
        self.portraitPlotWidget.setObjectName(_fromUtf8("portraitPlotWidget"))
        self.projectionsPlotWidget = ProjectionsMplCanvas(self.centralwidget, width=10, height=10, dpi=75)
        self.projectionsPlotWidget.setGeometry(QtCore.QRect(10, 580, 571, 181))
        self.projectionsPlotWidget.setObjectName(_fromUtf8("projectionsPlotWidget"))
        self.yarpPushButton = QtGui.QPushButton(self.centralwidget)
        self.yarpPushButton.setGeometry(QtCore.QRect(890, 670, 101, 101))
        self.yarpPushButton.setAutoFillBackground(False)
        self.yarpPushButton.setStyleSheet(_fromUtf8("background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 rgba(255, 255, 255, 255), stop:1 rgba(100, 255, 100, 255));\n"
"font-weight: bold;\n"
"font-size: 20pt;\n"
"color: green;\n"
""))
        self.yarpPushButton.setObjectName(_fromUtf8("yarpPushButton"))
        self.narpPushButton = QtGui.QPushButton(self.centralwidget)
        self.narpPushButton.setGeometry(QtCore.QRect(1070, 670, 101, 101))
        self.narpPushButton.setAutoFillBackground(False)
        self.narpPushButton.setStyleSheet(_fromUtf8("background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 rgba(255, 255, 255, 255), stop:1 rgba(255, 100, 100, 255));\n"
"font-weight: bold;\n"
"font-size: 20pt;\n"
"color: red;\n"
""))
        self.narpPushButton.setObjectName(_fromUtf8("narpPushButton"))
        self.ProjectionsPlotToolsWidget = NavigationToolbar(self.projectionsPlotWidget, self.centralwidget)
        self.ProjectionsPlotToolsWidget.setGeometry(QtCore.QRect(10, 760, 571, 31))
        self.ProjectionsPlotToolsWidget.setObjectName(_fromUtf8("ProjectionsPlotToolsWidget"))
        self.ETmodelPlotWidget = ETmodelMplCanvas(self.centralwidget, width=10, height=10, dpi=75)
        self.ETmodelPlotWidget.setGeometry(QtCore.QRect(660, 580, 191, 191))
        self.ETmodelPlotWidget.setObjectName(_fromUtf8("ETmodelPlotWidget"))
        self.spotcoveragePlotWidget = spotcoverageMplCanvas(self.centralwidget, width=10, height=10, dpi=75)
        self.spotcoveragePlotWidget.setGeometry(QtCore.QRect(860, 220, 331, 331))
        self.spotcoveragePlotWidget.setObjectName(_fromUtf8("spotcoveragePlotWidget"))
        self.meanintPushButton = QtGui.QPushButton(self.spotcoveragePlotWidget)
        self.meanintPushButton.setGeometry(QtCore.QRect(0, 0, 71, 21))
        self.meanintPushButton.setStyleSheet(_fromUtf8(""))
        self.meanintPushButton.setCheckable(True)
        self.meanintPushButton.setChecked(True)
        self.meanintPushButton.setObjectName(_fromUtf8("meanintPushButton"))
        self.coveragePushButton = QtGui.QPushButton(self.spotcoveragePlotWidget)
        self.coveragePushButton.setGeometry(QtCore.QRect(260, 0, 71, 21))
        self.coveragePushButton.setStyleSheet(_fromUtf8(""))
        self.coveragePushButton.setCheckable(True)
        self.coveragePushButton.setChecked(True)
        self.coveragePushButton.setObjectName(_fromUtf8("coveragePushButton"))
        self.spotcoveragePlotToolsWidget = NavigationToolbar(self.spotcoveragePlotWidget, self.centralwidget)
        self.spotcoveragePlotToolsWidget.setGeometry(QtCore.QRect(860, 550, 291, 31))
        self.spotcoveragePlotToolsWidget.setObjectName(_fromUtf8("spotcoveragePlotToolsWidget"))
        self.imlostPushButton = QtGui.QPushButton(self.centralwidget)
        self.imlostPushButton.setGeometry(QtCore.QRect(890, 610, 281, 51))
        self.imlostPushButton.setAutoFillBackground(True)
        self.imlostPushButton.setStyleSheet(_fromUtf8(""))
        self.imlostPushButton.setCheckable(True)
        self.imlostPushButton.setObjectName(_fromUtf8("imlostPushButton"))
        self.label = QtGui.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(360, 0, 831, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label.setFont(font)
        self.label.setObjectName(_fromUtf8("label"))
        self.residualsPlotWidget = residualsMplCanvas(self.centralwidget, width=10, height=10, dpi=75)
        self.residualsPlotWidget.setGeometry(QtCore.QRect(10, 220, 281, 351))
        self.residualsPlotWidget.setObjectName(_fromUtf8("residualsPlotWidget"))
        self.numericPortraitsPushButton = QtGui.QPushButton(self.centralwidget)
        self.numericPortraitsPushButton.setGeometry(QtCore.QRect(520, 198, 91, 21))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.numericPortraitsPushButton.setFont(font)
        self.numericPortraitsPushButton.setAutoFillBackground(True)
        self.numericPortraitsPushButton.setStyleSheet(_fromUtf8(""))
        self.numericPortraitsPushButton.setCheckable(True)
        self.numericPortraitsPushButton.setChecked(False)
        self.numericPortraitsPushButton.setObjectName(_fromUtf8("numericPortraitsPushButton"))
        self.interpolatedPortraitsPushButton = QtGui.QPushButton(self.centralwidget)
        self.interpolatedPortraitsPushButton.setGeometry(QtCore.QRect(610, 198, 111, 21))
        font = QtGui.QFont()
        font.setKerning(True)
        self.interpolatedPortraitsPushButton.setFont(font)
        self.interpolatedPortraitsPushButton.setAutoFillBackground(True)
        self.interpolatedPortraitsPushButton.setStyleSheet(_fromUtf8(""))
        self.interpolatedPortraitsPushButton.setCheckable(True)
        self.interpolatedPortraitsPushButton.setChecked(True)
        self.interpolatedPortraitsPushButton.setObjectName(_fromUtf8("interpolatedPortraitsPushButton"))
        self.label_2 = QtGui.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(350, 200, 161, 17))
        self.label_2.setObjectName(_fromUtf8("label_2"))
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
        self.meanintPushButton.setText(QtGui.QApplication.translate("MainWindow", "mean int", None, QtGui.QApplication.UnicodeUTF8))
        self.coveragePushButton.setText(QtGui.QApplication.translate("MainWindow", "coverage", None, QtGui.QApplication.UnicodeUTF8))
        self.imlostPushButton.setText(QtGui.QApplication.translate("MainWindow", "help me, I\'m lost with all these graphs", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("MainWindow", "No files found.", None, QtGui.QApplication.UnicodeUTF8))
        self.numericPortraitsPushButton.setText(QtGui.QApplication.translate("MainWindow", "numeric", None, QtGui.QApplication.UnicodeUTF8))
        self.interpolatedPortraitsPushButton.setText(QtGui.QApplication.translate("MainWindow", "interpolated", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("MainWindow", "Type of portait shown: ", None, QtGui.QApplication.UnicodeUTF8))

