# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'analyse_sm.ui'
#
# Created: Mon Mar  3 09:38:35 2014
#      by: PyQt4 UI code generator 4.9.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(384, 365)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.formLayoutWidget = QtGui.QWidget(self.centralwidget)
        self.formLayoutWidget.setGeometry(QtCore.QRect(10, 10, 361, 333))
        self.formLayoutWidget.setObjectName(_fromUtf8("formLayoutWidget"))
        self.formLayout = QtGui.QFormLayout(self.formLayoutWidget)
        self.formLayout.setSizeConstraint(QtGui.QLayout.SetDefaultConstraint)
        self.formLayout.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow)
        self.formLayout.setLabelAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.formLayout.setMargin(5)
        self.formLayout.setMargin(0)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.selectDirectoryPushButton = QtGui.QPushButton(self.formLayoutWidget)
        self.selectDirectoryPushButton.setObjectName(_fromUtf8("selectDirectoryPushButton"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.SpanningRole, self.selectDirectoryPushButton)
        self.selectSPEComboBox = QtGui.QComboBox(self.formLayoutWidget)
        self.selectSPEComboBox.setObjectName(_fromUtf8("selectSPEComboBox"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.SpanningRole, self.selectSPEComboBox)
        self.motorFileLabel = QtGui.QLabel(self.formLayoutWidget)
        self.motorFileLabel.setText(_fromUtf8(""))
        self.motorFileLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.motorFileLabel.setObjectName(_fromUtf8("motorFileLabel"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.SpanningRole, self.motorFileLabel)
        self.blankFileLabel = QtGui.QLabel(self.formLayoutWidget)
        self.blankFileLabel.setText(_fromUtf8(""))
        self.blankFileLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.blankFileLabel.setObjectName(_fromUtf8("blankFileLabel"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.SpanningRole, self.blankFileLabel)
        self.label_3 = QtGui.QLabel(self.formLayoutWidget)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.formLayout.setWidget(5, QtGui.QFormLayout.LabelRole, self.label_3)
        self.label_4 = QtGui.QLabel(self.formLayoutWidget)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.formLayout.setWidget(6, QtGui.QFormLayout.LabelRole, self.label_4)
        self.label_5 = QtGui.QLabel(self.formLayoutWidget)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.formLayout.setWidget(7, QtGui.QFormLayout.LabelRole, self.label_5)
        self.spotTypeComboBox = QtGui.QComboBox(self.formLayoutWidget)
        self.spotTypeComboBox.setObjectName(_fromUtf8("spotTypeComboBox"))
        self.spotTypeComboBox.addItem(_fromUtf8(""))
        self.spotTypeComboBox.addItem(_fromUtf8(""))
        self.formLayout.setWidget(7, QtGui.QFormLayout.FieldRole, self.spotTypeComboBox)
        self.spotSizeLabel = QtGui.QLabel(self.formLayoutWidget)
        self.spotSizeLabel.setObjectName(_fromUtf8("spotSizeLabel"))
        self.formLayout.setWidget(8, QtGui.QFormLayout.LabelRole, self.spotSizeLabel)
        self.spotSizeLabel_2 = QtGui.QLabel(self.formLayoutWidget)
        self.spotSizeLabel_2.setText(_fromUtf8(""))
        self.spotSizeLabel_2.setObjectName(_fromUtf8("spotSizeLabel_2"))
        self.formLayout.setWidget(9, QtGui.QFormLayout.FieldRole, self.spotSizeLabel_2)
        self.runAnalysisPushButton = QtGui.QPushButton(self.formLayoutWidget)
        self.runAnalysisPushButton.setObjectName(_fromUtf8("runAnalysisPushButton"))
        self.formLayout.setWidget(10, QtGui.QFormLayout.FieldRole, self.runAnalysisPushButton)
        self.spotCoordsFileLabel = QtGui.QLabel(self.formLayoutWidget)
        self.spotCoordsFileLabel.setText(_fromUtf8(""))
        self.spotCoordsFileLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.spotCoordsFileLabel.setObjectName(_fromUtf8("spotCoordsFileLabel"))
        self.formLayout.setWidget(4, QtGui.QFormLayout.SpanningRole, self.spotCoordsFileLabel)
        self.spotSizeSpinBox = QtGui.QDoubleSpinBox(self.formLayoutWidget)
        self.spotSizeSpinBox.setProperty("value", 6.0)
        self.spotSizeSpinBox.setObjectName(_fromUtf8("spotSizeSpinBox"))
        self.formLayout.setWidget(8, QtGui.QFormLayout.FieldRole, self.spotSizeSpinBox)
        self.SNRSpinBox = QtGui.QDoubleSpinBox(self.formLayoutWidget)
        self.SNRSpinBox.setSingleStep(0.1)
        self.SNRSpinBox.setProperty("value", 2.0)
        self.SNRSpinBox.setObjectName(_fromUtf8("SNRSpinBox"))
        self.formLayout.setWidget(5, QtGui.QFormLayout.FieldRole, self.SNRSpinBox)
        self.VFRSpinBox = QtGui.QDoubleSpinBox(self.formLayoutWidget)
        self.VFRSpinBox.setSingleStep(0.1)
        self.VFRSpinBox.setProperty("value", 0.5)
        self.VFRSpinBox.setObjectName(_fromUtf8("VFRSpinBox"))
        self.formLayout.setWidget(6, QtGui.QFormLayout.FieldRole, self.VFRSpinBox)
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.selectDirectoryPushButton.setText(QtGui.QApplication.translate("MainWindow", "select data directory", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("MainWindow", "SNR", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("MainWindow", "VFR", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("MainWindow", " spot type", None, QtGui.QApplication.UnicodeUTF8))
        self.spotTypeComboBox.setItemText(0, QtGui.QApplication.translate("MainWindow", "circle", None, QtGui.QApplication.UnicodeUTF8))
        self.spotTypeComboBox.setItemText(1, QtGui.QApplication.translate("MainWindow", "square", None, QtGui.QApplication.UnicodeUTF8))
        self.spotSizeLabel.setText(QtGui.QApplication.translate("MainWindow", "edge length", None, QtGui.QApplication.UnicodeUTF8))
        self.runAnalysisPushButton.setText(QtGui.QApplication.translate("MainWindow", "run analysis", None, QtGui.QApplication.UnicodeUTF8))

