# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'customizeMeshSizeDialog.ui'
#
# Created by: PyQt5 UI code generator 5.8.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(400, 300)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog)
        self.buttonBox.setGeometry(QtCore.QRect(30, 240, 341, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.xMinCoord_LEdit = QtWidgets.QLineEdit(Dialog)
        self.xMinCoord_LEdit.setGeometry(QtCore.QRect(58, 51, 113, 25))
        self.xMinCoord_LEdit.setObjectName("xMinCoord_LEdit")
        self.xMaxCoord_LEdit = QtWidgets.QLineEdit(Dialog)
        self.xMaxCoord_LEdit.setGeometry(QtCore.QRect(230, 50, 113, 25))
        self.xMaxCoord_LEdit.setObjectName("xMaxCoord_LEdit")
        self.yMinCoord_LEdit = QtWidgets.QLineEdit(Dialog)
        self.yMinCoord_LEdit.setGeometry(QtCore.QRect(60, 110, 113, 25))
        self.yMinCoord_LEdit.setObjectName("yMinCoord_LEdit")
        self.yMaxCoord_LEdit = QtWidgets.QLineEdit(Dialog)
        self.yMaxCoord_LEdit.setGeometry(QtCore.QRect(230, 110, 113, 25))
        self.yMaxCoord_LEdit.setObjectName("yMaxCoord_LEdit")
        self.zMinCoord_LEdit = QtWidgets.QLineEdit(Dialog)
        self.zMinCoord_LEdit.setGeometry(QtCore.QRect(60, 170, 113, 25))
        self.zMinCoord_LEdit.setObjectName("zMinCoord_LEdit")
        self.zMaxCoord_LEdit = QtWidgets.QLineEdit(Dialog)
        self.zMaxCoord_LEdit.setGeometry(QtCore.QRect(230, 170, 113, 25))
        self.zMaxCoord_LEdit.setObjectName("zMaxCoord_LEdit")
        self.label = QtWidgets.QLabel(Dialog)
        self.label.setGeometry(QtCore.QRect(100, 20, 31, 17))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(Dialog)
        self.label_2.setGeometry(QtCore.QRect(270, 20, 41, 17))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(Dialog)
        self.label_3.setGeometry(QtCore.QRect(20, 50, 16, 17))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(Dialog)
        self.label_4.setGeometry(QtCore.QRect(20, 110, 16, 21))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.label_5 = QtWidgets.QLabel(Dialog)
        self.label_5.setGeometry(QtCore.QRect(20, 170, 16, 17))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")

        self.retranslateUi(Dialog)
        self.buttonBox.accepted.connect(Dialog.accept)
        self.buttonBox.rejected.connect(Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "customize Mesh Size"))
        self.label.setText(_translate("Dialog", "Min"))
        self.label_2.setText(_translate("Dialog", "Max"))
        self.label_3.setText(_translate("Dialog", "x"))
        self.label_4.setText(_translate("Dialog", "y"))
        self.label_5.setText(_translate("Dialog", "z"))

