# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gfmeshgui.ui'
#
# Created by: PyQt5 UI code generator 5.8.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from gfmesh_20170306 import makeGFGrid

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 600)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.Generate_btn = QtWidgets.QPushButton(self.centralwidget)
        self.Generate_btn.setGeometry(QtCore.QRect(570, 460, 89, 25))
        self.Generate_btn.setObjectName("Generate_btn")
        self.xnum_tedit = QtWidgets.QLineEdit(self.centralwidget)
        self.xnum_tedit.setGeometry(QtCore.QRect(220, 70, 113, 25))
        self.xnum_tedit.setObjectName("xnum_tedit")
        self.ynum_tedit = QtWidgets.QLineEdit(self.centralwidget)
        self.ynum_tedit.setGeometry(QtCore.QRect(220, 130, 113, 25))
        self.ynum_tedit.setObjectName("ynum_tedit")
        self.znum_tedit = QtWidgets.QLineEdit(self.centralwidget)
        self.znum_tedit.setGeometry(QtCore.QRect(220, 190, 113, 25))
        self.znum_tedit.setObjectName("znum_tedit")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(180, 70, 31, 21))
        font = QtGui.QFont()
        font.setPointSize(16)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(180, 120, 31, 41))
        font = QtGui.QFont()
        font.setPointSize(16)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(180, 190, 31, 21))
        font = QtGui.QFont()
        font.setPointSize(16)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 22))
        self.menubar.setObjectName("menubar")
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.Generate_btn.setText(_translate("MainWindow", "Generate"))
        self.label.setText(_translate("MainWindow", "x="))
        self.label_2.setText(_translate("MainWindow", "y="))
        self.label_3.setText(_translate("MainWindow", "z="))



class mywindow(QtWidgets.QWidget, Ui_MainWindow):
    global mStr

    def __init__(self):
        super(mywindow, self).__init__()
        self.setupUi(self)
        self.Generate_btn.clicked.connect(self.makegfmesh)

    def makegfmesh(self):
        xnum = int(self.xnum_tedit.text())
        ynum = int(self.ynum_tedit.text())
        znum = int(self.znum_tedit.text())
        print("x=%d, y=%d, z=%d" %(xnum,ynum,znum))
        makeGFGrid("./bCtl.stl", xnum, ynum, znum)

if __name__=="__main__":
    import sys
    app=QtWidgets.QApplication(sys.argv)
    gfmeshgui=mywindow()
    gfmeshgui.show()
    sys.exit(app.exec_())

