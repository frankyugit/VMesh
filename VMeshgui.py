# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'VMeshgui.ui'
#
# Created by: PyQt5 UI code generator 5.8.2
#
# WARNING! All changes made in this file will be lost!

import re
from PyQt5 import QtCore, QtGui, QtWidgets
import vtk
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from VMesh import Geometry, Mesh, VTKScene
import customizeMeshSizeDialog

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1268, 960)
        font = QtGui.QFont()
        font.setPointSize(11)
        MainWindow.setFont(font)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.Generate_btn = QtWidgets.QPushButton(self.centralwidget)
        self.Generate_btn.setGeometry(QtCore.QRect(20, 430, 89, 25))
        self.Generate_btn.setObjectName("Generate_btn")
        self.xnum_tedit = QtWidgets.QLineEdit(self.centralwidget)
        self.xnum_tedit.setGeometry(QtCore.QRect(60, 70, 41, 25))
        self.xnum_tedit.setObjectName("xnum_tedit")
        self.ynum_tedit = QtWidgets.QLineEdit(self.centralwidget)
        self.ynum_tedit.setGeometry(QtCore.QRect(60, 130, 41, 25))
        self.ynum_tedit.setObjectName("ynum_tedit")
        self.znum_tedit = QtWidgets.QLineEdit(self.centralwidget)
        self.znum_tedit.setGeometry(QtCore.QRect(60, 190, 41, 25))
        self.znum_tedit.setObjectName("znum_tedit")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(20, 70, 31, 21))
        font = QtGui.QFont()
        font.setPointSize(16)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(20, 120, 31, 41))
        font = QtGui.QFont()
        font.setPointSize(16)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(20, 190, 31, 21))
        font = QtGui.QFont()
        font.setPointSize(16)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.qvtkWidget = QVTKRenderWindowInteractor(self.centralwidget)
        self.qvtkWidget.setGeometry(QtCore.QRect(400, 30, 600, 600))
        self.qvtkWidget.setObjectName("qvtkWidget")
        self.PreView_btn = QtWidgets.QPushButton(self.centralwidget)
        self.PreView_btn.setGeometry(QtCore.QRect(130, 430, 89, 25))
        self.PreView_btn.setObjectName("PreView_btn")
        self.show_CAD_btn = QtWidgets.QPushButton(self.centralwidget)
        self.show_CAD_btn.setGeometry(QtCore.QRect(20, 470, 89, 25))
        self.show_CAD_btn.setObjectName("show_CAD_btn")
        self.writeframe = QtWidgets.QFrame(self.centralwidget)
        self.writeframe.setGeometry(QtCore.QRect(20, 510, 321, 80))
        self.writeframe.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.writeframe.setFrameShadow(QtWidgets.QFrame.Raised)
        self.writeframe.setLineWidth(2)
        self.writeframe.setMidLineWidth(0)
        self.writeframe.setObjectName("writeframe")
        self.write_btn = QtWidgets.QPushButton(self.writeframe)
        self.write_btn.setGeometry(QtCore.QRect(10, 30, 89, 25))
        self.write_btn.setObjectName("write_btn")
        self.of_type_rdbtn = QtWidgets.QRadioButton(self.writeframe)
        self.of_type_rdbtn.setGeometry(QtCore.QRect(130, 10, 112, 23))
        self.of_type_rdbtn.setObjectName("of_type_rdbtn")
        self.netcdf_type_rdbtn = QtWidgets.QRadioButton(self.writeframe)
        self.netcdf_type_rdbtn.setGeometry(QtCore.QRect(130, 50, 112, 23))
        self.netcdf_type_rdbtn.setObjectName("netcdf_type_rdbtn")
        self.x_REcheck_text = QtWidgets.QLabel(self.centralwidget)
        self.x_REcheck_text.setGeometry(QtCore.QRect(130, 70, 180, 20))
        self.x_REcheck_text.setText("")
        self.x_REcheck_text.setObjectName("x_REcheck_text")
        self.y_REcheck_text = QtWidgets.QLabel(self.centralwidget)
        self.y_REcheck_text.setGeometry(QtCore.QRect(130, 130, 180, 20))
        self.y_REcheck_text.setText("")
        self.y_REcheck_text.setObjectName("y_REcheck_text")
        self.z_REcheck_text = QtWidgets.QLabel(self.centralwidget)
        self.z_REcheck_text.setGeometry(QtCore.QRect(140, 190, 180, 20))
        self.z_REcheck_text.setText("")
        self.z_REcheck_text.setObjectName("z_REcheck_text")
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(20, 240, 131, 17))
        self.label_4.setObjectName("label_4")
        self.cad_file_dir_text = QtWidgets.QLineEdit(self.centralwidget)
        self.cad_file_dir_text.setGeometry(QtCore.QRect(20, 270, 351, 25))
        self.cad_file_dir_text.setObjectName("cad_file_dir_text")
        self.cad_browse_btn = QtWidgets.QPushButton(self.centralwidget)
        self.cad_browse_btn.setGeometry(QtCore.QRect(20, 310, 89, 25))
        self.cad_browse_btn.setObjectName("cad_browse_btn")
        self.py_terminal_textBrowser = QtWidgets.QTextBrowser(self.centralwidget)
        self.py_terminal_textBrowser.setGeometry(QtCore.QRect(20, 710, 981, 171))
        self.py_terminal_textBrowser.setObjectName("py_terminal_textBrowser")
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(260, 50, 91, 17))
        self.label_5.setObjectName("label_5")
        self.coord_check_x_LEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.coord_check_x_LEdit.setGeometry(QtCore.QRect(1100, 70, 50, 25))
        self.coord_check_x_LEdit.setObjectName("coord_check_x_LEdit")
        self.index_check_x_LEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.index_check_x_LEdit.setGeometry(QtCore.QRect(1200, 70, 50, 25))
        self.index_check_x_LEdit.setObjectName("index_check_x_LEdit")
        self.coord_check_y_LEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.coord_check_y_LEdit.setGeometry(QtCore.QRect(1100, 110, 50, 25))
        self.coord_check_y_LEdit.setObjectName("coord_check_y_LEdit")
        self.index_check_y_LEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.index_check_y_LEdit.setGeometry(QtCore.QRect(1200, 110, 50, 25))
        self.index_check_y_LEdit.setObjectName("index_check_y_LEdit")
        self.coord_check_z_LEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.coord_check_z_LEdit.setGeometry(QtCore.QRect(1100, 150, 50, 25))
        self.coord_check_z_LEdit.setObjectName("coord_check_z_LEdit")
        self.index_check_z_LEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.index_check_z_LEdit.setGeometry(QtCore.QRect(1200, 150, 50, 25))
        self.index_check_z_LEdit.setObjectName("index_check_z_LEdit")
        self.label_6 = QtWidgets.QLabel(self.centralwidget)
        self.label_6.setGeometry(QtCore.QRect(1030, 70, 21, 17))
        self.label_6.setObjectName("label_6")
        self.label_7 = QtWidgets.QLabel(self.centralwidget)
        self.label_7.setGeometry(QtCore.QRect(1030, 110, 21, 17))
        self.label_7.setObjectName("label_7")
        self.label_8 = QtWidgets.QLabel(self.centralwidget)
        self.label_8.setGeometry(QtCore.QRect(1030, 150, 21, 17))
        self.label_8.setObjectName("label_8")
        self.label_9 = QtWidgets.QLabel(self.centralwidget)
        self.label_9.setGeometry(QtCore.QRect(1100, 40, 51, 17))
        self.label_9.setObjectName("label_9")
        self.label_10 = QtWidgets.QLabel(self.centralwidget)
        self.label_10.setGeometry(QtCore.QRect(1200, 40, 41, 17))
        self.label_10.setObjectName("label_10")
        self.coord_to_index_btn = QtWidgets.QPushButton(self.centralwidget)
        self.coord_to_index_btn.setGeometry(QtCore.QRect(1130, 200, 89, 25))
        self.coord_to_index_btn.setObjectName("coord_to_index_btn")
        self.index_to_coord_btn = QtWidgets.QPushButton(self.centralwidget)
        self.index_to_coord_btn.setGeometry(QtCore.QRect(1130, 230, 89, 25))
        self.index_to_coord_btn.setObjectName("index_to_coord_btn")
        self.check_inner_btn = QtWidgets.QPushButton(self.centralwidget)
        self.check_inner_btn.setGeometry(QtCore.QRect(1180, 280, 71, 25))
        self.check_inner_btn.setObjectName("check_inner_btn")
        self.label_11 = QtWidgets.QLabel(self.centralwidget)
        self.label_11.setGeometry(QtCore.QRect(270, 110, 67, 17))
        self.label_11.setObjectName("label_11")
        self.customize_btn = QtWidgets.QPushButton(self.centralwidget)
        self.customize_btn.setGeometry(QtCore.QRect(250, 190, 89, 25))
        self.customize_btn.setObjectName("customize_btn")
        self.locate_inner_btn = QtWidgets.QPushButton(self.centralwidget)
        self.locate_inner_btn.setGeometry(QtCore.QRect(1090, 280, 71, 25))
        self.locate_inner_btn.setObjectName("locate_inner_btn")
        self.cad_scaling_Slider = QtWidgets.QSlider(self.centralwidget)
        self.cad_scaling_Slider.setGeometry(QtCore.QRect(230, 70, 91, 20))
        self.cad_scaling_Slider.setOrientation(QtCore.Qt.Horizontal)
        self.cad_scaling_Slider.setObjectName("cad_scaling_Slider")
        self.cad_scaling_SpinBox = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.cad_scaling_SpinBox.setGeometry(QtCore.QRect(330, 70, 51, 26))
        self.cad_scaling_SpinBox.setObjectName("cad_scaling_SpinBox")
        self.margin_Slider = QtWidgets.QSlider(self.centralwidget)
        self.margin_Slider.setGeometry(QtCore.QRect(230, 140, 91, 16))
        self.margin_Slider.setOrientation(QtCore.Qt.Horizontal)
        self.margin_Slider.setObjectName("margin_Slider")
        self.margin_SpinBox = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.margin_SpinBox.setGeometry(QtCore.QRect(330, 130, 51, 26))
        self.margin_SpinBox.setObjectName("margin_SpinBox")
        self.label_12 = QtWidgets.QLabel(self.centralwidget)
        self.label_12.setGeometry(QtCore.QRect(400, 650, 101, 17))
        self.label_12.setObjectName("label_12")
        self.mesh_opacity_Slider = QtWidgets.QSlider(self.centralwidget)
        self.mesh_opacity_Slider.setGeometry(QtCore.QRect(510, 650, 160, 16))
        self.mesh_opacity_Slider.setOrientation(QtCore.Qt.Horizontal)
        self.mesh_opacity_Slider.setObjectName("mesh_opacity_Slider")
        self.mesh_opacity_SpinBox = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.mesh_opacity_SpinBox.setGeometry(QtCore.QRect(690, 640, 69, 26))
        self.mesh_opacity_SpinBox.setObjectName("mesh_opacity_SpinBox")
        self.scan_cad_btn = QtWidgets.QPushButton(self.centralwidget)
        self.scan_cad_btn.setGeometry(QtCore.QRect(130, 310, 89, 25))
        self.scan_cad_btn.setObjectName("scan_cad_btn")
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1268, 22))
        self.menubar.setObjectName("menubar")
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "VMesh"))
        self.Generate_btn.setText(_translate("MainWindow", "Generate"))
        self.label.setText(_translate("MainWindow", "x="))
        self.label_2.setText(_translate("MainWindow", "y="))
        self.label_3.setText(_translate("MainWindow", "z="))
        self.PreView_btn.setText(_translate("MainWindow", "PreView"))
        self.show_CAD_btn.setText(_translate("MainWindow", "Show CAD"))
        self.write_btn.setText(_translate("MainWindow", "Write"))
        self.of_type_rdbtn.setText(_translate("MainWindow", "0-1 array"))
        self.netcdf_type_rdbtn.setText(_translate("MainWindow", "NetCDF type"))
        self.label_4.setText(_translate("MainWindow", "CAD file directory"))
        self.cad_browse_btn.setText(_translate("MainWindow", "Browse"))
        self.label_5.setText(_translate("MainWindow", "CAD Scaling"))
        self.label_6.setText(_translate("MainWindow", "x"))
        self.label_7.setText(_translate("MainWindow", "y"))
        self.label_8.setText(_translate("MainWindow", "z"))
        self.label_9.setText(_translate("MainWindow", "Coord"))
        self.label_10.setText(_translate("MainWindow", "Index"))
        self.coord_to_index_btn.setText(_translate("MainWindow", ">>"))
        self.index_to_coord_btn.setText(_translate("MainWindow", "<<"))
        self.check_inner_btn.setText(_translate("MainWindow", "Check"))
        self.label_11.setText(_translate("MainWindow", "Margin"))
        self.customize_btn.setText(_translate("MainWindow", "Customize"))
        self.locate_inner_btn.setText(_translate("MainWindow", "Locate"))
        self.label_12.setText(_translate("MainWindow", "Mesh Opacity"))
        self.scan_cad_btn.setText(_translate("MainWindow", "Scan"))

class mywindow(QtWidgets.QWidget, Ui_MainWindow):
    global mStr

    def __init__(self):
        super(mywindow, self).__init__()
        self.setupUi(self)
        self.geo = None
        self.gfmesh = None
        self.scaling = 1.0
        self.margin = 0.0
        self.showorhideCAD = True
        self.customgeodim = {}

        self.customize_btn.clicked.connect(self.slotCustomizeMeshSize)
        self.customizeMeshSizeUi = customizeMeshSizeDialog.Ui_Dialog()

        self.cad_scaling_Slider.setMinimum(0)
        self.cad_scaling_Slider.setMaximum(1000)
        self.cad_scaling_Slider.valueChanged.connect(self.cad_scaling_SpinBox_setValue)
        self.cad_scaling_SpinBox.setDecimals(3)
        self.cad_scaling_SpinBox.setMinimum(0.000)
        self.cad_scaling_SpinBox.setMaximum(1.000)
        self.cad_scaling_SpinBox.setSingleStep(0.001)
        self.cad_scaling_SpinBox.valueChanged.connect(self.cad_scaling_Slider_setValue)

        self.margin_Slider.setMinimum(0)
        self.margin_Slider.setMaximum(100)
        self.margin_Slider.valueChanged.connect(self.margin_SpinBox_setValue)
        self.margin_SpinBox.setDecimals(2)
        self.margin_SpinBox.setMinimum(0.00)
        self.margin_SpinBox.setMaximum(1.00)
        self.margin_SpinBox.setSingleStep(0.01)
        self.margin_SpinBox.valueChanged.connect(self.margin_Slider_setValue)

        self.mesh_opacity_Slider.setMinimum(0)
        self.mesh_opacity_Slider.setMaximum(10)
        self.mesh_opacity_Slider.valueChanged.connect(self.mesh_opacity_SpinBox_setValue)
        self.mesh_opacity_SpinBox.setDecimals(1)
        self.mesh_opacity_SpinBox.setMinimum(0.0)
        self.mesh_opacity_SpinBox.setMaximum(1.0)
        self.mesh_opacity_SpinBox.setSingleStep(0.1)
        self.mesh_opacity_SpinBox.valueChanged.connect(self.mesh_opacity_Slider_setValue)

        self.coord_to_index_btn.clicked.connect(self.coordtoindex)
        self.locate_inner_btn.clicked.connect(self.renderpickupcell)
        self.check_inner_btn.clicked.connect(self.checkinnercell)

        self.Generate_btn.clicked.connect(self.makegfmesh)
        self.PreView_btn.clicked.connect(self.renderMesh)
        self.show_CAD_btn.clicked.connect(self.renderCAD)
        self.write_btn.clicked.connect(self.writemeshfile)
        self.cad_browse_btn.clicked.connect(self.openFile)
        self.qvtkWidget.Initialize()
        ren = vtk.vtkRenderer()
        self.qvtkWidget.GetRenderWindow().AddRenderer(ren)
        # sys.stdout = EmittingStream(textWritten=self.normalOutputWritten)
        # self.qvtkWidget.Start()


    def slotCustomizeMeshSize(self):
        dlg=QtWidgets.QDialog()
        self.customizeMeshSizeUi.setupUi(dlg)
        if dlg.exec_():
            self.customgeodim["xmin"] = float(self.customizeMeshSizeUi.xMinCoord_LEdit.text())
            self.customgeodim["xmax"] = float(self.customizeMeshSizeUi.xMaxCoord_LEdit.text())
            self.customgeodim["ymin"] = float(self.customizeMeshSizeUi.yMinCoord_LEdit.text())
            self.customgeodim["ymax"] = float(self.customizeMeshSizeUi.yMaxCoord_LEdit.text())
            self.customgeodim["zmin"] = float(self.customizeMeshSizeUi.zMinCoord_LEdit.text())
            self.customgeodim["zmax"] = float(self.customizeMeshSizeUi.zMaxCoord_LEdit.text())
        else:
            pass
        print(self.customgeodim)

    def cad_scaling_SpinBox_setValue(self):
        scale = self.cad_scaling_Slider.value()/1000
        self.cad_scaling_SpinBox.setValue(scale)

    def cad_scaling_Slider_setValue(self):
        scale = self.cad_scaling_SpinBox.value()*1000
        self.cad_scaling_Slider.setValue(scale)

    def margin_SpinBox_setValue(self):
        scale = self.margin_Slider.value()/100
        self.margin_SpinBox.setValue(scale)

    def margin_Slider_setValue(self):
        scale = self.margin_SpinBox.value()*100
        self.margin_Slider.setValue(scale)

    def mesh_opacity_SpinBox_setValue(self):
        scale = self.mesh_opacity_Slider.value()/10
        self.mesh_opacity_SpinBox.setValue(scale)

    def mesh_opacity_Slider_setValue(self):
        scale = self.mesh_opacity_SpinBox.value()*10
        self.mesh_opacity_Slider.setValue(scale)

    def GUIprint(self, text=""):
        """Append text to the QTextEdit."""
        # Maybe QTextEdit.append() works as well, but this is how I do it:
        cursor = self.py_terminal_textBrowser.textCursor()
        cursor.movePosition(QtGui.QTextCursor.End)
        cursor.insertText(text)
        self.py_terminal_textBrowser.setTextCursor(cursor)
        self.py_terminal_textBrowser.ensureCursorVisible()

    def makegfmesh(self):
        self.scaling = self.cad_scaling_SpinBox.value()
        self.margin = self.margin_SpinBox.value()
        print("scaling = %f" % self.scaling)
        print(" *** CAD file parsed.....")
        self.geo = Geometry(self.cad_file_dir_text.text(), self.scaling)
        self.geo.parsebinarySTLFile()
        print("  done\n")
        self.geo.printMetrics()
        nodecoords = self.genNodeCoords(self.geo)
        self.gfmesh = Mesh(nodecoords, self.geo)
        self.gfmesh.setVoxelising()

    def genNodeCoords(self, geo):
        nodecoords = [[],[],[]]
        if self.customgeodim:
            geodim = self.customgeodim
        else:
            geodim = self.geo.getgeodim()
        Nx = int(self.xnum_tedit.text())
        Ny = int(self.ynum_tedit.text())
        Nz = int(self.znum_tedit.text())
        print("Cell num: x=%d, y=%d, z=%d" %(Nx,Ny,Nz))
        NNodex = Nx + 1
        NNodey = Ny + 1
        NNodez = Nz + 1

        lengthx = geodim['xmax'] - geodim['xmin']
        meshlengthx = (1 + 2*self.margin)*lengthx
        voxelxwidth = meshlengthx/Nx
        for i in range(NNodex):
            nodecoords[0].append(geodim['xmin'] - self.margin*lengthx + voxelxwidth*i)

        lengthy = geodim['ymax'] - geodim['ymin']
        meshlengthy = (1 + 2*self.margin)*lengthy
        voxelywidth = meshlengthy/Ny
        for j in range(NNodey):
            nodecoords[1].append(geodim['ymin'] - self.margin*lengthy + voxelywidth*j)

        lengthz = geodim['zmax'] - geodim['zmin']
        meshlengthz = (1 + 2*self.margin)*lengthz
        voxelzwidth = meshlengthz/Nz
        for k in range(NNodez):
            nodecoords[2].append(geodim['zmin'] - self.margin*lengthz + voxelzwidth*k)
        return nodecoords

    def genVoxelCoords(self, geo):
        voxelcoords = [[],[],[]]
        if self.customgeodim:
            geodim = self.customgeodim
        else:
            geodim = self.geo.getgeodim()
        Nx = int(self.xnum_tedit.text())
        Ny = int(self.ynum_tedit.text())
        Nz = int(self.znum_tedit.text())
        print("x=%d, y=%d, z=%d" %(Nx,Ny,Nz))

        lengthx = geodim['xmax'] - geodim['xmin']
        meshlengthx = (1 + 2*self.margin)*lengthx
        voxelxwidth = meshlengthx/Nx
        for i in range(Nx):
            voxelcoords[0].append(geodim['xmin'] + voxelxwidth*(i+0.5) - self.margin*lengthx)

        lengthy = geodim['ymax'] - geodim['ymin']
        meshlengthy = (1 + 2*self.margin)*lengthy
        voxelywidth = meshlengthy/Ny
        for j in range(Ny):
            voxelcoords[1].append(geodim['ymin'] + voxelywidth*(j+0.5) - self.margin*lengthy)

        lengthz = geodim['zmax'] - geodim['zmin']
        meshlengthz = (1 + 2*self.margin)*lengthz
        voxelzwidth = meshlengthz/Nz
        for k in range(Nz):
            voxelcoords[2].append(geodim['zmin'] + voxelzwidth*(k+0.5) - self.margin*lengthz)
        return voxelcoords


    def renderMesh(self):
        print(self.qvtkWidget.GetRenderWindow().GetRenderers().GetFirstRenderer().GetActors())
        if self.qvtkWidget.GetRenderWindow().GetRenderers().GetFirstRenderer().GetActors().GetNumberOfItems():
            for actor in self.qvtkWidget.GetRenderWindow().GetRenderers().GetFirstRenderer().GetActors():
                self.qvtkWidget.GetRenderWindow().GetRenderers().GetFirstRenderer().RemoveActor(actor)
        print(self.qvtkWidget.GetRenderWindow().GetRenderers().GetFirstRenderer().GetActors())
        opacity = self.mesh_opacity_SpinBox.value()
        fieldlist = {"start": 1, "end": 5}
        self.vtkscene = VTKScene(self.gfmesh, fieldlist, opacity, self.qvtkWidget)
        self.vtkscene.vtkRender()
        # vtkscene.renderpickupcell(32,11,42)
        self.vtkscene.vtkwrite()

    def renderpickupcell(self):
        Xindex = int(self.index_check_x_LEdit.text())
        Yindex = int(self.index_check_y_LEdit.text())
        Zindex = int(self.index_check_z_LEdit.text())
        self.vtkscene.renderpickupcell(Xindex, Yindex, Zindex)
        self.vtkscene.vtkwrite()

    def checkinnercell(self):
        Xindex = int(self.index_check_x_LEdit.text())
        Yindex = int(self.index_check_y_LEdit.text())
        Zindex = int(self.index_check_z_LEdit.text())
        self.vtkscene.innercheck(Xindex, Yindex, Zindex)


    def renderCAD(self):
        ren = self.qvtkWidget.GetRenderWindow().GetRenderers().GetFirstRenderer()
        if self.showorhideCAD:
            reader = vtk.vtkSTLReader()
            reader.SetFileName(self.cad_file_dir_text.text())
            reader.Update()
            print(reader.GetFileName())
            scale = vtk.vtkTransform()
            scale.Scale(self.cad_scaling_SpinBox.value(), self.cad_scaling_SpinBox.value(), self.cad_scaling_SpinBox.value())
            scaledstl = vtk.vtkTransformFilter()
            scaledstl.SetInputConnection(reader.GetOutputPort())
            scaledstl.SetTransform(scale)
            sMapper = vtk.vtkDataSetMapper()
            sMapper.SetInputConnection(scaledstl.GetOutputPort())
            sActor = vtk.vtkLODActor()
            sActor.SetMapper(sMapper)
            sActor.GetProperty().SetRepresentationToSurface()
            sActor.GetProperty().SetColor(0.2, 0.8, 0.6)
            sActor.GetProperty().SetOpacity(0.9)
            ren.AddActor(sActor)
            self.cadactor = sActor
            self.show_CAD_btn.setText("Hide CAD")
        else:
            if self.cadactor:
                ren.RemoveActor(self.cadactor)
            self.show_CAD_btn.setText("Show CAD")
            self.cadactor = None
        self.showorhideCAD = not self.showorhideCAD
        self.qvtkWidget.Initialize()
        self.qvtkWidget.Render()
        self.qvtkWidget.Start()

    def writemeshfile(self):
        if (self.of_type_rdbtn.isChecked()):
            self.gfmesh.writeOFFile("Cell.txt", "./")
        if (self.netcdf_type_rdbtn.isChecked()):
            self.gfmesh.writeNetCDFFile("Cell.cdl","./")


    def on_xnum_tedit_textChanged(self):
        xnum = re.match(r'\d+',self.xnum_tedit.text())
        if xnum:
            self.x_REcheck_text.setText("")
        else:
            self.x_REcheck_text.setText("Xnum must be positive integer")

    def on_ynum_tedit_textChanged(self):
        ynum = re.match(r'\d+',self.ynum_tedit.text())
        if ynum:
            self.y_REcheck_text.setText("")
        else:
            self.y_REcheck_text.setText("Ynum must be positive integer")

    def on_znum_tedit_textChanged(self):
        znum = re.match(r'\d+',self.znum_tedit.text())
        if znum:
            self.z_REcheck_text.setText("")
        else:
            self.z_REcheck_text.setText("Znum must be positive integer")

    def openFile(self):
        s = QtWidgets.QFileDialog.getOpenFileName(self, "Open file dialog", "./", "CAD files(*.stl)")
        dir = str(s).split()[0][2:-2]
        self.cad_file_dir_text.setText(dir)

    def coordtoindex(self):
        coordx = float(self.coord_check_x_LEdit.text())
        coordy = float(self.coord_check_y_LEdit.text())
        coordz = float(self.coord_check_z_LEdit.text())
        voxindices = [-1,-1,-1]
        voxindices = self.gfmesh.getVoxel3DIndexFromCoord(coordx, coordy, coordz)

        Xindex = voxindices[0]
        Yindex = voxindices[1]
        Zindex = voxindices[2]

        self.index_check_x_LEdit.setText(str(Xindex))
        self.index_check_y_LEdit.setText(str(Yindex))
        self.index_check_z_LEdit.setText(str(Zindex))
        # QtWidgets.QMessageBox.critical(self, "Critical",
        #                      self.tr("Out of range, please check x,y and z"))


if __name__=="__main__":
    import sys
    app=QtWidgets.QApplication(sys.argv)
    gfmeshgui=mywindow()
    gfmeshgui.show()
    sys.exit(app.exec_())
