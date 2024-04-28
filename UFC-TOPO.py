'''
Script para:
Download dos arquivos .tif da base de dados do Topodata/INPE com em uma quadrícula de referência
Extração das elevações de cada ponto de grade do arquivo .tif baixado
Geração das curvas de nível com base na elevação de cada ponto de grade
'''

'''
Autor: Carlos Eduardo Sousa Lima
Departamento de Engenharia Hidráulica e Ambiental/ Universidade Federal do Ceará (DEHA/UFC)
Doutorando em Engenharia Civil (Recursos Hídricos) - DEHA/UFC
Mestre em Engenharia Civil (Recursos Hídricos) - DEHA/UFC
Engenheiro Civil com ênfase em Meio Ambiente - Universidade Estadual Vale do Acaraú (UVA)
'''
'''
Padrão URL 'http://www.dsr.inpe.br/topodata/data/geotiff/'
03S435ZN.zip
LAHLON - LA = Latitude; H = Hemisfério [Norte (N) e Sul (S)]; LON = Longitude
LON = nn_ quando inteio; nn quando
Lat varia em 1° de quadrícula para quadrícula
Long varia em 1,5° de quadrícula para quadrícula
Latitude e Longitude do canto superior esquerdo da quadrícula
'''
import os
import requests
import numpy as np
import geopandas as gpd
import pandas as pd
from glob import glob
import zipfile as zp
import rasterio as rio
from rasterio import Affine
from rasterio.merge import merge
from rasterio.mask import mask
from rasterio.plot import show
from rasterio.warp import reproject, calculate_default_transform, Resampling
from shapely.geometry import box
import pyproj
import matplotlib
import matplotlib.pyplot as plot
import ezdxf as dxf
import PyQt5
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QDialog, QApplication, QFileDialog, QMessageBox

if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    PyQt5.QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    PyQt5.QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

   
class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(720, 600)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setEnabled(True)
        self.centralwidget.setObjectName("centralwidget")
        #############---------------------############-------------------------######### 

        self.dir_shp_line = QtWidgets.QLineEdit(self.centralwidget)
        self.dir_shp_line.setGeometry(QtCore.QRect(10, 30, 611, 20))
        self.dir_shp_line.setText("")
        self.dir_shp_line.setObjectName("dir_shp_line")
        self.dir_shp_button = QtWidgets.QPushButton(self.centralwidget)
        self.dir_shp_button.setGeometry(QtCore.QRect(630, 30, 75, 20))
        self.dir_shp_button.setObjectName("dir_shp_button")
        self.dir_shp_button.clicked.connect(self.buscar_shp)
        self.dir_shp_label = QtWidgets.QLabel(self.centralwidget)
        self.dir_shp_label.setGeometry(QtCore.QRect(10, 10, 161, 16))
        self.dir_shp_label.setObjectName("dir_shp_label")

        #############---------------------############-------------------------######### 

        self.obs_label = QtWidgets.QLabel(self.centralwidget)
        self.obs_label.setGeometry(QtCore.QRect(10, 60, 421, 16))
        self.obs_label.setObjectName("obs_label")

        self.check_coord = QtWidgets.QCheckBox(self.centralwidget)
        self.check_coord.setGeometry(QtCore.QRect(10, 90, 161, 17))
        self.check_coord.setObjectName("check_coord")
        self.check_coord.stateChanged.connect(self.ativar_coord)

        self.long_esq_line = QtWidgets.QLineEdit(self.centralwidget)
        self.long_esq_line.setEnabled(False)
        self.long_esq_line.setGeometry(QtCore.QRect(10, 120, 90, 20))
        self.long_esq_line.setObjectName("long_esq_line")
        self.long_esq_label = QtWidgets.QLabel(self.centralwidget)
        self.long_esq_label.setGeometry(QtCore.QRect(110, 120, 101, 16))
        self.long_esq_label.setObjectName("long_esq_label")
        self.lat_inf_label = QtWidgets.QLabel(self.centralwidget)
        self.lat_inf_label.setGeometry(QtCore.QRect(110, 150, 101, 16))
        self.lat_inf_label.setObjectName("lat_inf_label")
        self.lat_inf_line = QtWidgets.QLineEdit(self.centralwidget)
        self.lat_inf_line.setEnabled(False)
        self.lat_inf_line.setGeometry(QtCore.QRect(10, 150, 90, 20))
        self.lat_inf_line.setObjectName("lat_inf_line")
        self.long_dir_label = QtWidgets.QLabel(self.centralwidget)
        self.long_dir_label.setGeometry(QtCore.QRect(320, 120, 101, 16))
        self.long_dir_label.setObjectName("long_dir_label")
        self.long_line_dir = QtWidgets.QLineEdit(self.centralwidget)
        self.long_line_dir.setEnabled(False)
        self.long_line_dir.setGeometry(QtCore.QRect(220, 120, 90, 20))
        self.long_line_dir.setText("")
        self.long_line_dir.setObjectName("long_line_dir")
        self.lat_sup_line = QtWidgets.QLineEdit(self.centralwidget)
        self.lat_sup_line.setEnabled(False)
        self.lat_sup_line.setGeometry(QtCore.QRect(220, 150, 90, 20))
        self.lat_sup_line.setText("")
        self.lat_sup_line.setObjectName("lat_sup_line")
        self.lat_sup_label = QtWidgets.QLabel(self.centralwidget)
        self.lat_sup_label.setGeometry(QtCore.QRect(320, 150, 101, 16))
        self.lat_sup_label.setObjectName("lat_sup_label")
        self.proj_label = QtWidgets.QLabel(self.centralwidget)
        self.proj_label.setGeometry(QtCore.QRect(280, 90, 361, 16))
        self.proj_label.setObjectName("proj_label")
        self.proj_line = QtWidgets.QLineEdit(self.centralwidget)
        self.proj_line.setEnabled(False)
        self.proj_line.setGeometry(QtCore.QRect(180, 90, 90, 20))
        self.proj_line.setText("")
        self.proj_line.setObjectName("proj_line")

        #############---------------------############-------------------------#########

        self.Scale_Horizontal_Line = QtWidgets.QLineEdit(self.centralwidget)
        self.Scale_Horizontal_Line.setEnabled(True)
        self.Scale_Horizontal_Line.setGeometry(QtCore.QRect(10, 200, 90, 20))
        self.Scale_Horizontal_Line.setObjectName("Scale_Horizontal_Line")
        self.Scale_Horizontal_Label = QtWidgets.QLabel(self.centralwidget)
        self.Scale_Horizontal_Label.setGeometry(QtCore.QRect(10, 180, 101, 16))
        self.Scale_Horizontal_Label.setObjectName("Scale_Horizontal_Label")
        self.Scale_Vertical_Line = QtWidgets.QLineEdit(self.centralwidget)
        self.Scale_Vertical_Line.setEnabled(True)
        self.Scale_Vertical_Line.setGeometry(QtCore.QRect(150, 200, 90, 20))
        self.Scale_Vertical_Line.setObjectName("Scale_Vertical_Line")
        self.Scale_Vertical_Label = QtWidgets.QLabel(self.centralwidget)
        self.Scale_Vertical_Label.setGeometry(QtCore.QRect(150, 180, 101, 16))
        self.Scale_Vertical_Label.setObjectName("Scale_Vertical_Label")
        self.Scale_Horizontal_Line.setText("1.2")
        self.Scale_Vertical_Line.setText("1.2")

        #############---------------------############-------------------------######### 
        
        self.check_getraster = QtWidgets.QCheckBox(self.centralwidget)
        self.check_getraster.setGeometry(QtCore.QRect(10, 230, 291, 17))
        self.check_getraster.setObjectName("check_getraster")
        self.check_getraster.stateChanged.connect(self.ativar_getraster)

        self.dir_getraster_label = QtWidgets.QLabel(self.centralwidget)
        self.dir_getraster_label.setGeometry(QtCore.QRect(10, 250, 300, 20))
        self.dir_getraster_label.setObjectName("dir_getraster_label")
        self.dir_getraster_line = QtWidgets.QLineEdit(self.centralwidget)
        self.dir_getraster_line.setGeometry(QtCore.QRect(10, 270, 611, 20))
        self.dir_getraster_line.setText("")
        self.dir_getraster_line.setObjectName("dir_getraster_line")
        self.dir_getraster_button = QtWidgets.QPushButton(self.centralwidget)
        self.dir_getraster_button.setGeometry(QtCore.QRect(630, 270, 75, 20))
        self.dir_getraster_button.setObjectName("dir_getraster_button")
        self.dir_getraster_button.clicked.connect(self.buscar_input_raster)

        #############---------------------############-------------------------######### 
        
        self.dir_saveraster_line = QtWidgets.QLineEdit(self.centralwidget)
        self.dir_saveraster_line.setGeometry(QtCore.QRect(10, 320, 611, 20))
        self.dir_saveraster_line.setText("")
        self.dir_saveraster_line.setObjectName("dir_saveraster_line")
        self.dir_saveraster_button = QtWidgets.QPushButton(self.centralwidget)
        self.dir_saveraster_button.setGeometry(QtCore.QRect(630, 320, 75, 20))
        self.dir_saveraster_button.setObjectName("dir_saveraster_button")
        self.dir_saveraster_button.clicked.connect(self.buscar_output_raster)
        self.dir_save_raster_label = QtWidgets.QLabel(self.centralwidget)
        self.dir_save_raster_label.setGeometry(QtCore.QRect(10, 300, 201, 16))
        self.dir_save_raster_label.setObjectName("dir_save_raster_label")


        #############---------------------############-------------------------######### 

        self.check_xyz = QtWidgets.QCheckBox(self.centralwidget)
        self.check_xyz.setGeometry(QtCore.QRect(10, 350, 211, 17))
        self.check_xyz.setObjectName("check_xyz")
        self.check_xyz.stateChanged.connect(self.ativar_xyz)

        self.dir_xyz_label = QtWidgets.QLabel(self.centralwidget)
        self.dir_xyz_label.setGeometry(QtCore.QRect(10, 370, 221, 16))
        self.dir_xyz_label.setObjectName("dir_xyz_label")
        self.dir_xyz_line = QtWidgets.QLineEdit(self.centralwidget)
        self.dir_xyz_line.setEnabled(False)
        self.dir_xyz_line.setGeometry(QtCore.QRect(10, 390, 611, 20))
        self.dir_xyz_line.setText("")
        self.dir_xyz_line.setObjectName("dir_xyz_line")
        self.dir_xyz_button = QtWidgets.QPushButton(self.centralwidget)
        self.dir_xyz_button.setGeometry(QtCore.QRect(630, 390, 75, 20))
        self.dir_xyz_button.setObjectName("dir_xyz_button")
        self.dir_xyz_button.clicked.connect(self.buscar_dir_xyz)

        #############---------------------############-------------------------#########
        
        self.check_curva = QtWidgets.QCheckBox(self.centralwidget)
        self.check_curva.setGeometry(QtCore.QRect(10, 420, 231, 17))
        self.check_curva.setChecked(False)
        self.check_curva.setObjectName("check_curva")
        self.check_curva.stateChanged.connect(self.ativar_dxf)

        self.dir_dxf_label = QtWidgets.QLabel(self.centralwidget)
        self.dir_dxf_label.setGeometry(QtCore.QRect(10, 480, 251, 16))
        self.dir_dxf_label.setObjectName("dir_dxf_label")
        self.dir_dxf_line = QtWidgets.QLineEdit(self.centralwidget)
        self.dir_dxf_line.setEnabled(False)
        self.dir_dxf_line.setGeometry(QtCore.QRect(10, 500, 611, 20))
        self.dir_dxf_line.setText("")
        self.dir_dxf_line.setObjectName("dir_dxf_line")
        self.dir_dxf_button = QtWidgets.QPushButton(self.centralwidget)
        self.dir_dxf_button.setGeometry(QtCore.QRect(630, 500, 75, 20))
        self.dir_dxf_button.setObjectName("dir_dxf_button")
        self.dir_dxf_button.clicked.connect(self.buscar_dir_dxf)

        self.curva_intervalo_label = QtWidgets.QLabel(self.centralwidget)
        self.curva_intervalo_label.setGeometry(QtCore.QRect(110, 450, 331, 16))
        self.curva_intervalo_label.setObjectName("curva_intervalo_label")
        self.curva_intervalo_line = QtWidgets.QLineEdit(self.centralwidget)
        self.curva_intervalo_line.setEnabled(False)
        self.curva_intervalo_line.setGeometry(QtCore.QRect(10, 450, 90, 20))
        self.curva_intervalo_line.setObjectName("curva_intervalo_line")
        self.curva_intervalo_line.setText("1")

        #############---------------------############-------------------------#########
        self.exec_button = QtWidgets.QPushButton(self.centralwidget)
        self.exec_button.setGeometry(QtCore.QRect(10, 530, 141, 41))
        font = QtGui.QFont()
        font.setBold(True)
        font.setUnderline(False)
        font.setWeight(75)
        self.exec_button.setFont(font)
        self.exec_button.setObjectName("exec_button")
        self.exec_button.clicked.connect(self.executar)


        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.long_esq_label.setText(_translate("MainWindow", "Longitude Esquerda"))
        self.lat_inf_label.setText(_translate("MainWindow", "Latitude Inferior"))
        self.long_dir_label.setText(_translate("MainWindow", "Longitude Direita"))
        self.lat_sup_label.setText(_translate("MainWindow", "Latitude Superior"))
        self.proj_label.setText(_translate("MainWindow", "Projeção (Fomarto EPSG) // Ex: EPSG:31984 = SIRGAS 2000 UTM 24S"))
        self.dir_shp_button.setText(_translate("MainWindow", "Buscar"))
        self.dir_shp_label.setText(_translate("MainWindow", "Shapefile (Projeção em UTM):"))
        self.dir_saveraster_button.setText(_translate("MainWindow", "Buscar"))
        self.dir_save_raster_label.setText(_translate("MainWindow", "Diretório de saída dos Rasters:"))
        self.dir_xyz_label.setText(_translate("MainWindow", "Diretório de saída do arquivo PontoCota.txt:"))
        self.dir_xyz_button.setText(_translate("MainWindow", "Buscar"))
        self.dir_dxf_label.setText(_translate("MainWindow", "Diretório de saída do arquivo Curva_de_Nível.dxf:"))
        self.dir_dxf_button.setText(_translate("MainWindow", "Buscar"))
        self.curva_intervalo_label.setText(_translate("MainWindow", "Variação de elevação entre as curvas de nível (Default = 1 metro)"))
        self.exec_button.setText(_translate("MainWindow", "Executar UFC-TOPO"))
        self.obs_label.setText(_translate("MainWindow", "Caso queira, ou não haja shapefile, informe as coordenadas da quadrícula desejada:"))
        self.check_coord.setText(_translate("MainWindow", "Informar Coordenadas (UTM):"))
        self.check_xyz.setText(_translate("MainWindow", "Gerar arquivo PontoCota.txt (X, Y, Z)"))
        self.check_curva.setText(_translate("MainWindow", "Gerar Curvas de Nível em um arquivo .dxf:"))
        self.Scale_Horizontal_Label.setText(_translate("MainWindow", "Scale Horizontal:"))
        self.Scale_Vertical_Label.setText(_translate("MainWindow", "Scale Vertical:"))
        self.check_getraster.setText(_translate("MainWindow", "Download das informações Topográficas do TOPODATA"))
        self.dir_getraster_label.setText(_translate("MainWindow", "Diretório do arquivo Raster com informações altimétricas:"))
        self.dir_getraster_button.setText(_translate("MainWindow", "Buscar"))

    def ativar_coord(self):
        if self.check_coord.isChecked():
            self.proj_line.setEnabled(True)
            self.lat_inf_line.setEnabled(True)
            self.lat_sup_line.setEnabled(True)
            self.long_esq_line.setEnabled(True)
            self.long_line_dir.setEnabled(True)
            self.dir_shp_line.setEnabled(False)
            self.dir_shp_button.setEnabled(False)
            self.dir_shp_line.setText("")
        else:
            self.proj_line.setEnabled(False)
            self.lat_inf_line.setEnabled(False)
            self.lat_sup_line.setEnabled(False)
            self.long_esq_line.setEnabled(False)
            self.long_line_dir.setEnabled(False)
            self.dir_shp_line.setEnabled(True)
            self.dir_shp_button.setEnabled(True)

    def ativar_getraster(self):
        if self.check_getraster.isChecked():
            self.dir_getraster_label.setEnabled(False)
            self.dir_getraster_button.setEnabled(False)
            self.dir_getraster_line.setText("")
            self.dir_getraster_line.setEnabled(False)
        else:
            self.dir_getraster_label.setEnabled(True)
            self.dir_getraster_button.setEnabled(True)
            self.dir_getraster_line.setEnabled(True)

    def ativar_xyz(self):
        if self.check_xyz.isChecked():
            self.dir_xyz_line.setEnabled(True)
            self.dir_xyz_button.setEnabled(True)
        else:
            self.dir_xyz_line.setEnabled(False)
            self.dir_xyz_button.setEnabled(False)


    def ativar_dxf(self):
        if self.check_curva.isChecked():
            self.dir_dxf_line.setEnabled(True)
            self.curva_intervalo_line.setEnabled(True)
            self.dir_dxf_button.setEnabled(True)

        else:
            self.dir_dxf_line.setEnabled(False)
            self.dir_dxf_line.setEnabled(False)
            self.dir_dxf_button.setEnabled(False)
            self.curva_intervalo_line.setEnabled(False)


    def buscar_shp(self):
        file_name = QFileDialog.getOpenFileName(None, "Buscar Shape", os.getcwd(), "Shapefile (*.shp)")
        self.dir_shp_line.setText(file_name[0])

    def buscar_input_raster(self):
        dir_name = QFileDialog.getExistingDirectory(None, "Diretórios dos Rasters de entrada", os.getcwd())
        self.dir_getraster_line.setText(dir_name)

    def buscar_output_raster(self):
        dir_name = QFileDialog.getExistingDirectory(None, "Diretório de saída dos Rasters", os.getcwd())
        self.dir_saveraster_line.setText(dir_name)

    def buscar_dir_xyz(self):
        dir_name = QFileDialog.getExistingDirectory(None, "Diretório PontoCota (x, y, z)", os.getcwd())
        self.dir_xyz_line.setText(dir_name+"/PontoCota.txt")

    def buscar_dir_dxf(self):
        dir_name = QFileDialog.getExistingDirectory(None, "Diretório .dxf", os.getcwd())
        self.dir_dxf_line.setText(dir_name+"/Curva_de_Nível.dxf")


    def executar(self):
        if self.check_coord.isChecked():
            read_shp = False
            coords = [float(self.long_esq_line.text()),
                      float(self.lat_inf_line.text()),
                      float(self.long_line_dir.text()),
                      float(self.lat_sup_line.text())]
            proj = "EPSG:"+self.proj_line.text()
            dir_shp = None
        else:
            read_shp = True
            dir_shp = self.dir_shp_line.text()
            coords = None
            proj = None
        #dir_save_raster onde será criada a pasta SRTM e baixado os rasters
        #path_raster dir_save_raster + SRTM

        scale_factor = [float(self.Scale_Horizontal_Line.text()), float(self.Scale_Vertical_Line.text())]
        dir_save_raster = self.dir_saveraster_line.text()

        if self.check_getraster.isChecked():
            path_raster = Download_Topodata(dir_shp, read_shp, coords, proj, dir_save_raster)
            Descompactar_zips(path_raster)
            dir_output_merge = path_raster + "/merge_WGS.tif"
            merge_raster(path_raster, dir_output_merge)
            dir_output_utm = path_raster + "/merge_UTM.tif"
            project_rater(dir_output_merge, dir_output_utm, dir_shp, read_shp, proj)
            dir_output_mask = path_raster + "/mask_UTM.tif"
            mask_raster(dir_output_utm, dir_output_mask, dir_shp, read_shp, coords, proj, scale_factor)
            if self.check_xyz.isChecked():
                dir_output_csv = self.dir_xyz_line.text()
                extract_xyz(dir_output_utm, dir_shp, dir_output_csv, read_shp, coords, proj, scale_factor)
                if self.check_curva.isChecked():
                    delta_z = float(self.curva_intervalo_line.text())
                    dir_output_dxf = self.dir_dxf_line.text()
                    extract_curves(dir_output_csv, dir_output_dxf, delta_z)

            exlcuir_raster(path_raster)
        else:

            dir_get_raster = self.dir_getraster_line.text()

            dir_output_merge = dir_save_raster + '/merge_WGS.tif'
            dir_output_utm = dir_save_raster + '/merge_UTM.tif'
            dir_output_mask = dir_save_raster + '/mask_UTM.tif'

            if os.path.exists(dir_output_utm):
                os.remove(dir_output_utm)
            if os.path.exists(dir_output_merge):
                os.remove(dir_output_merge)
            if os.path.exists(dir_output_mask):
                os.remove(dir_output_mask)
            
            merge_raster(dir_get_raster, dir_output_merge)
            project_rater(dir_output_merge, dir_output_utm, dir_shp, read_shp, proj)
            mask_raster(dir_output_utm, dir_output_mask, dir_shp, read_shp, coords, proj, scale_factor)

            if self.check_xyz.isChecked():
                dir_output_csv = self.dir_xyz_line.text()
                extract_xyz(dir_output_utm, dir_shp, dir_output_csv, read_shp, coords, proj, scale_factor)

                if self.check_curva.isChecked():
                    delta_z = float(self.curva_intervalo_line.text())
                    dir_output_dxf = self.dir_dxf_line.text()
                    extract_curves(dir_output_csv, dir_output_dxf, delta_z)

            if os.path.exists(dir_output_utm):
                os.remove(dir_output_utm)
            if os.path.exists(dir_output_merge):
                os.remove(dir_output_merge)

        msg = QMessageBox()
        msg.setWindowTitle("UFC-TOPO")
        msg.setText("Finalizado!")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec()


def Download_Topodata(dir_shp, read_shp, coords, proj, dir_save_raster):  #Como entrar, o diretório do Shape
    output_path = dir_save_raster + "/SRTM"  #Diretório de saída
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    else:
        list_remove = glob(output_path + "/*")
        for path_remove in list_remove:
            os.remove(path_remove)
        os.rmdir(output_path)
        os.makedirs(output_path)

    if read_shp == True:
        print("Lendo shapefile..........")
        shp = gpd.read_file(dir_shp)
        if shp.crs != "EPSG:4326":
            shp = shp.to_crs("EPSG:4326")
        extents = shp.geometry.total_bounds
        long_esq, lat_inf, long_dir, lat_sup = extents
    else:
        print ("Lendo Coordenadas..........")
        #coords = [long_esq, lat_inf, long_dir, lat_sup]
        transformer = pyproj.Transformer.from_crs(proj, "EPSG:4326", always_xy=True)
        long_esq, lat_inf = transformer.transform(coords[0], coords[1])
        long_dir, lat_sup = transformer.transform(coords[2], coords[3])
    #Entrada

    H = "S"  # If H = S, lat_sup e lat_inf é negativa... Olhar isso depois.

    # Formação da quadrícula
    sup_esq = (lat_sup, long_esq)
    inf_esq = (lat_inf, long_esq)
    sup_dir = (lat_sup, long_dir)
    inf_dir = (lat_inf, long_dir)

    # Seleção das quardículas
    #Na positiva é a parte interira +1 na negativa só a parte inteira
    lat_ini = int(lat_sup)
    lat_end = int(lat_inf)
    long_ini = ((int(long_dir / 1.5) - 1) * 1.5)
    long_end = ((int(long_esq / 1.5) - 1) * 1.5)

    # Nomeclatura das quadrículas
    #if para saber ajeitar os sinais da iteração
    lat_array = np.arange(abs(lat_ini), abs(lat_end) + 1, 1)
    long_array = np.arange(abs(long_ini), abs(long_end) + 1, 1.5)
    qua_list = []
    for i in lat_array:
        for j in long_array:
            qua_list.append((i, j))
    qua_array = np.array(qua_list)

    # Criando a URL

    # url = "http://www.dsr.inpe.br/topodata/data/geotiff/03S435ZN.zip"
    url1 = "http://www.dsr.inpe.br/topodata/data/geotiff/"
    url2_list = []
    for k in qua_array:
        #if para saber se a latitude é norte ou sul
        if ((k[0] < 10) and (str(k[1])[3] != "0")):
            url2_aux = ("0" + str(int(k[0])) + H + str((k[1])).replace(".", '') + "ZN.zip")
            url2_list.append(url2_aux)
            print("Link Gerado..........")

        if ((k[0] < 10) and (str(k[1])[3] == "0")):
            url2_aux = ("0" + str(int(k[0])) + H + str((k[1])).split(".")[0] + "_ZN.zip")
            url2_list.append(url2_aux)
            print("Link Gerado..........")

        if ((k[0] >= 10) and (str(k[1])[3] != "0")):
            url2_aux = (str(int(k[0])) + H + str((k[1])).replace(".", '') + "ZN.zip")
            url2_list.append(url2_aux)
            print("Link Gerado..........")

        if ((k[0] >= 10) and (str(k[1])[3] == "0")):
            url2_aux = (str(int(k[0])) + H + str((k[1])).split(".")[0] + "_ZN.zip")
            url2_list.append(url2_aux)
            print("Link Gerado..........")

    for url2 in url2_list:
        print("Baixando Quadrícula..........")
        url = url1 + url2
        server_INPE = requests.get(url, stream=True)
        arq = (output_path + "/" + url2)
        with open(arq, 'wb') as arq:
            arq.write(server_INPE.content)
        print("Quadrículas Baixadas..........")
    return output_path


def Descompactar_zips (path):
    zip_list = glob(path + "/" + "*zip")
    for zip in zip_list:
        zip_file = zp.ZipFile(zip)
        zip_file.extractall(path=path)
        zip_file.close()
        os.remove(zip)
    print("Arquivos Descompactados..........")
    return


def merge_raster(dir_raster, dir_output_merge):
    raster_list = glob(dir_raster + "/" + "*.tif")
    mosaic_raster = []
    for path_raster in raster_list: #mudar o nome dessa variável, dir_save_raster ou path_raster
        raster = rio.open(path_raster)
        mosaic_raster.append(raster)
    raster_merge, out_merge = merge(mosaic_raster)
    out_meta = raster.meta.copy()
    out_meta.update({"driver": "GTiff",
                        "height": raster_merge.shape[1],
                        "width": raster_merge.shape[2],
                        "transform": out_merge,
                        "crs": pyproj.CRS.to_wkt(pyproj.crs.CRS("EPSG:4326"))})
    with rio.open(dir_output_merge, "w", **out_meta) as dst:
        dst.write(raster_merge)
    print("Merge das quadrículas concluído..........")

    return


def project_rater(dir_save_raster, dir_output_utm, dir_shp, read_shp, proj):
    if read_shp == True:
        shp = gpd.read_file(dir_shp)
        output_crs = pyproj.CRS.to_wkt(shp.crs)
    else:
        output_crs = pyproj.CRS.to_wkt(pyproj.crs.CRS(proj))
    raster = rio.open(dir_save_raster)
    out_meta = raster.meta.copy()
    transform, width, height = calculate_default_transform(
        raster.crs, output_crs, raster.width, raster.height, *raster.bounds)
    out_meta.update({
        'crs': output_crs,
        'transform': transform,
        'width': width,
        'heigth': height})
    with rio.open(dir_output_utm, 'w', **out_meta) as dst:
        reproject(
            source=rio.band(raster, 1),
            destination=rio.band(dst, 1),
            src_transform=raster.transform,
            src_crs=raster.crs,
            dst_transform=transform,
            dst_crs=output_crs,
            resampling=Resampling.nearest)
    print("Projeção para UTM realizada..........")

    return


def mask_raster(dir_save_raster, dir_output_mask, dir_shp, read_shp, coords, proj, scale):
    #"crs": prj.crs.CRS.to_proj4(prj.CRS(raster.crs)) retorna a projeção como:
    #"+proj=longlat +datum=WGS84 +no_defs +type=crs". só funcionou esse formato no out_meta.upadate
    raster = rio.open(dir_save_raster)
    if read_shp == True:
        shp = gpd.read_file(dir_shp)
        shp_clip = shp.envelope.scale(scale[0], scale[1]) #scale[scale_x, scale_y]
        raster_mask, out_mask = mask(raster, shp_clip.geometry, crop=True, filled=True, nodata=np.nan)
        out_meta = raster.meta.copy()
        out_meta.update({"driver": "GTiff",
                            "height": raster_mask.shape[1],
                            "width": raster_mask.shape[2],
                            "transform": out_mask,
                            "crs": pyproj.CRS.to_wkt(shp.crs),
                            "nodata": np.nan})
    else:
        bbox = box(coords[0], coords[1], coords[2], coords[3])
        gdf_bbox = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs=proj)
        bbox_clip = gdf_bbox.envelope.scale(scale[0], scale[1]) #scale[scale_x, scale_y]
        raster_mask, out_mask = mask(raster, bbox_clip.geometry, crop=True, filled=True, nodata=np.nan)
        out_meta = raster.meta.copy()
        out_meta.update({"driver": "GTiff",
                            "height": raster_mask.shape[1],
                            "width": raster_mask.shape[2],
                            "transform": out_mask,
                            "crs": pyproj.CRS.to_wkt(pyproj.crs.CRS(proj)),
                            "nodata": np.nan})
    with rio.open(dir_output_mask, "w", **out_meta) as dst:
        dst.write(raster_mask)
    print("Clip concluído..........")
    return


def extract_xyz(dir_save_raster, dir_shp, dir_output_csv, read_shp, coords, proj, scale):
    #Refaço o Mask aqui, pois preciso de informações que preciso dentro da variáveil out_mask
    #Não consegui fazer puxando do raster Mask_UTM
    #Outro ponto, é que preciso de funções diferentes, pois caso não queira o arquivo PontoCota, posso obter só o raster recortado
    dir_output_csv = dir_output_csv
    raster = rio.open(dir_save_raster)
    if read_shp == True:
        shp = gpd.read_file(dir_shp) #Lendo Shapefile
        shp_clip = shp.envelope.scale(scale[0], scale[1]) #Fator de escala → scale[scale_x, scale_y]
        raster_mask, out_mask = mask(raster, shp_clip.geometry, crop=True, filled=True, nodata=-999)
    else:
        bbox = box(coords[0], coords[1], coords[2], coords[3])
        gdf_bbox = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs=proj)
        bbox_clip = gdf_bbox.envelope.scale(scale[0], scale[1]) #scale[scale_x, scale_y]
        raster_mask, out_mask = mask(raster, bbox_clip.geometry, crop=True, filled=True, nodata=-999)
    data_elev = raster_mask[0]
    row, col = np.where(data_elev != -999)  #Linhas e colunas onde a elevação não é NaN
    elev = np.extract(data_elev != -999, data_elev)  #Elevações diferentes de NaN
    aff_center = out_mask * Affine.translation(0.5, 0.5)  #Referencia o centro de cada pixel, Ler transformação Affine
    long, lat = aff_center * (col, row)
    df_gpd = gpd.GeoDataFrame({'Long': long, 'Lat': lat, 'Elev': elev})
    df_gpd.to_csv(dir_output_csv, index=False, header=False)
    print("Pontos Extraídos..........")
    return


def extract_curves(dir_output_csv,  dir_output_dxf, delta_z):
    dir_output_dxf = dir_output_dxf
    xyz_table = pd.read_csv(dir_output_csv, header=None, names=["x", "y", "z"])
    z_matrix = xyz_table.pivot_table(index="x", columns="y", values="z").T.values
    # values torna as linhas e colunas em n° [0, 1, 2, 3], sem ele, as linhas seriam as coordenadas.
    # Não referenciaria com o grid // T Transpõe a matriz
    x_unique = np.sort(xyz_table.x.unique())
    y_unique = np.sort(xyz_table.y.unique())
    ## Unique exclui os valores repetidos, mostra apenas os valores contidos no array
    X, Y = np.meshgrid(x_unique, y_unique)
    #Cria o grid, retornando as coordenadas X e Y de cada nó
    z_min = round(z_matrix.min(), 0)
    z_max = round(z_matrix.max(), 0)
    lvl = np.arange(z_min, z_max + 1, delta_z) #Espaçamento entre as curvas de nível
    plot.ioff()
    contours = plot.contour(X, Y, z_matrix, levels=lvl)
    doc = dxf.new(dxfversion='AC1032', setup=True)
    msp = doc.modelspace()
    doc.layers.new(name="Curva_de_Nivel", dxfattribs={'linetype': "CONTINUOUS", 'color': 7})
    for i in range(0, len(contours.allsegs)):
        segments = contours.allsegs[i]
        label = contours.levels[i]
        if (len(segments) == 0):
            print("Segmento Zerado!.....")
        else:
            for segment in segments:
                fit_points = []
                if (len(segment) > 1):
                    for vertice in segment:
                        x = vertice[0]
                        y = vertice[1]
                        fit_points.append([x, y, label])
                    spline = msp.add_spline(fit_points, dxfattribs={'layer': 'Curva_de_Nivel'})
                else:
                    print("Segmento com um ponto!....")
    print("Salvando Arquivo!")
    doc.saveas(dir_output_dxf)

def exlcuir_raster(dir_save_raster):
    list_raster = glob(dir_save_raster + "/" + "*.tif")
    for i in list_raster:
        if i != "{}\mask_UTM.tif".format(dir_save_raster):#Deixar apenas o mask_UTM
            os.remove(i)
if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
