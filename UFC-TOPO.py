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
from glob import glob
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

def Download_Topodata(dir_shp, read_shp, coords, proj, dir_raster):  #Como entrar, o diretório do Shape
    output_path = dir_raster + '/SRTM'  #Diretório de saída
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    else:
        list_remove = glob(output_path + '/*')
        for path_remove in list_remove:
            os.remove(path_remove)
        os.rmdir(output_path)
        os.makedirs(output_path)

    if read_shp == True:
        print('Lendo shapefile..........')
        shp = gpd.read_file(dir_shp)
        if shp.crs != 'EPSG:4326':
            shp = shp.to_crs('EPSG:4326')
        extents = shp.geometry.total_bounds
        long_esq, lat_inf, long_dir, lat_sup = extents
    else:
        print ('Lendo Coordenadas..........')
        #coords = [long_esq, lat_inf, long_dir, lat_sup]
        transformer = pyproj.Transformer.from_crs(proj, 'EPSG:4326', always_xy=True)
        long_esq, lat_inf = transformer.transform(coords[0], coords[1])
        long_dir, lat_sup = transformer.transform(coords[2], coords[3])
    #Entrada

    H = 'S'  # If H = S, lat_sup e lat_inf é negativa... Olhar isso depois.

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

    # url = 'http://www.dsr.inpe.br/topodata/data/geotiff/03S435ZN.zip'
    url1 = 'http://www.dsr.inpe.br/topodata/data/geotiff/'
    url2_list = []
    for k in qua_array:
        #if para saber se a latitude é norte ou sul
        if ((k[0] < 10) and (str(k[1])[3] != '0')):
            url2_aux = ('0' + str(int(k[0])) + H + str((k[1])).replace('.', '') + 'ZN.zip')
            url2_list.append(url2_aux)
            print('Link Gerado..........')

        if ((k[0] < 10) and (str(k[1])[3] == '0')):
            url2_aux = ('0' + str(int(k[0])) + H + str((k[1])).split('.')[0] + '_ZN.zip')
            url2_list.append(url2_aux)
            print('Link Gerado..........')

        if ((k[0] >= 10) and (str(k[1])[3] != '0')):
            url2_aux = (str(int(k[0])) + H + str((k[1])).replace('.', '') + 'ZN.zip')
            url2_list.append(url2_aux)
            print('Link Gerado..........')

        if ((k[0] >= 10) and (str(k[1])[3] == '0')):
            url2_aux = (str(int(k[0])) + H + str((k[1])).split('.')[0] + '_ZN.zip')
            url2_list.append(url2_aux)
            print('Link Gerado..........')

    for url2 in url2_list:
        print('Baixando Quadrícula..........')
        url = url1 + url2
        server_INPE = requests.get(url, stream=True)
        arq = (output_path + '/' + url2)
        with open(arq, 'wb') as arq:
            arq.write(server_INPE.content)
        print('Quadrículas Baixadas..........')
    return output_path


def Descompactar_zips (path):
    zip_list = glob(path + '/' + '*zip')
    for zip in zip_list:
        zip_file = zp.ZipFile(zip)
        zip_file.extractall(path=path)
        zip_file.close()
        os.remove(zip)
    print('Arquivos Descompactados..........')
    return


def merge_raster(dir_raster, dir_output_merge):

    raster_list = glob(dir_raster + '/' + '*.tif')

    mosaic_raster = []
    for output_path in raster_list: #mudar o nome dessa variável, dir_raster ou output_path
        raster = rio.open(output_path)
        mosaic_raster.append(raster)
    raster_merge, out_merge = merge(mosaic_raster)
    out_meta = raster.meta.copy()
    out_meta.update({'driver': 'GTiff',
                        'height': raster_merge.shape[1],
                        'width': raster_merge.shape[2],
                        'transform': out_merge,
                        'crs': pyproj.CRS.to_wkt(pyproj.crs.CRS('EPSG:4326'))})
    with rio.open(dir_output_merge, 'w', **out_meta) as dst:
        dst.write(raster_merge)
    print('Merge das quadrículas concluído..........')

    return dir_output_merge


def project_rater(dir_raster, dir_output_utm, dir_shp, read_shp, proj):
    if read_shp == True:
        shp = gpd.read_file(dir_shp)
        output_crs = pyproj.CRS.to_wkt(shp.crs)
    else:
        output_crs = pyproj.CRS.to_wkt(pyproj.crs.CRS(proj))
    raster = rio.open(dir_raster)
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
    print('Projeção para UTM realizada..........')

    return


def mask_raster(dir_raster, dir_output_mask, dir_shp, read_shp, coords, proj, scale, ):
    #'crs': prj.crs.CRS.to_proj4(prj.CRS(raster.crs)) retorna a projeção como:
    #'+proj=longlat +datum=WGS84 +no_defs +type=crs'. só funcionou esse formato no out_meta.upadate
    raster = rio.open(dir_raster)
    if read_shp == True:
        shp = gpd.read_file(dir_shp)
        shp_clip = shp.envelope.scale(scale[0], scale[1]) #scale[scale_x, scale_y]
        raster_mask, out_mask = mask(raster, shp_clip.geometry, crop=True, filled=True, nodata=np.nan)
        out_meta = raster.meta.copy()
        out_meta.update({'driver': 'GTiff',
                            'height': raster_mask.shape[1],
                            'width': raster_mask.shape[2],
                            'transform': out_mask,
                            'crs': pyproj.CRS.to_wkt(shp.crs),
                            'nodata': np.nan})
    else:
        bbox = box(coords[0], coords[1], coords[2], coords[3])
        gdf_bbox = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs=proj)
        bbox_clip = gdf_bbox.envelope.scale(scale[0], scale[1]) #scale[scale_x, scale_y]
        raster_mask, out_mask = mask(raster, bbox_clip.geometry, crop=True, filled=True, nodata=np.nan)
        out_meta = raster.meta.copy()
        out_meta.update({'driver': 'GTiff',
                            'height': raster_mask.shape[1],
                            'width': raster_mask.shape[2],
                            'transform': out_mask,
                            'crs': pyproj.CRS.to_wkt(pyproj.crs.CRS(proj)),
                            'nodata': np.nan})
    with rio.open(dir_output_mask, 'w', **out_meta) as dst:
        dst.write(raster_mask)
    print('Clip concluído..........')
    return


def extract_xyz(dir_raster, dir_shp, dir_output_csv, read_shp, coords, proj, scale):
    #Refaço o Mask aqui, pois preciso de informações que preciso dentro da variáveil out_mask
    #Não consegui fazer puxando do raster Mask_UTM
    #Outro ponto, é que preciso de funções diferentes, pois caso não queira o arquivo PontoCota, posso obter só o raster recortado
    dir_output_csv = dir_output_csv
    raster = rio.open(dir_raster)
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
    print('Pontos Extraídos..........')
    return


def extract_curves(dir_output_csv,  dir_output_dxf, delta_z):
    dir_output_dxf = dir_output_dxf
    xyz_table = pd.read_csv(dir_output_csv, header=None, names=['x', 'y', 'z'])
    z_matrix = xyz_table.pivot_table(index='x', columns='y', values='z').T.values
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
    doc.layers.new(name='Curva_de_Nivel', dxfattribs={'linetype': 'CONTINUOUS', 'color': 7})
    for i in range(0, len(contours.allsegs)):
        segments = contours.allsegs[i]
        label = contours.levels[i]
        if (len(segments) == 0):
            print('Segmento Zerado!.....')
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
                    print('Segmento com um ponto!....')
    print('Salvando Arquivo!')
    doc.saveas(dir_output_dxf)
    return

def exlcuir_raster(dir_raster):
    list_raster = glob(dir_raster + '/' + '*.tif')
    for i in list_raster:
        if i != '{}\mask_UTM.tif'.format(dir_raster):#Deixar apenas o mask_UTM
            os.remove(i)
    return

if __name__ == '__main__':
    
    print('Como informar a região de estudo?')

    ref_mod = input('1 - Shapefile\n2 - Coordenadas\nDigite uma Opção: ')
    proj = 'EPSG:32723'
    #Mudar proj caso for passado coordenadas
    #Se passar shapefile, ele lê a projeção do shapefile e atribui a variável proj

    if ref_mod == '1':
        read_shp = True
        dir_shp = glob('shp/*.shp')[0]
        coords = list()

    if ref_mod == '2':
        read_shp = False
        dir_shp = ''
        coords_df = pd.read_csv('coords.csv', sep = ',')
        #coords = [long_esq, lat_inf, long_dir, lat_sup]
        coords = [coords_df['lon_esq'][0], coords_df['lat_inf'][0], coords_df['lon_dir'][0], coords_df['lat_sup'][0]]


    print('\nModo de aquisição das informações topográficas?')
    input_mode = input('1 - Arquivo Existente\n2 - Download Topodata (Temporariamente Inativa)\nDigite uma Opção: ')

    if input_mode == '1':
        
        path_raster = 'raster'
        dir_raster =  glob(path_raster + '/*.tif')[0]
        
        dir_output_merge = path_raster + '/merge_WGS.tif'
        dir_output_utm = path_raster + '/merge_UTM.tif'
        dir_output_mask = path_raster + '/mask_UTM.tif'

        if os.path.exists(dir_output_utm):
            os.remove(dir_output_utm)
        if os.path.exists(dir_output_merge):
            os.remove(dir_output_merge)
        if os.path.exists(dir_output_mask):
            os.remove(dir_output_mask)

        merge_raster(path_raster, dir_output_merge)
        project_rater(dir_output_merge, dir_output_utm, dir_shp, read_shp, proj)

        scale_x = float(input('\nInforme um fator de escala longitudinal: '))
        scale_y = float(input('Informe um fator de escala latitudinal: '))
        scale_factor = [scale_x, scale_y]

        mask_raster(dir_output_utm, dir_output_mask, dir_shp, read_shp, coords, proj, scale_factor)

        make_XYZ = input('\nGerar arquivo XYZ?\n1 - Sim\n2 - Não\nDigite uma Opção: ')

        if make_XYZ == '1':
            dir_output_csv = os.getcwd() + '/PontoCota.txt'
            extract_xyz(dir_output_utm, dir_shp, dir_output_csv, read_shp, coords, proj, scale_factor)
        
        make_curve = input('\nGerar Curvas de Nível?\n1 - Sim\n2 - Não\nDigite uma Opção: ')

        if make_curve == '1':
            delta_z = float(input('\nInforme o intervalo das curvas de nível: '))
            dir_output_dxf = os.getcwd() + '/Curva_de_Nível.dxf'
            extract_curves(dir_output_csv, dir_output_dxf, delta_z)

        if os.path.exists(dir_output_utm):
            os.remove(dir_output_utm)
        if os.path.exists(dir_output_merge):
            os.remove(dir_output_merge)
            
    if input_mode == '2':

        print('Função Temporariamente Desativada!')
        # path_raster = 'raster' 
        # exlcuir_raster(dir_raster)
        # output_path = Download_Topodata(dir_shp, read_shp, coords, proj, dir_raster)
        # Descompactar_zips(output_path)
        # dir_output_merge = output_path + '/merge_WGS.tif'
        # merge_raster(output_path, dir_output_merge)
        # dir_output_utm = output_path + '/merge_UTM.tif'
        # project_rater(dir_output_merge, dir_output_utm, dir_shp, read_shp, proj)

os.system('pause')