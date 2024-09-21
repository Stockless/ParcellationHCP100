#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Agosto 10 20:09:15 2020

@author: fondecyt-1190701
"""
import numpy as np
import os
import operator
import subprocess as sproc
import networkx as nx
from collections import defaultdict
from time import time
from itertools import combinations
from random import random
from copy import deepcopy
from statistics import mode
import re

import bundleTools as bt
import visualizationTools as vt
# from BundleTools import bundleTools3 as bt3

# from Parcellation import visualizationTools as vt

import pickle

def get_bundle_names(bundle):
    start = bundle.find("lh")
    if start == -1:
        start = bundle.find("rh")
    splitted = re.split("[._-]",bundle[start+3:])
    splitted.pop(-1)
    init_region = splitted[0]
    end_region = None
    if len(splitted) > 1:
        end_region = splitted[1]
        return [init_region,end_region]
    return [init_region]

def load_restricted_triangles():
    #Existen ciertos triángulos que, por definición, no pueden representar parcelas. Estos corresponden a regiones posterior-inferior de la línea media...
    #...donde ambos hemisferios se unen, ya que en estricto rigor no corresponden a corteza, si no que a materia blanca/gris interna.
    
    #La parcelación obtenida por Lefranc indica cuales triángulos están "restringidos", los cuales están guardados en el archivo .pkl:
    #Lrestricted.pkl para el hemisferio izquierdo.
    #Rrestricted.pkl para el hemisferio derecho.
    
    with open('Parcellation/Lnope.pkl', 'rb') as f:
        Lrestricted = set(pickle.load(f))
    
    with open('Parcellation/Rnope.pkl', 'rb') as f:
        Rrestricted = set(pickle.load(f))
    
    return Lrestricted, Rrestricted

#Carga parcelas ('hard','cc','final')    
def load_parcels (parcels, parcels_path): #hard cc final
    Lparcels_path =  parcels_path + '/'+ parcels +'_parcels/left/'
    Rparcels_path = parcels_path + '/'+ parcels +'_parcels/right/'

    Lparcels = dict()
    Rparcels = dict()

    for Dir in os.listdir(Lparcels_path):
        Lparcels[Dir.split('.')[0]] = bt.read_parcels(Lparcels_path + Dir)
        
    for Dir in os.listdir(Rparcels_path):
        Rparcels[Dir.split('.')[0]] = bt.read_parcels(Rparcels_path + Dir)
            
    return Lparcels, Rparcels

#Formula del coeficiente DICE
def DSC(A,B):
    set_A = set(A)
    set_B = set(B)
    
    interx = set_A.intersection(set_B)
    
    return (2*len(interx))/((len(set_A)+len(set_B)))

#Calculo de coeficiente DICE para comparar parcelas con atlas
def dise_comparision (atlas_comparision, parcels_path, dice_thr):
    # Lectura de parcelas propias.    
    lh_mine_dict, rh_mine_dict = load_parcels('final', parcels_path)
        
    # Cargar parcelación del atlas seleccionado.
    with open(atlas_comparision + '_Lparcels.pkl', 'rb') as f:   
        lh_atlas_dict = pickle.load(f)
    
    with open(atlas_comparision + '_Rparcels.pkl', 'rb') as f:   
        rh_atlas_dict = pickle.load(f)    
        
    # Cálculo de DSC entre parcelas.
    # Diccionario parcela_propia:parcela_atlas, para aquellos casos que superan el umbral de DSC.
    lh_dice_dict = defaultdict(int)
    rh_dice_dict = defaultdict(int)

    # Lista que almacena las parcelas ya utilizadas.
    used = []

    # Iteración sobre todas las parcelas propias obtenidas.
    for parcel, tris in lh_mine_dict.items():
        #Se registra el mayor DSC obtenido entre todas las comparaciones, y el nombre de la parcela asociada.
        dice_max = 0 
        winner = ''
      
        #Se itera sobre las parcelas del atlas.
        for ref, rtris in lh_atlas_dict.items():
            
            #Si la parcela ya fue utilizada, ignorar.
            if ref in used:
                continue
            
            #Cálculo de DSC entre parcela propia y parcela del atlas.
            dice = DSC(tris,rtris)
            
            #Si es mayor al máximo obtenido, actualizar.
            if dice > dice_max:
                dice_max = dice
                winner = ref
        
        #Si el mayor valor obtenido supera el umbral de similitud, se agrega al diccionario de parcelas similares, y a la lista de parcelas ya utilizadas.
        if dice_max >= dice_thr:
            lh_dice_dict[parcel] = winner
            used.append(winner)
            
        # Si no supera el umbral, continuar con la siguiente parcela.
        else:
            continue
    
    #Se repite lo mismo para el hemisferio derecho.
    used = []
    
    for parcel, tris in rh_mine_dict.items():
        dice_max = 0
        winner = ''
      
        for ref, rtris in rh_atlas_dict.items():
            if ref in used:
                continue
            dice = DSC(tris,rtris)
            if dice > dice_max:
                dice_max = dice
                winner = ref
            
        if dice_max >= dice_thr:
            rh_dice_dict[parcel] = winner
            used.append(winner)
            
        else:
            continue  
    
    #Diccionario de parcelas similares, considerando los triángulos asociados en el caso de la parcelación propia.
    lh_mine_common = {k:lh_mine_dict[k] for k in list(lh_dice_dict.keys()) if k in lh_mine_dict}
    rh_mine_common = {k:rh_mine_dict[k] for k in list(rh_dice_dict.keys()) if k in rh_mine_dict}
    
    #Diccionario de parcelas similares, considerando los triángulos asociados en el caso de la parcelación del atlas.
    lh_atlas_common = {k:lh_atlas_dict[k] for k in list(lh_dice_dict.values()) if k in lh_atlas_dict}
    rh_atlas_common = {k:rh_atlas_dict[k] for k in list(rh_dice_dict.values()) if k in rh_atlas_dict}
    
    if not os.path.exists(parcels_path+'Dise/'):
        os.makedirs(parcels_path+'Dise/')
        
    save_txt(parcels_path+'Dise/lh_list_'+atlas_comparision.replace("/",""),lh_mine_common, dice_thr)
    save_txt(parcels_path+'Dise/rh_list_'+atlas_comparision.replace("/",""), rh_mine_common, dice_thr)
    
    return lh_mine_common, rh_mine_common, lh_atlas_common, rh_atlas_common

def save_txt (name, diccionary, dice_thr):
    count=0
    with open(name+'_'+str(int(dice_thr*100))+'.txt', 'w') as f:
        for line in diccionary.keys():
            f.write(line)
            f.write('\n')
            count+=1
        f.write('\nTotal=%d' % count)
   #%%     

def visualize_parcellation(meshes_path, L_sp, R_sp, sub, seed = False):
    #Cargar trianguos restringidos
    Lrestricted, Rrestricted= load_restricted_triangles()
    
    #Semilla para colores aleatorios
    if seed != False:
        np.random.seed(seed)
    
    #Parcelas finales a graficar
    final_parcels = set()

    for k in L_sp.keys():
        final_parcels.add(k)
        
    for k in R_sp.keys():
        final_parcels.add(k)

    fp = list(final_parcels)    

    #Paleta de colores según cantidad de parcelas
    paleta = [(np.random.random(), np.random.random(), np.random.random()) for i in range(len(fp))]

    #Directorios de los mallados corticales
    Lhemi_path = meshes_path + sub + '/lh.obj'; # left hemisphere path
    Rhemi_path = meshes_path + sub + '/rh.obj'; # right hemisphere path

    #Lectura de mallados
    Lvertex, Lpolygons = bt.read_mesh_obj(Lhemi_path)
    Rvertex, Rpolygons = bt.read_mesh_obj(Rhemi_path)
    
    Lhemi = vt.Polygon(Lvertex, Lpolygons);
    Rhemi = vt.Polygon(Rvertex, Rpolygons);
    
    Lhemi.setOpacity(1);
    Rhemi.setOpacity(1);
    
    #Creación del render a visualizar
    render = vt.Render();
    
    #Se renderizan los mallados
    render.AddActor(Lhemi);
    render.AddActor(Rhemi);
    
    #Para cada parcela del hemisferio izquierdo...
    for k, v in L_sp.items():
        #Se selecciona un color de la paleta
        color = paleta[fp.index(k)]
        
        if len(v) == 0:
            continue
      
        #De todos los triángulos con parcelas, se eliminan aquellos que pertenecen al conjunto de triángulos restringidos.
        v_restricted = list(set(v).difference(Lrestricted))

        #Se renderizan los triángulos y polígonos.
        sp_tri = vt.Polygon(Lvertex, Lpolygons[v_restricted]);
        sp_tri.setColor((color[0], color[1], color[2]));
        render.AddActor(sp_tri);

    #Ídem para el derecho
    for k, v in R_sp.items():
        color = paleta[fp.index(k)]
    
        if len(v) == 0:
            continue
    
        v_restricted = list(set(v).difference(Rrestricted))
    
        sp_tri = vt.Polygon(Rvertex, Rpolygons[v_restricted]);
        sp_tri.setColor((color[0], color[1], color[2]));
        render.AddActor(sp_tri);
    
    render.Start();
    del render


def read_parcel_names(names_path):
    parcel_names = []
    anatomic_parcels = []
    with open(names_path, 'r') as file:
        for i,line in enumerate(file):
            parcel_names.append(str(line.rstrip("\n\r")))
    return parcel_names
t0 = time()
#Sujeto base. Puede ser cualquiera, ya que todos los mallados tienen triángulos correspondientes.
sub = '001'
theta_QB=20
intersection_path= 'inputs/intersection/'
meshes_path= '../../Data/meshes_obj/' 
#0.1idc - 0,2dc
idc_thr = [0.10, 0.20, 0.30, 0.40]
dc_thr = [0.10, 0.15, 0.20, 0.25, 0.30]
dice_thr=0.5
ero = 1
dil = 6
output_parcellation= 'Parcellation/IDC_'+str(idc_thr*100)+'%_DC_'+str(dc_thr*100)+'%/' #revisar
# Selección de atlas a comparar.
atlases = ['atlas/Lefranc','atlas/Brainnetome','atlas/Narciso']
semilla_visualizacion = 48
parcel_names = read_parcel_names('inputs/dk_names.txt')
start = 1
end = 3
#Parcelas preliminares
# preliminar_subparcels (start, end, sub, intersection_path, meshes_path,parcel_names)

#Parcelas duras
# Lparcels_hard, Rparcels_hard= hard_parcels (start, end, sub, intersection_path, meshes_path, parcel_names, dc_thr, idc_thr, output_parcellation)

#Parcelas finales
# Lparcel_cc, Rparcel_cc, Rparcels_final, Lparcels_final=  final_parcels (output_parcellation, meshes_path, sub, ero, dil)

#Comparacion Dice
lh_mine_common, rh_mine_common, lh_atlas_common, rh_atlas_common= dise_comparision(atlases[0], output_parcellation, dice_thr)

#Visualizacion de parcelas
# visualize_parcellation(meshes_path, Lparcels_hard, Rparcels_hard, '001', seed = semilla_visualizacion)
# visualize_parcellation(meshes_path, Lparcel_cc, Rparcel_cc, '001', seed = semilla_visualizacion)
#visualize_parcellation(meshes_path, Lparcels_final, Rparcels_final, '001', seed = semilla_visualizacion)

#Visualizacion para DISE
#
#visualize_parcellation(meshes_path, lh_atlas_common, rh_atlas_common, '001', seed = 47)
#visualize_parcellation(meshes_path, lh_mine_common, rh_mine_common, '001', seed = 47)

#Lparcels_final, Rparcels_final= load_parcels('final', output_parcellation)
#Lparcels_hard, Rparcels_hard= load_parcels('hard', output_parcellation)
# Lparcels_cc, Rparcels_cc= load_parcels('cc', output_parcellation)
#visualize_parcellation(meshes_path, Lparcels_final, Rparcels_final, sub, seed = semilla_visualizacion)
#visualize_parcellation(meshes_path, Lparcels_hard, Rparcels_hard, '001', seed = semilla_visualizacion)
#visualize_parcellation(meshes_path, Lparcels_cc, Rparcels_cc, '001', seed = semilla_visualizacion)
