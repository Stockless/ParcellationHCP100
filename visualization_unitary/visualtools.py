#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Agosto 10 20:09:15 2020

@author: fondecyt-1190701
"""
import numpy as np
import os
from collections import defaultdict

import bundleTools as bt
import visual_tools as vt
import vtk

import pickle

fasciculo = "D:/documentos/universidad/TESIS/HCP100/Tractografias/101006/aligned/left/aligned_lh_AR_ANT0.bundles"
fasciculo1 = "D:/documentos/universidad/TESIS/HCP100/Tractografias/101006/aligned/left/aligned_lh_PoCi-RAC_0.bundles"
fasciculo2 = "D:/documentos/universidad/TESIS/HCP100/Tractografias/101006/aligned/left/aligned_lh_PoCi-SF_0.bundles"
L_parcel_names = "D:/documentos/universidad/TESIS/parcellation/parcellation-master/visualization_unitary/Parcellation/ps_parcels/leftparcel_names.txt"
R_parcel_names = "D:/documentos/universidad/TESIS/parcellation/parcellation-master/visualization_unitary/Parcellation/ps_parcels/rightparcel_names.txt"

def AddBundle(fasciculo, fibras):
    tractografia = bt.read_bundle(fasciculo)
    all_points = [tractografia]

    for i in range(len(all_points)):
        a=0
        for fiber in all_points[i]:
            if(a<5):
                body = vt.Line(fiber);
                body.setJetScalar();
                fibras.append(body);
                a=a+1
            else:
                a=0
    return fibras

def load_restricted_triangles():
    #Existen ciertos triángulos que, por definición, no pueden representar parcelas. Estos corresponden a regiones posterior-inferior de la línea media...
    #...donde ambos hemisferios se unen, ya que en estricto rigor no corresponden a corteza, si no que a materia blanca/gris interna.
    
    #La parcelación obtenida por Lefranc indica cuales triángulos están "restringidos", los cuales están guardados en el archivo .pkl:
    #Lrestricted.pkl para el hemisferio izquierdo.
    #Rrestricted.pkl para el hemisferio derecho.
    
    with open('pickles/Lnope.pkl', 'rb') as f:
        Lrestricted = set(pickle.load(f))
    
    with open('pickles/Rnope.pkl', 'rb') as f:
        Rrestricted = set(pickle.load(f))
    
    return Lrestricted, Rrestricted

#Carga parcelas ('hard','cc','final')    
def load_parcels (parcels, parcels_path): #hard cc final
    Lparcels_path =  parcels_path + parcels +'_parcels/left/'
    Rparcels_path = parcels_path + parcels +'_parcels/right/'

    Lparcels = dict()
    Rparcels = dict()

    for Dir in os.listdir(Lparcels_path):
        Lparcels[Dir.split('.')[0]] = bt.read_parcels(Lparcels_path + Dir)
        
    for Dir in os.listdir(Rparcels_path):
        Rparcels[Dir.split('.')[0]] = bt.read_parcels(Rparcels_path + Dir)
    
    R_sp_empty = {k: v for k, v in Rparcels.items() if len(v) == 0}
    Rparcels = {k: v for k, v in Rparcels.items() if len(v) > 0}
    L_sp_empty = {k: v for k, v in Lparcels.items() if len(v) == 0}
    Lparcels = {k: v for k, v in Lparcels.items() if len(v) > 0}
    #Save list of empty parcels and delete them
    deleted_files = []
    folder_path = Rparcels_path
    for parcel_name in R_sp_empty.keys():
        file_path_bundle = os.path.join(folder_path,parcel_name+".parcelsdata")
        deleted_files.append(parcel_name)
        if os.path.exists(file_path_bundle):
            os.remove(file_path_bundle)
    df_path = parcels_path + parcels + '_parcels/R_deleted.txt'
    with open(df_path,'w') as f:
        for file_name in deleted_files:
            f.write(f'{file_name}\n')    #Save list of empty parcels and delete them
    deleted_files = []
    folder_path = Lparcels_path
    for parcel_name in L_sp_empty.keys():
        file_path_bundle = os.path.join(folder_path,parcel_name+".parcelsdata")
        deleted_files.append(parcel_name)
        if os.path.exists(file_path_bundle):
            os.remove(file_path_bundle)
    df_path = parcels_path + parcels + '_parcels/L_deleted.txt'
    with open(df_path,'w') as f:
        for file_name in deleted_files:
            f.write(f'{file_name}\n')
    return Lparcels, Rparcels

def create_vtk_polygon(points, polygons):
    # Create a vtkPoints object
    vtk_points = vtk.vtkPoints()
    for point in points:
        vtk_points.InsertNextPoint(point)

    # Create a vtkCellArray object for polygons
    vtk_polygons = vtk.vtkCellArray()
    for polygon in polygons:
        vtk_polygon = vtk.vtkPolygon()
        for vertex_index in polygon:
            vtk_polygon.GetPointIds().InsertNextId(vertex_index)
        vtk_polygons.InsertNextCell(vtk_polygon)

    # Create a vtkPolyData object
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.SetPolys(vtk_polygons)

    return polydata

def visualize_parcellation(meshes_path, L_sp, R_sp, sub, fibras, seed=False):
    # Cargar triángulos restringidos
    Lrestricted, Rrestricted = load_restricted_triangles()
    
    # Semilla para colores aleatorios
    if seed:
        np.random.seed(seed)
    
    # Parcelas finales a graficar
    final_parcels = set()
    for k in L_sp.keys():
        final_parcels.add(k)
    for k in R_sp.keys():
        final_parcels.add(k)
    fp = list(final_parcels)

    # Paleta de colores según cantidad de parcelas
    paleta = [(np.random.random(), np.random.random(), np.random.random()) for _ in range(len(fp))]

    # Directorios de los mallados corticales
    Lhemi_path = os.path.join(meshes_path, sub, 'lh.obj')  # left hemisphere path
    Rhemi_path = os.path.join(meshes_path, sub, 'rh.obj')  # right hemisphere path

    # Lectura de mallados
    Lvertex, Lpolygons = bt.read_mesh_obj(Lhemi_path)
    Rvertex, Rpolygons = bt.read_mesh_obj(Rhemi_path)
    
    Lhemi = create_vtk_polygon(Lvertex, Lpolygons)
    Rhemi = create_vtk_polygon(Rvertex, Rpolygons)

    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window.SetSize(1920, 1080)
    render_window.AddRenderer(renderer)
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    L_sp_biggest = {k: v for k, v in L_sp.items() if len(v)>20}
    num_items = len(L_sp_biggest)
    num_columns = int((num_items * (16 / 9)) ** 0.5)
    num_rows = int((num_items + num_columns - 1) // num_columns)
    i = 0


    # Para cada parcela del hemisferio izquierdo...
    for k, v in L_sp_biggest.items(): #Reemplazar por L_sp en lugar de L_SP_biggest para ver todas
        color = paleta[fp.index(k)]
        if len(v) == 0:
            continue
    
        v_restricted = list(set(v).difference(Lrestricted))
        sp_tri = create_vtk_polygon(Lvertex, Lpolygons[v_restricted])

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(sp_tri)
        
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)

        mapper2 = vtk.vtkPolyDataMapper()
        mapper2.SetInputData(Lhemi)
        
        actor2 = vtk.vtkActor()
        actor2.SetMapper(mapper2)
        actor2.GetProperty().SetOpacity(0.3)
        
        x = i % num_columns
        y = (i // num_columns) % num_rows
        z = i // (num_columns * num_rows)
        actor.SetPosition(x * 200, y * 200, z * 200)
        actor2.SetPosition(x * 200, y * 200, z * 200)
        actor.RotateX(90)
        actor.RotateZ(90)
        actor2.RotateX(90)
        actor2.RotateZ(90)

        # Añadir las fibras al renderizador
        fibra_actor = []
        j=0
        for fibra in fibras:
            fibra_actor.append(fibra._myPolygonActor)
            fibra_actor[j].SetPosition(x*200,y*200,z*200)
            fibra_actor[j].RotateX(90)
            fibra_actor[j].RotateZ(90)
            j += 1
        for act in fibra_actor:
            renderer.AddActor(act)
        # Create and position the 3D text actor
        text_source = vtk.vtkVectorText()
        file = open(L_parcel_names,"r")
        data = file.read()
        data_list = data.replace('\n',' ').split(' ')
        text_source.SetText(k)#.split('_')[1]+": "+data_list[int(k.split('_')[1])])
        text_mapper = vtk.vtkPolyDataMapper()
        text_mapper.SetInputConnection(text_source.GetOutputPort())
        text_actor = vtk.vtkFollower()
        text_actor.SetMapper(text_mapper)
        text_actor.SetScale(10)
        text_actor.SetPosition(x * 200-150, y * 200, z * 200+100)  # position above the actor
        text_actor.GetProperty().SetColor(1, 1, 1)  # black color
        text_actor.SetCamera(renderer.GetActiveCamera())  # Ensure the text faces the camera
        renderer.AddActor(actor)
        renderer.AddActor(actor2)
        renderer.AddActor(text_actor)
        i += 1
    renderer.SetBackground(0, 0, 0)
    render_window.Render()
    render_window_interactor.Start()
    del render_window_interactor, render_window, renderer

    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window.SetSize(1920, 1080)
    render_window.AddRenderer(renderer)
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    R_sp_biggest = {k: v for k, v in R_sp.items() if len(v)>20}
    #R_sp_biggest = {k: v for k, v in R_sp.items() if k == "parcel_PoCi0_PrCu0" or k == "parcel_PrCu0_PoCi0"}
    num_items = len(R_sp_biggest)
    num_columns = int((num_items * (16 / 9)) ** 0.5)
    num_rows = int((num_items + num_columns - 1) // num_columns)
    i = 0

    # Para cada parcela del hemisferio derecho...
    for k, v in R_sp_biggest.items():#Reemplazar por R_sp en lugar de R_SP_biggest para ver todas
        color = paleta[fp.index(k)]
        if len(v) == 0:
            continue

        v_restricted = list(set(v).difference(Rrestricted))
        sp_tri = create_vtk_polygon(Rvertex, Rpolygons[v_restricted])

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(sp_tri)
        
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)

        mapper2 = vtk.vtkPolyDataMapper()
        mapper2.SetInputData(Rhemi)
        
        actor2 = vtk.vtkActor()
        actor2.SetMapper(mapper2)
        actor2.GetProperty().SetOpacity(0.3)
        
        x = i % num_columns
        y = (i // num_columns) % num_rows
        z = i // (num_columns * num_rows)
        actor.SetPosition(x * 200, y * 200, z * 200)
        actor2.SetPosition(x * 200, y * 200, z * 200)
        actor.RotateX(90)
        actor.RotateZ(-90)
        actor2.RotateX(90)
        actor2.RotateZ(-90)

        # Create and position the 3D text actor
        text_source = vtk.vtkVectorText()
        file = open(R_parcel_names,"r")
        data = file.read()
        data_list = data.replace('\n',' ').split(' ')
        text_source.SetText(k)#.split('_')[1]+": "+data_list[int(k.split('_')[1])])
#
        text_mapper = vtk.vtkPolyDataMapper()
        text_mapper.SetInputConnection(text_source.GetOutputPort())
#
        text_actor = vtk.vtkFollower()
        text_actor.SetMapper(text_mapper)
        text_actor.SetScale(10)
        text_actor.SetPosition(x * 200+20, y * 200, z * 200)  # position above the actor
        text_actor.GetProperty().SetColor(1, 1, 1)  # black color
        text_actor.SetCamera(renderer.GetActiveCamera())  # Ensure the text faces the camera

        renderer.AddActor(actor)
        renderer.AddActor(actor2)
        renderer.AddActor(text_actor)
        i += 1

    renderer.SetBackground(0, 0, 0)
    render_window.Render()
    render_window_interactor.Start()
    del render_window_interactor, render_window, renderer




#Sujeto base. Puede ser cualquiera, ya que todos los mallados tienen triángulos correspondientes.
sub = '101006'
meshes_path= '../../../HCP100/Mallados/'
dice_thr=0.6

semilla_visualizacion = 48
final_parcels = 'Parcellation/'

fibras = []
#fibras = AddBundle(fasciculo,fibras)
#fibras = AddBundle(fasciculo1,fibras)
#fibras = AddBundle(fasciculo2,fibras)

#Lparcels_ps, Rparcels_ps= load_parcels('ps', final_parcels)
#Lparcels_fp, Rparcels_fp= load_parcels('fp', final_parcels)
#Lparcels_hard, Rparcels_hard= load_parcels('hard', final_parcels)
#Lparcels_cc, Rparcels_cc= load_parcels('cc', final_parcels)
Lparcels_final, Rparcels_final= load_parcels('final', final_parcels)

#visualize_parcellation(meshes_path, Lparcels_ps, Rparcels_ps, sub, fibras, seed = semilla_visualizacion)
#visualize_parcellation(meshes_path, Lparcels_fp, Rparcels_fp, sub,fibras, seed = semilla_visualizacion)
#visualize_parcellation(meshes_path, Lparcels_hard, Rparcels_hard, sub,fibras, seed = semilla_visualizacion)
#visualize_parcellation(meshes_path, Lparcels_cc, Rparcels_cc, sub, fibras, seed = semilla_visualizacion)
#visualize_parcellation(meshes_path, Lparcels_final, Rparcels_final, sub,fibras, seed = semilla_visualizacion)
