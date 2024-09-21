#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 21:32:28 2020

@author: christopher
"""

import os;
import numpy as np;
import subprocess as sp;
from reverse import reverse_bundle;
from shutil import copyfile, move;
import glob;
import utils.bundleTools as bt;
import utils.bundleTools3 as bt3;
import utils.visualizationTools as vt;

def read_intersection( infile ):
    
    f = open(infile, 'rb');

    total_triangles = np.frombuffer( f.read( 4 ), np.uint32 )[ 0 ];
    
    InTri = np.frombuffer( f.read( 4 * total_triangles), np.uint32 );
    FnTri = np.frombuffer( f.read( 4 * total_triangles), np.uint32 );
    
    InPoints = np.frombuffer( f.read( 4 * 3 * total_triangles), np.float32 ).reshape(-1,3);
    FnPoints = np.frombuffer( f.read( 4 * 3 * total_triangles), np.float32 ).reshape(-1,3);
    
    fib_index = np.frombuffer( f.read( 4 * total_triangles), np.uint32 );
    
    f.close();
    return InTri, FnTri, InPoints, FnPoints, fib_index;
    
#========================== Cálculo de intersección ==========================================    


ix_path = 'intersection/QB_20/001/R-R_direct/cluster_60_T1.intersectiondata'

ix = read_intersection(ix_path)

#%%
#
clustering_method = 'QB_20'
#
##---------- Sujeto a evaluar ----------#    
sub = '001';
#
##---------- Setup (mallado y creación de carpetas) ----------#
Lhemi_path = os.getcwd() + '/../../../../../ARCHI/Meshes_OBJ/' + sub + '/lh.obj'; # left hemisphere path
Rhemi_path = os.getcwd() + '/../../../../../ARCHI/Meshes_OBJ/' + sub + '/rh.obj'; # right hemisphere path

Lvertex, Lpolygons = bt.read_mesh_obj(Lhemi_path)
Rvertex, Rpolygons = bt.read_mesh_obj(Rhemi_path)

Lhemi = vt.Polygon(Lvertex, Lpolygons);
Rhemi = vt.Polygon(Rvertex, Rpolygons);

Lhemi.setOpacity(0.1);
Rhemi.setOpacity(0.1);

render = vt.Render();

render.AddActor(Lhemi);
render.AddActor(Rhemi);

int_type = ['L-L','L-R','R-L','R-R'];

for tipo in int_type:

    tipo = 'Inter'    

#    intersection_path = os.getcwd() + '/intersection/'+ clustering_method + '/' + sub + '/L-L_direct/'
    intersection_path = os.getcwd() + '/intersection/' + clustering_method + '/' + sub + '/Final/' + tipo + '/'; # intersection path
    
    intersection_list = os.listdir(intersection_path);

    int_class = intersection_path.split('/')[-2];
    
    int_class = 'L-R'
    
    for idx in range(len(intersection_list)):
    
        idx = 16;
#        idx = np.random.randint(0,len(intersection_list)-1);
        
        cluster_name = os.path.splitext(intersection_list[idx])[0];
        cluster_name = 'cluster_9_T1'
        bundles_path = os.getcwd() + '/../isaias_clustering/Brain_fiber_clustering/data/79subjects/inter-clustered_' + clustering_method + '/' + sub + '/individual_aligned_clusters/T1/' + cluster_name + '.bundles'; # bundles path
        cluster = bt.read_bundle(bundles_path);
        
        InTri, FnTri, InPoints, FnPoints, fib_index = read_intersection(intersection_path + cluster_name + '.intersectiondata');
        
        for fiber in cluster[0:int(0.2*len(cluster))]:
            fib = vt.Line(fiber);
            fib.setJetColormap(); # mapa de color Jet
            fib.setOpacity(1);
            render.AddActor(fib); # se agrega el objeto fib
        
        if int_class == 'L-L':
            in_tri = vt.Polygon(Lvertex, Lpolygons[InTri]);
            fin_tri = vt.Polygon(Lvertex, Lpolygons[FnTri]);
            in_tri.setColor((0.0, 1.0, 0.0));
            fin_tri.setColor((1.0, 1.0, 0.0));
            render.AddActor(in_tri);
            render.AddActor(fin_tri);
        
        elif int_class == 'L-R':
            in_tri = vt.Polygon(Lvertex, Lpolygons[InTri]);
            fin_tri = vt.Polygon(Rvertex, Rpolygons[FnTri]);
            in_tri.setColor((0.0, 1.0, 0.0));
            fin_tri.setColor((1.0, 1.0, 0.0));
            render.AddActor(in_tri);
            render.AddActor(fin_tri);
            
        elif int_class == 'R-L':
            in_tri = vt.Polygon(Rvertex, Rpolygons[InTri]);
            fin_tri = vt.Polygon(Lvertex, Lpolygons[FnTri]);
            in_tri.setColor((0.0, 1.0, 0.0));
            fin_tri.setColor((1.0, 1.0, 0.0));
            render.AddActor(in_tri);
            render.AddActor(fin_tri);
        
        elif int_class == 'R-R':
            in_tri = vt.Polygon(Rvertex, Rpolygons[InTri]);
            fin_tri = vt.Polygon(Rvertex, Rpolygons[FnTri]);
            in_tri.setColor((0.0, 1.0, 0.0));
            fin_tri.setColor((1.0, 1.0, 0.0));
            render.AddActor(in_tri);
            render.AddActor(fin_tri);
        break
    break
    
render.Start();

del render
    