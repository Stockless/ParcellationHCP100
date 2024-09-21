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
import bundleTools as bt;
import bundleTools3 as bt3;
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

sub = '037';
#
##---------- Setup (mallado y creación de carpetas) ----------#
Lhemi_path = os.getcwd() + '/../../../../../ARCHI/Meshes_OBJ/' + sub + '/lh.obj'; # left hemisphere path
Rhemi_path = os.getcwd() + '/../../../../../ARCHI/Meshes_OBJ/' + sub + '/rh.obj'; # right hemisphere path

Lvertex, Lpolygons = bt.read_mesh_obj(Lhemi_path)
Rvertex, Rpolygons = bt.read_mesh_obj(Rhemi_path)

Lnbors = bt.mesh_neighbors(Lpolygons)
Rnbors = bt.mesh_neighbors(Rpolygons)

Tri = 70273

Tri_nbors = Lnbors[Tri]

render = vt.Render()

tri = vt.Polygon(Lvertex,Lpolygons[[Tri]])
nbors = vt.Polygon(Lvertex,Lpolygons[Tri_nbors])

tri.setColor((1.0,0.0,0.0))
nbors.setColor((1.0,0.6,0.6))

Lhemi = vt.Polygon(Lvertex, Lpolygons);
Rhemi = vt.Polygon(Rvertex, Rpolygons);

Lhemi.setOpacity(1);
Rhemi.setOpacity(1);

render.AddActor(Lhemi);

render.AddActor(tri)
render.AddActor(nbors)

render.Start()
del render




#%%

sub = '037'

Lhemi_path = os.getcwd() + '/../../../../../ARCHI/Meshes_OBJ/' + sub + '/lh.obj'; # left hemisphere path
Rhemi_path = os.getcwd() + '/../../../../../ARCHI/Meshes_OBJ/' + sub + '/rh.obj'; # right hemisphere path

Lvertex, Lpolygons = bt.read_mesh_obj(Lhemi_path)
Rvertex, Rpolygons = bt.read_mesh_obj(Rhemi_path)

Lnbors = bt.mesh_neighbors(Lpolygons)
Rnbors = bt.mesh_neighbors(Rpolygons)


InTri, FnTri, InPoints, FnPoints, fib_index = read_intersection('intersection/HDB_1_652/' + sub + '/Final/L-L/cluster_113_T1.intersectiondata');

cluster = bt.read_bundle('bundle_classification/HDB_1_652/' + sub + '/L-L/cluster_113_T1.bundles')

if 70273 in InTri:
    print("YEAH!")
#%%
Lhemi = vt.Polygon(Lvertex, Lpolygons);
Rhemi = vt.Polygon(Rvertex, Rpolygons);

Lhemi.setOpacity(0.1);
Rhemi.setOpacity(1);

render = vt.Render();

#render.AddActor(Lhemi);
#render.AddActor(Rhemi);

punto = vt.Points(InPoints)
punto.setPointSize(5)
punto.setColor((1.0,0.5,0.0))
#render.AddActor(punto)

punto = vt.Points(FnPoints)
punto.setPointSize(5)
punto.setColor((0.0,0.5,1.0))
#render.AddActor(punto)

tri = vt.Polygon(Lvertex, Lpolygons[InTri])
tri.setColor((0.0,0.3,1.0))
#render.AddActor(tri)

tri = vt.Polygon(Lvertex, Lpolygons[FnTri])
tri.setColor((1.0,0.8,0.0))
#render.AddActor(tri)

for fiber in cluster:
    fib = vt.Line(fiber)
    fib.setJetColormap()
    render.AddActor(fib)
    break

render.Start()

del render



#%%
clustering_method = 'QB_202'
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

Lhemi.setOpacity(0)
Rhemi.setOpacity(0);

render = vt.Render();

render.AddActor(Lhemi);
#render.AddActor(Rhemi);

int_type = ['L-L','Inter','R-R'];

for tipo in int_type:

    
#    intersection_path = os.getcwd() + '/intersection/'+ clustering_method + '/' + sub + '/L-L_direct/'
    intersection_path = os.getcwd() + '/intersection/' + clustering_method + '/' + sub + '/Final/' + tipo + '/'; # intersection path
    
    intersection_list = os.listdir(intersection_path);

    int_class = intersection_path.split('/')[-2];
    
    for idx in range(len(intersection_list)):
    
        idx = 75;
#        idx = np.random.randint(0,len(intersection_list)-1);
        
        cluster_name = os.path.splitext(intersection_list[idx])[0];
#        cluster_name = 'cluster_9_T1'
        bundles_path = os.getcwd() + '/../isaias_clustering/Brain_fiber_clustering/data/79subjects/inter-clustered_' + clustering_method + '/' + sub + '/individual_aligned_clusters/T1/' + cluster_name + '.bundles'; # bundles path
        cluster = bt.read_bundle(bundles_path);
        
        InTri, FnTri, InPoints, FnPoints, fib_index = read_intersection(intersection_path + cluster_name + '.intersectiondata');
        
        fibra_idx = fib_index[0]
        
        for fiber in cluster[0:int(0.3*len(cluster))]:
#        for fiber in cluster[fibra_idx:fibra_idx+1]:    
            fib = vt.Line(fiber);
            fib.setOpacity(0.5)
#            fib.setJetColormap(); # mapa de color Jet
            fib.setColor((0.2,0.35,0.73))
            fib.setOpacity(1);
            render.AddActor(fib); # se agrega el objeto fib
    
#        for i, vertex in enumerate(Lpolygons[0:10000]):
#            v_coord = Lvertex[vertex]
#            side = vt.Line(v_coord)
#            side.setColor((0.0,0.0,0.0))
#            side.setOpacity(1)
#            render.AddActor(side)
        
        
        if tipo != 'Inter':
#            in_points = vt.Points(InPoints[0:int(0.3*len(cluster))])
            in_tri = vt.Polygon(Lvertex, Lpolygons[InTri])
            fn_tri = vt.Polygon(Lvertex, Lpolygons[FnTri])
            in_tri.setColor((0.31,1.0,0.0));
            fn_tri.setColor((1.0,0.94,0.0));
#            fn_points = vt.Points(FnPoints[0:int(0.3*len(cluster))])
#            in_points.setPointSize(5);
#            fn_points.setPointSize(5);
#            in_points.setColor((0.31, 1.0, 0.0));
#            fn_points.setColor((1.0, 0.94, 0.0));
#            render.AddActor(in_tri);
#            render.AddActor(fn_tri);
       
        else:
            in_points = vt.Points(InPoints)
            fn_points = vt.Points(FnPoints)
            in_points.setPointSize(3);
            fn_points.setPointSize(3);
            in_points.setColor((0.0, 1.0, 1.0));
            fn_points.setColor((1.0, 0.0, 1.0));
            render.AddActor(in_points);
            render.AddActor(fn_points);
            
        break
    break
            
render.Start();

del render
    
#%%

#bt3.write_bundle('216_onefiber.bundles',cluster[2:3])