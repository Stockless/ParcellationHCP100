import sys
from tkinter import N
from unicodedata import name
import bundleTools as bt
import numpy as np
from queue import PriorityQueue as pq
from classes import *
import os

"""
Lee un fasciculo grande (e.g lh_AR) y lo segmenta en las regiones que 
este conecta según sus vertices y las etiquetas de desikan-killiany
"""

def label_triangles(triangles):

    for tri in triangles:
        label0 = tri.v1.label_parcel
        label1 = tri.v2.label_parcel
        label2 = tri.v3.label_parcel

        if (label0 == label1 == label2):
            tri.label_parcel = label2
        elif (label0 == label1 and label0!=label2):
            tri.label_parcel = label0
        elif (label1 == label2 and (label1!= label0)):
            tri.label_parcel = label1
        elif (label0 == label2) and (label0!=label1):
            tri.label_parcel = label0
        else:
            tri.label_parcel = label0

    return triangles

def read_intersection(file):
    with open(file,'r') as f:
        numTri = f.readline()
        InTri = list(map(int,f.readline().split()))
        FnTri = list(map(int,f.readline().split()))
        prob = list(map(float,f.readline().split()))
        prob_in = [[prob[i],prob[i+1],prob[i+2]] for i in range(0,len(prob),3)]
        prob = list(map(float,f.readline().split()))
        prob_fn = [[prob[i],prob[i+1],prob[i+2]] for i in range(0,len(prob),3)]
        id_fib = f.readline().split()
        f.close()
    return numTri, InTri, FnTri, prob_in, prob_fn, id_fib

def read_labels(labels_path):
    labels = []
    with open(labels_path, 'r') as file:
        for line in file:
            label = int(line)
            if label == -1:
                label = 0
            labels.append(label)
    return labels

def read_mesh_vtk(mesh_path,vertex_labels):
    vertices, triangles, vertex_coord, triangles_vertices = [],[],[],[]
    with open(mesh_path) as file:
        lines = file.readlines()
        vertex_index = 0
        triangle_index = 0
        for i,line in enumerate(lines):
            splitted = line.split(" ")
            id = str(splitted[0])
            if (id == 'v'):
                label = vertex_labels[vertex_index]
                vertex = Vertex(vertex_index,float(splitted[1]),float(splitted[2]),float(splitted[3]),label,[])
                vertex_index+=1
                vertices.append(vertex)
                coords = np.asarray([float(splitted[1]), float(splitted[2]), float(splitted[3])])
                vertex_coord.append(coords)
            if (id == 'f'):
                id_v1 = int(splitted[1].split("//")[0])
                id_v2 = int(splitted[2].split("//")[0])
                id_v3 = int(splitted[3].split("//")[0])
                v1 = vertices[id_v1-1]    
                v2 = vertices[id_v2-1]
                v3 = vertices[id_v3-1]
                triangle = Triangle(triangle_index,v1,v2,v3,-1,[],[])  
                triangles.append(triangle)  
                v1.triangles.append(triangle)   
                v2.triangles.append(triangle)
                v3.triangles.append(triangle)
                triangle_index+=1
                vs = np.asarray([int(splitted[1].split("//")[0]) - 1, int(splitted[2].split("//")[0]) - 1,
                                int(splitted[3].split("//")[0]) - 1])
                triangles_vertices.append(vs)
    file.close()
    return np.asarray(triangles),np.asarray(vertices)

def read_names(names):
    with open(names,'r') as f:
        parcel_names = []
        for line in f.readlines():
            parcel_names.append(line)
        f.close()
    return parcel_names

def get_triangle_fibers(intersection,triangles,names,bundle,thr=100):
    """
    Diccionario de conexión de parcelas para cada par de triangulos inicial y final
    contiene el nombre de las parcelas que conectan y el índice de la fibra correspondiente
    """
    conn_dict = {}
    name1 = bundle.split("_")[1]
    name2 = bundle.split("_")[2]
    if not name1+"_"+name2 in conn_dict:
        conn_dict[name1+"_"+name2] = []
    for i,tri in enumerate(intersection.InTri):
        conn_dict[name1+"_"+name2].append(intersection.id_fib[i])
    fibs_path = "fibers/"
    if not os.path.exists(fibs_path):
        os.makedirs(fibs_path)

    #Writes the fibers lists to a .txt file
    for name,fibs in conn_dict.items():
        # if len(fibs) < thr: #minimum number of fibers on bundle
        #     continue        #will only write bundles with significant fibers
        if not os.path.exists(fibs_path+bundle):
            os.makedirs(fibs_path+bundle)
        file = open(fibs_path+bundle+"/fibs_"+name+".txt",'w')
        for idx in fibs:
            file.write(str(idx)+" ")
        file.close()


#python segment_large_fascicles.py rh.obj large_fascicles/intersections/ rh_labels.txt dk_names.txt        

hemi_file = sys.argv[1] #rh.obj
hemi_intersection = sys.argv[2] #intersection files folder
hemi_labels = sys.argv[3] #rh_labels.txt
parcel_names = read_names(sys.argv[4]) #dk_names.txt

HemiVtxlabels = read_labels(hemi_labels)
HemiTriangles, Lvertex = read_mesh_vtk(hemi_file,HemiVtxlabels)

HemiTriangles = label_triangles(HemiTriangles)

for file in os.listdir(hemi_intersection):
    name = file.split(".")[:-1][0]
    print(name)
    numTri, InTri, FnTri, inter_in, inter_fn, fib_idx = read_intersection(hemi_intersection+file)
    inter = Intersection(numTri, InTri, FnTri, inter_in, inter_fn, fib_idx)

    get_triangle_fibers(inter,HemiTriangles,parcel_names,name)




