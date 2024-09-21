import sys
from unicodedata import name
import numpy as np
import os
import sys

#sys.path.append('../bundleTools')
import bundleTools as bt
from classes import *

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

def get_triangle_fibers(intersection,triangles,names,bundle,output_path,hemi,thr=100):
    """
    Diccionario de conexión de parcelas para cada par de triangulos inicial y final
    contiene el nombre de las parcelas que conectan y el índice de la fibra correspondiente
    """
    conn_dict = {}
    for i,tri in enumerate(intersection.InTri):
        name1 = names[triangles[tri].label_parcel].rsplit("\n")[0] #Initial region name
        name2 = names[triangles[intersection.FnTri[i]].label_parcel].rsplit("\n")[0] #Final region name
        if name2+"_"+name1 in conn_dict:
            name1, name2 = name2, name1
        if not name1+"_"+name2 in conn_dict:
            conn_dict[name1+"_"+name2] = []
        conn_dict[name1+"_"+name2].append(i)

    for name,inter in conn_dict.items():
        if len(inter) < thr: #minimum number of fibers on bundle
            continue        #will only write bundles with significant fibers
        file = open(output_path+"aligned_"+hemi+"_"+name+".txt",'w')
        #print(str(file))
        file.write(str(len(inter))+"\n")
        # print(len(inter))
        for idx in inter:
            file.write(str(intersection.InTri[idx])+" ")
        file.write("\n")
        for idx in inter:
            file.write(str(intersection.FnTri[idx])+" ")
        file.write("\n")
        for idx in inter:
            for i in range(3):
                file.write(str(intersection.inter_in[idx][i])+" ")
        file.write("\n")
        for idx in inter:
            for i in range(3):
                file.write(str(intersection.inter_fn[idx][i])+" ")
        file.write("\n")
        for idx in inter:
            file.write(str(intersection.id_fib[idx])+" ")
        file.write("\n")
        file.close()


#python segment_large_fascicles.py rh.obj large_fascicles/intersections/ rh_labels.txt dk_names.txt        
meshes_path = sys.argv[1] #Data/meshes_obj/
intersection_folder = sys.argv[2] #Data/intersections/
labels_folder = sys.argv[3] #labels/
parcel_names = read_names(sys.argv[4]) #atlas_parcel_names.txt

for sub in os.listdir(intersection_folder):    
    print(sub)
    # For right hemisphere
    HemiVtxlabels = read_labels(labels_folder+"rh_labels.txt")
    HemiTriangles, Lvertex = read_mesh_vtk(meshes_path+sub+'/rh.obj',HemiVtxlabels)

    HemiTriangles = label_triangles(HemiTriangles)
    intersection_path = intersection_folder+sub+'/'+sub+'_right-hemi/'
    output_path = (os.getcwd()).replace('\\','/')+"/segmentation_filter/"+sub+'/'+sub+'_right-hemi/'
    path_existence = os.path.exists(output_path)
    if not path_existence:
        os.makedirs(output_path)
    for file in os.listdir(intersection_path):
        name = file.split(".")[:-1][0]
        if len(name.split("_")) == 3:
            numTri, InTri, FnTri, inter_in, inter_fn, fib_idx = read_intersection(intersection_path+file)
            inter = Intersection(numTri, InTri, FnTri, inter_in, inter_fn, fib_idx)
            get_triangle_fibers(inter,HemiTriangles,parcel_names,name,output_path,'rh')
            os.remove(output_path+file)
    
    # For left hemisphere
    HemiVtxlabels = read_labels(labels_folder+"lh_labels.txt")
    HemiTriangles, Lvertex = read_mesh_vtk(meshes_path+sub+'/lh.obj',HemiVtxlabels)

    HemiTriangles = label_triangles(HemiTriangles)
    intersection_path = intersection_folder+sub+'/'+sub+'_left-hemi/'
    output_path = os.getcwd().replace('\\','/')+"/segmentation_filter/"+sub+'/'+sub+'_left-hemi/'
    path_existence = os.path.exists(output_path)
    if not path_existence:
        os.makedirs(output_path)
    for file in os.listdir(intersection_path):
        name = file.split(".")[:-1][0]
        if len(name.split("_")) == 3:
            numTri, InTri, FnTri, inter_in, inter_fn, fib_idx = read_intersection(intersection_path+file)
            inter = Intersection(numTri, InTri, FnTri, inter_in, inter_fn, fib_idx)
            get_triangle_fibers(inter,HemiTriangles,parcel_names,name,output_path,'lh')
            os.remove(output_path+file)
    


