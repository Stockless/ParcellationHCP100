import bundleTools3 as bt
import os
import visual_tools as vt
import numpy as np
from classes import *

def distance_point_to_vertex(point, vertex):
    return np.linalg.norm(point - np.array([vertex.x, vertex.y, vertex.z]))

def distance_point_to_triangle(point, triangle):
    distances = [
        distance_point_to_vertex(point, triangle.v1),
        distance_point_to_vertex(point, triangle.v2),
        distance_point_to_vertex(point, triangle.v3)
    ]
    return min(distances)

def find_closest_triangle(point, triangles):
    min_distance = float('inf')
    closest_triangle = None
    
    for triangle in triangles:
        distance = distance_point_to_triangle(point, triangle)
        if distance < min_distance:
            min_distance = distance
            closest_triangle = triangle
    
    return closest_triangle.label_parcel

def find_fiber_triangle_labels(fiber, triangles):
    labels = []
    start_point = np.array(fiber[0])
    end_point = np.array(fiber[-1])
    
    start_label = find_closest_triangle(start_point, triangles)
    end_label = find_closest_triangle(end_point, triangles)
    
    labels.append(start_label)
    labels.append(end_label)
    return labels
    

def label_triangles(triangles,parcel_names):

    for tri in triangles:
        label0 = tri.v1.label_parcel
        label1 = tri.v2.label_parcel
        label2 = tri.v3.label_parcel

        if (label0 == label1 == label2):
            tri.label_parcel = parcel_names[label2]
        elif (label0 == label1 and label0!=label2):
            tri.label_parcel = parcel_names[label0]
        elif (label1 == label2 and (label1!= label0)):
            tri.label_parcel = parcel_names[label1]
        elif (label0 == label2) and (label0!=label1):
            tri.label_parcel = parcel_names[label0]
        else:
            tri.label_parcel = parcel_names[label0]

    return triangles

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
            parcel_names.append(line.replace("\n",""))
        f.close()
    return parcel_names

labels_folder = "D:/documentos/universidad/TESIS/parcellation/parcellation-master/11-atlasDWM_segmentation/labels/"
parcel_names = read_names("D:/documentos/universidad/TESIS/parcellation/parcellation-master/11-atlasDWM_segmentation/dk_names.txt")
mesh_path = "D:/documentos/universidad/TESIS/parcellation/parcellation-master/11-atlasDWM_segmentation/rh.obj"

# For left hemisphere
HemiVtxlabels = read_labels(labels_folder+"rh_labels.txt")
HemiTriangles, Lvertex = read_mesh_vtk(mesh_path,HemiVtxlabels)
HemiTriangles = label_triangles(HemiTriangles,parcel_names)

bundle_path = "D:/documentos/universidad/TESIS/parcellation/parcellation-master/11-atlasDWM_segmentation/atlasMNI/bundles_DWM/MNI_atlas_rh_AR.bundles"
bundle = bt.read_bundle(bundle_path)


fiber_bundles_dict = {}
for fiber in bundle:
    labels = find_fiber_triangle_labels(fiber,HemiTriangles)
    key1 = labels[0] + "_" + labels[1]
    key2 = labels[1] + "_" + labels[0]
    
    if key1 in fiber_bundles_dict:
        fiber_bundles_dict[key1].append(fiber)
        print(key1)
    elif key2 in fiber_bundles_dict:
        fiber_bundles_dict[key2].append(fiber)
        print(key2)
    else:
        fiber_bundles_dict[key1] = [fiber]
        print(key1)

# Print the bundles
for key, fibers in fiber_bundles_dict.items():
    if len(fibers) > 10:
        bt.write_bundle( "D:/documentos/universidad/TESIS/parcellation/parcellation-master/11-atlasDWM_segmentation/new_atlasMNI/bundles0/MNI_atlas_rh_AR"+key+".bundles", fibers)
        print(f"Bundle {key} saved with {len(fibers)} fibers")