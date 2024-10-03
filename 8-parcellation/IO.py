
# Copyright (C) 2019  Andrea Vázquez Varela

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#Authors:
# Narciso López López
# Andrea Vázquez Varela
#Creation date: 31/01/2019
#Last update: 20/07/2019


import numpy as np
import os
from classes import *
import utils
import shutil
import bundleTools as bt
import itertools
import sys
import json
import pickle
import copy


def create_atlas_dirs(path,trac):
    trac_path = ""
    if os.path.exists(path):
        shutil.rmtree(path)
    os.mkdir(path)
    atlas_path = path+"/final_atlas"
    os.mkdir(atlas_path)
    if trac.lower() == "y":
        trac_path = path+"/trac_atlas"
        os.mkdir(trac_path)
    hp_path = path+"/hard_parcels"
    os.mkdir(hp_path)
    cc_path = path+"/cc_parcels"
    os.mkdir(cc_path)
    fp_path = path+"/fp_parcels"
    os.mkdir(fp_path)
    ps_path = path+"/ps_parcels"
    os.mkdir(ps_path)
    return atlas_path,trac_path,cc_path,hp_path,fp_path,ps_path


def read_parcel_names(names_path):
    parcel_names = []
    anatomic_parcels = []
    with open(names_path, 'r') as file:
        for i,line in enumerate(file):
            parcel_names.append(str(line.rstrip("\n\r")))
            anatomicParcel = AnatomicParcel(i)
            anatomic_parcels.append(anatomicParcel)
    return(parcel_names, anatomic_parcels)

def read_vertex_labels(labels_path):
    labels = []
    with open(labels_path, 'r') as file:
        for line in file:
            label = int(line)
            if label == -1:
                label = 0
            labels.append(label)
    return(np.asarray(labels))


def read_mesh_obj(mesh_path,vertex_labels):
    vertices, triangles = [],[]
    with open(mesh_path) as file:
        lines = file.readlines()
        vertex_index = 0
        triangle_index = 0
        for line in lines:
            splitted = line.split(" ")
            id = str(splitted[0])
            if (id == 'v'):
                label = vertex_labels[vertex_index]
                vertex = Vertex(vertex_index,float(splitted[1]),float(splitted[2]),float(splitted[3]),label,[])
                vertex_index+=1
                vertices.append(vertex)
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
    file.close()
    return np.asarray(triangles)


def assign_preliminary_subparcels(infile,subject_name,name,parcel_names,triangles,anatomic_parcels,seen_tri):
    with open(infile) as f:
        content = f.readlines()
        if (content[0] != '0\n'):
            bundle_name = utils.get_bundle_names(name)
            if len(bundle_name) == 1:
                print("There are large fascicles. Segment them with the filtering tools and try again.")
                sys.exit()
            bundle_index = ""
            if len(bundle_name) == 3:
                bundle_index = "_"+str(bundle_name[2])

            """Gets the anatomic parcels based on bundle names"""
            anatomic_label1 = utils.find_label(bundle_name[0],parcel_names)
            anatomic_label2 = utils.find_label(bundle_name[1],parcel_names)
            anatomic_parcel1 = anatomic_parcels[anatomic_label1]
            anatomic_parcel2 = anatomic_parcels[anatomic_label2]

            """Generates names for preliminary subparcels 1 and 2"""
            if (bundle_name[0] != bundle_name[1]):
                parcel1_name = bundle_name[0]+"_"+bundle_name[1]+bundle_index
                parcel2_name = bundle_name[1]+"_"+bundle_name[0]+bundle_index
            else:
                parcel1_name = bundle_name[0] + "_0-" + bundle_name[1] + "_1" + bundle_index
                parcel2_name = bundle_name[1] + "_1-" + bundle_name[0] + "_0" + bundle_index

            """Creates or finds preliminary subparcels 1 and 2"""
            if parcel1_name in parcel_names:    #if subparcel is already created, finds it on the anatomic parcel
                parcel1_label = utils.find_label(parcel1_name,parcel_names)
                parcel1 = anatomic_parcel1.find_subparcel(parcel1_label)
            else:   #if not already created...
                parcel_names.append(parcel1_name) #adds it to the parcel names list
                """assigns an internal label based on the position where the subparcel is stored on the parcel names list"""
                parcel1_label = len(parcel_names) - 1
                parcel1 = SubParcel(parcel1_label,anatomic_label1,set(),[])
            """same for subparcel 2..."""
            if parcel2_name in parcel_names:
                parcel2_label = utils.find_label(parcel2_name,parcel_names)
                parcel2 = anatomic_parcel2.find_subparcel(parcel2_label)
            else:
                parcel_names.append(parcel2_name)
                parcel2_label = len(parcel_names) - 1
                parcel2 = SubParcel(parcel2_label,anatomic_label2,set(),[])

            """Adds bundle_label to the subparcel bundle_labels"""
            parcel1.add_bundle_label(parcel1_name)
            parcel2.add_bundle_label(parcel2_name)

            """Connects subparcel 1 to subparcel 2 and subparcel 2 to subparcel 1"""
            # parcel1.connected_parcels.append(parcel2)
            # parcel2.connected_parcels.append(parcel1)

            """Adds the subparcels to their corresponding anatomic parcel"""
            anatomic_parcel1.sub_parcels[parcel1.label] = parcel1
            anatomic_parcel2.sub_parcels[parcel2.label] = parcel2
#            
#            if parcel1.label not in subparcels_dict:
#                subparcels_dict[parcel1.label]={}
#            if parcel2.label not in subparcels_dict:
#                subparcels_dict[parcel2.label]={}

            """Gets intersection file data"""
            fibers = content[5].split(" ")[:-1]
            triangles_init = content[1].split(" ")[:-1]
            triangles_end = content[2].split(" ")[:-1]
            inter_inits = list(map(np.float32,content[3].split(" ")[:-1]))
            inter_ends = list(map(np.float32,content[4].split(" ")[:-1]))

            """Assigns triangles and fibers to the preliminary subparcels 1 and 2"""
            for i in range(len(fibers)):
                triangle_init = triangles[int(triangles_init[i])]
                triangle_end = triangles[int(triangles_end[i])]
                inter_init = utils.get_inter_coords(inter_inits,i)
                inter_end = utils.get_inter_coords(inter_ends,i)

                triangle_init.labels_subparcel.add(parcel1.label)
                triangle_end.labels_subparcel.add(parcel2.label)

                """We add the fiber to the fiber map on the triangle, adding it to the subparcel count"""
                if parcel1.label not in triangle_init.fibers_map:
                    triangle_init.fibers_map[parcel1.label] = 0
                triangle_init.fibers_map[parcel1.label] +=1
                if parcel2.label not in triangle_end.fibers_map:
                    triangle_end.fibers_map[parcel2.label] = 0
                triangle_end.fibers_map[parcel2.label] +=1

                """We add the fiber to the subjects map on the triangle, adding it to the subparcel cont
                only if it hasn't been seen previously"""
                if (triangle_init.index,parcel1.label) not in seen_tri:
                    seen_tri.add((triangle_init.index,parcel1.label))
                    if parcel1.label not in triangle_init.subjects_map:
                        triangle_init.subjects_map[parcel1.label]=0
                    triangle_init.subjects_map[parcel1.label] += 1
                if (triangle_end.index,parcel2.label) not in seen_tri:
                    seen_tri.add((triangle_end.index,parcel2.label))
                    if parcel2.label not in triangle_end.subjects_map:
                        triangle_end.subjects_map[parcel2.label]=0
                    triangle_end.subjects_map[parcel2.label] += 1

                if triangle_init not in parcel1.triangles:
                    parcel1.triangles.add(triangle_init)
                    parcel1.inter_points.append(inter_init)
                if triangle_end not in parcel2.triangles:
                    parcel2.triangles.add(triangle_end)
                    parcel2.inter_points.append(inter_end)


    return anatomic_parcels, parcel_names, seen_tri


def preliminary_subparcels(intersection_path,parcel_names,triangles,anatomic_parcels,hemi):
    if hemi == 'l':
        selected_hemi = "left-hemi"
    else:
        selected_hemi = "right-hemi"
    subject_dirs = os.listdir(intersection_path)
    subject_dirs.sort()
    triangles_dict = {}
    for subj_dir in subject_dirs:
        subject_name = subj_dir.split("-")[0]
        intersectionDir = os.listdir(intersection_path+"/"+subj_dir+"/"+subj_dir.split("_")[0]+"_"+selected_hemi)
        intersectionDir.sort()
        seen_tri = set()
        for file in intersectionDir:
            file_path = intersection_path+"/"+subj_dir+"/"+subj_dir.split("_")[0]+"_"+selected_hemi+"/"+file
            anatomic_parcels,parcel_names, seen_tri = assign_preliminary_subparcels(file_path,subject_name,file,parcel_names,triangles,anatomic_parcels,seen_tri)
    return anatomic_parcels,parcel_names



def rename_bundles(segmentation_path):
    subject_dirs = os.listdir(segmentation_path)
    subject_dirs.sort()
    for subj_dir in subject_dirs:
        for hemi in ["left-hemi", "right-hemi"]:
            segmentationDir = os.listdir(segmentation_path+"/"+subj_dir+"/"+hemi)
            segmentationDir.sort()
            for file in segmentationDir:
                file_split = file.split("-")
                if (len(file_split)>3):
                    index = file_split[3].split(".")[0]
                    extension = file_split[3].split(".")[1]
                    new_name = file_split[0]+"-"+file_split[1]+"-"+file_split[2]+"_"+index+"."+extension
                    original_path = segmentation_path+"/"+subj_dir+"/"+hemi+"/"+file
                    new_path = segmentation_path+"/"+subj_dir+"/"+hemi+"/"+new_name
                    os.rename(original_path, new_path)


def write_atlas_names(names,atlas_path,hemi):
    file_path = atlas_path+"/"+hemi+"atlas_names.txt"
    f = open(file_path, "w+")
    for i,n in enumerate(names):
        if n != "None":
            f.write(str(i)+" "+n+"\n")
    f.close()

def write_parcels(aparcels,atlas_path,hemi):
    file_path = atlas_path+"/"+hemi+"parcels.txt"
    f = open(file_path,"w+")
    for ap in aparcels:
        f.write("ap "+str(ap.label)+"\n")
        for k,sp in ap.sub_parcels.items():
            f.write("sp " + str(sp.label) + "\n")  
            # f.write("c")
            # for cp in sp.connected_parcels:
            #     f.write(" "+str(cp.label))
            # f.write("\n")
            f.write("t")
            for t in sp.triangles:  
                f.write(" "+str(t.index))
            f.write("\n")
            f.write("p")
            for t in sp.triangles:
                if sp.label in t.prob_map:
                    prob  = round(t.prob_map[sp.label],5)
                    f.write(" "+str(prob))
            f.write("\n")
    f.close()


def write_hparcels(aparcels,triangles,atlas_path,names,hemi):
    tri_map = utils.map_triangles(triangles)
    file_path = atlas_path+"/"+hemi+"hard_parcels.txt"
    f = open(file_path,"w+")
    new_hlabel = len(names)
    for ap in aparcels:
        ap_alltri = []
        ap_triangles = []
        if ap.label in tri_map:
            ap_alltri = tri_map[ap.label]
        f.write("ap "+str(ap.label)+"\n") 
        if len(ap.hard_parcels) > 0:
            for hp in ap.hard_parcels:
                if len(hp.triangles)>0:
                    f.write("hp " + str(hp.label) + "\n") 
                    f.write("t")
                    for t in hp.triangles:  
                        ap_triangles.append(t)
                        f.write(" "+str(t.index))
                    f.write("\n")
        else:
            f.write("hp " + str(new_hlabel) + "\n") 
            for tri in ap_alltri:
                f.write(" " + str(tri.index)) 
            new_hlabel+=1
            f.write("\n")

    f.close()

def write_preliminary_subparcels_base(aparcels,atlas_path,hemi,trac_path,parcel_names):
    "Save parcel_names"
    g = open(atlas_path+"/"+hemi+"parcel_names.txt","w")
    for name in parcel_names:
        g.write(name+"\n")
    g.close()
    "save preliminary subparcels"
    output_path = atlas_path+"/"+hemi
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    for ap in aparcels:
        if len(ap.sub_parcels) > 0:
            for sp in ap.sub_parcels.values():
                if len(sp.triangles)>0:
                    bt.write_parcels(output_path+"/subparcel_"+str(sp.label)+".parcelsdata",[tri.index for tri in sp.triangles])
                    bundle_labels = sp.get_all_bundle_labels()
                    bundles_path = trac_path+"/"+hemi+"_bundle_labels.txt"
                    f = open(bundles_path, "a")
                    f.write(str(sp.label)+" ")
                    f.write(", ".join(bundle_labels) + "\n")
    total_subparcels = 0
    aparcels_data = []
    for ap in aparcels:
        ap_data = {
            'label': int(ap.label),
            'sub_parcels': []
        }
        
        for subparcel in ap.sub_parcels.values():
            total_subparcels += 1
            subparcel_data = {
                'label': int(subparcel.label),
                'label_anatomic': int(subparcel.label_anatomic),
                'bundle_labels': list(subparcel.bundle_labels),
                'triangles': [],
                'inter_points': []
            }
            
            for inter_point in subparcel.inter_points:
                inter_point_list = inter_point.tolist()
                subparcel_data['inter_points'].append(inter_point_list)
            
            for triangle in subparcel.triangles:
                labels_subparcel_list = list(triangle.labels_subparcel)
                triangle_data = {
                    'index': int(triangle.index),
                    'v1':{
                            'index': int(triangle.v1.index),
                            'x': triangle.v1.x,
                            'y': triangle.v1.y,
                            'z': triangle.v1.z,
                            'label_parcel': int(triangle.v1.label_parcel)
                    },
                    'v2': {
                            'index': int(triangle.v2.index),
                            'x': triangle.v2.x,
                            'y': triangle.v2.y,
                            'z': triangle.v2.z,
                            'label_parcel': int(triangle.v2.label_parcel)
                    },
                    'v3': {
                            'index': int(triangle.v3.index),
                            'x': triangle.v3.x,
                            'y': triangle.v3.y,
                            'z': triangle.v3.z,
                            'label_parcel': int(triangle.v3.label_parcel)
                    },
                    'label_parcel': int(triangle.label_parcel),
                    'label_subparcel': labels_subparcel_list,
                    'fibers': triangle.fibers,
                    'neighbors_index' : [],
                    'neighbors_labels_subparcel': []
                }
                for tri in triangle.neighbors:
                    triangle_data['neighbors_index'].append(int(tri.index))
                    triangle_data['neighbors_labels_subparcel'].append(list(tri.labels_subparcel))
                subparcel_data['triangles'].append(triangle_data)
            
            ap_data['sub_parcels'].append(subparcel_data)
        
        aparcels_data.append(ap_data)
    json_string = json.dumps({'aparcels': aparcels_data})
    with open(atlas_path+"/"+hemi+".json","w") as p_file:
            p_file.write(json_string)

def triangle_copy(triangles,empty_vertex):
    v_triangles =[]
    for tri1 in triangles:
        new_tri = Triangle(tri1.index,empty_vertex,empty_vertex,empty_vertex,tri1.label_parcel,tri1.labels_subparcel,[])
        v_triangles.append(new_tri)
        del new_tri,tri1
    return v_triangles

def write_preliminary_subparcels(aparcels,atlas_path,hemi,trac_path,parcel_names):
    "Save parcel_names"
    g = open(atlas_path+"/"+hemi+"parcel_names.txt","w")
    for name in parcel_names:
        g.write(name+"\n")
    g.close()
    "save preliminary subparcels and create copy to save"
    output_path = atlas_path+"/"+hemi
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    empty_vertex = Vertex(0,0,0,0,"",[])
    for ap in aparcels:
        for i, sp in ap.sub_parcels.items():
            triangles = []
            bt.write_parcels(output_path+"/subparcel_"+str(sp.label)+".parcelsdata",[tri.index for tri in sp.triangles])
            bundle_labels = sp.get_all_bundle_labels()
            bundles_path = trac_path+"/"+hemi+"_bundle_labels.txt"
            f = open(bundles_path, "a")
            f.write(str(sp.label)+" ")
            f.write(", ".join(bundle_labels) + "\n")
    #        for tri in sp.triangles:
    #            v1_triangles = triangle_copy(tri.v1.triangles,empty_vertex)
    #            v2_triangles = triangle_copy(tri.v2.triangles,empty_vertex)
    #            v3_triangles = triangle_copy(tri.v3.triangles,empty_vertex)
    #            tri.v1.triangles = v1_triangles
    #            tri.v2.triangles = v2_triangles
    #            tri.v3.triangles = v3_triangles
    #            neighbors = triangle_copy(tri.neighbors,empty_vertex)
    #            tri.neighbors = set(neighbors)
    #            del v1_triangles, v2_triangles, v3_triangles, tri, neighbors
    #        del sp, triangles
    #    del ap
    #with open(atlas_path+"/"+hemi+".pickle","wb") as p_file:
    #        pickle.dump(aparcels,p_file)
            
def write_parcels_fussed(aparcels,atlas_path,hemi):
    output_path = atlas_path+"/"+hemi
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    for ap in aparcels:
        if len(ap.sub_parcels) > 0:
            for sp in ap.sub_parcels.values():
                if len(sp.triangles)>0:
                    bt.write_parcels(output_path+"/fussed_parcel_"+str(sp.label)+".parcelsdata",[tri.index for tri in sp.triangles])

def write_hparcels_triangles(aparcels,atlas_path,hemi):
    output_path = atlas_path+"/"+hemi
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    for ap in aparcels:
        if len(ap.hard_parcels) > 0:
            for hp in ap.hard_parcels:
                if len(hp.triangles)>0:
                    bt.write_parcels(output_path+"/hparcel_"+str(hp.label)+".parcelsdata",[tri.index for tri in hp.triangles])
        
def write_parcels_cc(cc_parcels,cc_path,hemi):
    if not os.path.exists(cc_path+"/"+hemi):
        os.makedirs(cc_path+"/"+hemi)
    for k, Tris in cc_parcels.items():
        bt.write_parcels(cc_path+"/"+hemi+"/parcel_"+str(k)+".parcelsdata",Tris)

def write_k(aparcels,path,hemi):
    file_path = path+"/"+hemi+"k.txt"
    file = open(file_path,"w+")
    for ap in aparcels:
        label = str(ap.label)
        k=0
        for hp in ap.hard_parcels:
            if len(hp.triangles)>0:
                k+=1
        file.write(label+" "+str(k)+"\n")
    file.close()


def write_atlas(anatomic_parcels,names,triangles,atlas_path,trac,hemi):
    if trac == False:
        write_k(anatomic_parcels,atlas_path,hemi)
        write_hparcels(anatomic_parcels,triangles,atlas_path,names,hemi)

    write_atlas_names(names,atlas_path,hemi)
    write_parcels(anatomic_parcels,atlas_path,hemi)
