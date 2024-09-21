
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


import argparse
from time import time
import IO
from utils import *
import networkx as nx
import subprocess as sproc
from classes import *
import os
import time
import sys
import math


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


def fusion(aparcels,anatomic_parcel,sub_parcel_list,parcel_names):
    while(len(sub_parcel_list)!=1):
        sub_parcel1 = sub_parcel_list[0]
        sub_parcel2 = sub_parcel_list[1]
        new_name = fusion_names(parcel_names[sub_parcel1.label],parcel_names[sub_parcel2.label])
        new_triangles = sub_parcel1.triangles.union(sub_parcel2.triangles)
        sub_parcel1.triangles = new_triangles
        anatomic_parcel.remove_subparcel(sub_parcel2) 
        for k,sub_parcel in anatomic_parcel.sub_parcels.items(): 
            for tri in sub_parcel.triangles:
                tri.replace_label(sub_parcel2.label,sub_parcel1.label)
                """Fuse fibers_map"""
                if sub_parcel1.label not in tri.fibers_map and sub_parcel2.label in tri.fibers_map:
                    tri.fibers_map[sub_parcel1.label] = 0
                if sub_parcel2.label in tri.fibers_map:
                    tri.fibers_map[sub_parcel1.label] += tri.fibers_map[sub_parcel2.label]
                    del tri.fibers_map[sub_parcel2.label]
        parcel_names[sub_parcel1.label] = new_name
        sub_parcel_list.remove(sub_parcel_list[1])  
    return parcel_names


def joinable_sparcels(clique, visited):
    subparcels = []
    for node in clique:
        if node.label < len(visited):
            if not visited[node.label]:
                subparcels.append(node)
                visited[node.label] = True
    return subparcels


def delete_reps(us,idc_matrix):
    cliques = []
    for clique in us:
        cliques.append(sorted(clique, key = lambda kv: kv.label))

    cliques = sorted(cliques, key = lambda kv: sum([sp.label for sp in kv]))

    rep_map = {}
    selected_cliques = {}
    for i,clique in enumerate(cliques):
        for sparcel in clique:
            if sparcel not in rep_map:
                rep_map[sparcel] = [i]
            else:
                rep_map[sparcel].append(i)
    for sparcel1,clique_list in rep_map.items():
        if len(clique_list) > 1:
            max_idc = 1
            selected_cliques[sparcel1] = -1
            for i in clique_list:
                clique = cliques[i]
                sum_idc = 0
                for sparcel2 in clique:
                    sum_idc+=idc_matrix[sparcel1.label][sparcel2.label]
                if sum_idc > max_idc:
                    max_idc = sum_idc
                    selected_cliques[sparcel1] = i

    for sparcel,index_clique in selected_cliques.items():
        for i in range(len(cliques)):
            if i !=index_clique and sparcel in cliques[i]:
                cliques[i].remove(sparcel)


    return (sorted(cliques, key = len, reverse = True))

def remove_less_representative_triangles(subparcel):
    triangles_to_remove = []
    thr = 1
    for triangle in subparcel.triangles:
        if subparcel.label in triangle.subjects_map and triangle.subjects_map[subparcel.label]< thr:
            triangles_to_remove.append(triangle.index)
            del triangle.fibers_map[subparcel.label]
            triangle.labels_subparcel.remove(subparcel.label)
    for triangle in triangles_to_remove:
        subparcel.remove_triangle(triangle)
        

def create_fusion_list_nuevo_thr(anatomic_parcel,dc_thr,thr_idc,names,hemi):
    fusion_list = []
    nsparcels = len(names)
    idc_matrix = np.zeros((nsparcels,nsparcels))
    connect_graph = nx.Graph()
    for k,sparcel in anatomic_parcel.sub_parcels.items():
        connect_graph.add_node(sparcel)
    for k,sparcel1 in anatomic_parcel.sub_parcels.items():
        for k,sparcel2  in anatomic_parcel.sub_parcels.items():
            if sparcel1!=sparcel2:
                triangles1 = [tri for tri in sparcel1.triangles if tri.prob_map[sparcel1.label]>0.1725*math.exp(-0.07*len(tri.prob_map))]
                triangles2 = [tri for tri in sparcel2.triangles if tri.prob_map[sparcel2.label]>0.1725*math.exp(-0.07*len(tri.prob_map))]
                inter = intersection(triangles1, triangles2)
                min_inter = min(len(sparcel1.triangles),len(sparcel2.triangles))
                idc = len(inter) / min_inter
                if idc >= thr_idc: #Si se solapan mucho las subparcelas...
                    connect_graph.add_edge(sparcel1,sparcel2)
                    idc_matrix[sparcel1.label][sparcel2.label] = idc
                    idc_matrix[sparcel2.label][sparcel1.label] = idc

    cliques = sorted(nx.find_cliques(connect_graph), key=len, reverse=True)
    cliques = delete_reps(cliques,idc_matrix)
    
    visited = np.zeros(len(names), bool)
    for clique in cliques:
        parcel_list = joinable_sparcels(clique,visited)
        fusion_list.append(parcel_list)
    return fusion_list

def create_fusion_list(anatomic_parcel,dc_thr,thr_idc,names,hemi):
    fusion_list = []
    nsparcels = len(names)
    idc_matrix = np.zeros((nsparcels,nsparcels))
    connect_graph = nx.Graph()
    for k,sparcel in anatomic_parcel.sub_parcels.items():
        connect_graph.add_node(sparcel)
    for k,sparcel1 in anatomic_parcel.sub_parcels.items():
        for k,sparcel2  in anatomic_parcel.sub_parcels.items():
            if sparcel1!=sparcel2:
                triangles1 = sparcel1.get_triangles_prob(dc_thr, sparcel1.label)
                triangles2 = sparcel2.get_triangles_prob(dc_thr, sparcel2.label)
                inter = intersection(triangles1, triangles2)
                min_inter = min(len(sparcel1.triangles),len(sparcel2.triangles))
                idc = len(inter) / min_inter
                if idc >= thr_idc: #Si se solapan mucho las subparcelas...
                    connect_graph.add_edge(sparcel1,sparcel2)
                    idc_matrix[sparcel1.label][sparcel2.label] = idc
                    idc_matrix[sparcel2.label][sparcel1.label] = idc

    cliques = sorted(nx.find_cliques(connect_graph), key=len, reverse=True)
    cliques = delete_reps(cliques,idc_matrix)
    
    visited = np.zeros(len(names), bool)
    for clique in cliques:
        parcel_list = joinable_sparcels(clique,visited)
        fusion_list.append(parcel_list)
    return fusion_list


def recalc_probability(anatomic_parcel):
    for k,sub_parcel in anatomic_parcel.sub_parcels.items():
        for tri in sub_parcel.triangles:
            tri.set_prob_map()


def remove_subparcel(aparcels,anatomic_parcel,subparcel):
    for k,sparcel in anatomic_parcel.sub_parcels.items():
        for tri in sparcel.triangles:
            tri.remove_labels([subparcel.label])
            if subparcel.label in tri.fibers_map:
                del tri.fibers_map[subparcel.label]
    anatomic_parcel.remove_subparcel(subparcel)


def remove_small_parcels(aparcels,anatomic_parcel,parcel_names,size_thr,trac,trac_path,hemi):
    avg_inter = 0
    parcel_sizes = [len(sparcel.triangles) for k,sparcel in anatomic_parcel.sub_parcels.items()]
    avg_size = sum(parcel_sizes) / (len(anatomic_parcel.sub_parcels))
    ed_parcels = np.std(parcel_sizes)
    thr = max(0,avg_size - size_thr*ed_parcels)
    if trac == "y":
        remove_file = open(trac_path+"/"+hemi+"remove.txt","a+")
    inter_sizes = []
    for k,sparcel in anatomic_parcel.sub_parcels.items():
        inter_sizes += [len(tri.labels_subparcel) for tri in sparcel.triangles]
        avg_inter +=sum([len(tri.labels_subparcel) for tri in sparcel.triangles])
    avg_inter = avg_inter / (len(anatomic_parcel.sub_parcels))
    ed_inter = np.std(inter_sizes)
    thr_inter = max(0,avg_inter - size_thr*ed_inter)
    n = 0
    nap_flag = 0
    sp_to_remove = []
    for k,sparcel in anatomic_parcel.sub_parcels.items():
        num_inters = sum([len(tri.labels_subparcel) for tri in sparcel.triangles])
        if (len(sparcel.triangles) < thr) or (num_inters < thr_inter):
            n+=1
            if trac == "y":
                if nap_flag == 0:
                    remove_file.write("ap "+str(anatomic_parcel.label)+"\n")
                    nap_flag = 1
                remove_file.write(str(sparcel.label)+"\n")
            sp_to_remove.append(sparcel)
    for sparcel in sp_to_remove:        
        remove_subparcel(aparcels,anatomic_parcel,sparcel)
    return n

def remove_small_parcels_original(aparcels,anatomic_parcel,parcel_names,size_thr,trac,trac_path,hemi):
    avg_inter = 0
    avg_size = sum([len(sparcel.triangles) for k,sparcel in anatomic_parcel.sub_parcels.items()]) / (len(anatomic_parcel.sub_parcels))
    if trac == "y":
        remove_file = open(trac_path+"/"+hemi+"remove.txt","a+")
    for k,sparcel in anatomic_parcel.sub_parcels.items():
        avg_inter +=sum([len(tri.labels_subparcel) for tri in sparcel.triangles])
    avg_inter = avg_inter / (len(anatomic_parcel.sub_parcels))
    thr = (avg_size*size_thr)
    thr_inter = (avg_inter*size_thr)
    n = 0
    nap_flag = 0
    sp_to_remove = []
    for k,sparcel in anatomic_parcel.sub_parcels.items():
        num_inters = sum([len(tri.labels_subparcel) for tri in sparcel.triangles])
        if (len(sparcel.triangles) < thr) or (num_inters < thr_inter):
            n+=1
            if trac == "y":
                if nap_flag == 0:
                    remove_file.write("ap "+str(anatomic_parcel.label)+"\n")
                    nap_flag = 1
                remove_file.write(str(sparcel.label)+"\n")
            sp_to_remove.append(sparcel)
    for sparcel in sp_to_remove:        
        remove_subparcel(aparcels,anatomic_parcel,sparcel)
    return n


def processing_parcels(aparcels,idc,dc_thr,size_thr,parcel_names,trac,trac_path,hemi):
    n_remove = 0
    for i,anatomic_parcel in enumerate(aparcels):
        if len(anatomic_parcel.sub_parcels) > 0:
            if trac == "y":
                fusion_file = open(trac_path+"/"+hemi+"fusion.txt","a+")
            """Remove less representative triangle from subparcels"""
            for i, subparcel in anatomic_parcel.sub_parcels.items():
                remove_less_representative_triangles(subparcel)
                if len(subparcel.triangles) == 0:
                    anatomic_parcel.remove_subparcel(subparcel) #debo eliminarlas después
            """Density center calculation"""
            recalc_probability(anatomic_parcel)
            n_remove += remove_small_parcels(aparcels,anatomic_parcel,parcel_names,size_thr,trac,trac_path,hemi)
            """Parcel overlapping"""
            fusion_list = create_fusion_list(anatomic_parcel,dc_thr,idc,parcel_names,hemi)
            if trac == "y":
                fusion_file.write("ap "+str(anatomic_parcel.label)+"\n")
            for list in fusion_list:
                if len(list)>0:
                    if trac == "y":
                        fusion_file.write(" ".join([str(parcel.label) for parcel in list])+"\n")
                if (len(list)>1):
                    parcel_names = fusion(aparcels,anatomic_parcel,list,parcel_names)
    print("Removed subparcels: ",n_remove)
    return aparcels, parcel_names

def get_hard_parcels(aparcels,parcel_names,trac,trac_path,hemi):
    if trac == "y":
        probmap_file = open(trac_path+"/"+hemi+"probmap.txt","a+")
    for ap in aparcels:
        recalc_probability(ap)
        hparcels_map = {}
        for k,sp in ap.sub_parcels.items():
            for tri in sp.triangles:
                selected_parcel = most_probable(tri.prob_map)
                if trac == "y":
                    probmap_file.write(str(tri.index)+" "+str(selected_parcel)+"\n")
                if selected_parcel != -1:
                    selected_parcel = ap.find_subparcel(selected_parcel)
                    if selected_parcel != None:
                        if not selected_parcel.label in hparcels_map:
                            hparcels_map[selected_parcel.label] = set()
                        hparcels_map[selected_parcel.label].add(tri)
        for label, triangles in hparcels_map.items():
            sp = ap.find_subparcel(label)
            hparcel = SubParcel(label, ap.label, triangles, [])
            ap.add_hparcel_triangles(parcel_names[int(label)],triangles)
            ap.hard_parcels.add(hparcel)
    return aparcels

def get_PCC(aparcels,triangles):
    parcel_cc = {}
    triangles_arr = np.asarray([[tri.v1.index,tri.v2.index,tri.v3.index] for tri in triangles])
    t1 = time.time()
    n = 0
    for ap in aparcels:
        n += len(ap.hard_parcels)
        for k,tris in ap.hparcels.items():
            #parcel_cc[k] = PCC(tris,triangles_arr)
            tris = np.array(list(tris))
            if len(tris) == 0:
                parcel_cc[k] = tris
            else:
                edges_poly = []
                for ind in tris:
                    edges_poly.append([triangles_arr[ind,0],triangles_arr[ind,1]])
                    edges_poly.append([triangles_arr[ind,0],triangles_arr[ind,2]])
                    edges_poly.append([triangles_arr[ind,1],triangles_arr[ind,2]])

                edges_poly = np.unique(edges_poly,axis=0)

                G = nx.Graph();
                G.add_edges_from(edges_poly);
                
                cc = list(nx.connected_components(G));
                len_cc = [len(comp) for comp in cc];
                len_cc_thr = int(np.max(len_cc))
                regions_cc = list(cc[np.argmax(len_cc)]);
                regions_cc_2 = [list(cc[i]) for i, comp in enumerate(cc) if len_cc[i] >= len_cc_thr]
                #print(len(regions_cc_2))
                if len(regions_cc_2) > 1:
                    for n_cc, region_cc in enumerate(regions_cc_2):
                        ind = [np.where(triangles_arr[tris] == vertex)[0] for vertex in region_cc]
                        ind = np.unique(np.concatenate(ind))
                        parcel_cc[k + "_" + str(n_cc)] = tris[ind]
                else:
                    ind = [np.where(triangles_arr[tris]==region_cc)[0] for region_cc in regions_cc];
                    ind = np.unique(np.concatenate(ind));
                    parcel_cc[k] = tris[ind]

    t2 = time.time()
    # print("PCC time: ",t2-t1)
    print("Total hard parcels:",n)
    return parcel_cc

def PCC(Tri,triangles):
    Tri = np.array(list(Tri))
    if len(Tri) == 0:
        return Tri
    edges_poly = []
    for ind in Tri:
        edges_poly.append([triangles[ind,0],triangles[ind,1]])
        edges_poly.append([triangles[ind,0],triangles[ind,2]])
        edges_poly.append([triangles[ind,1],triangles[ind,2]])

    edges_poly = np.unique(edges_poly,axis=0)

    G = nx.Graph();
    G.add_edges_from(edges_poly);
    
    cc = list(nx.connected_components(G));
    len_cc = [len(comp) for comp in cc];
    regions_cc = list(cc[np.argmax(len_cc)]);
    
    ind = [np.where(triangles[Tri]==region_cc)[0] for region_cc in regions_cc];
    print(regions_cc)
    ind = np.unique(np.concatenate(ind));
    return Tri[ind]

def main():
    sys.setrecursionlimit(1750)
    parser = argparse.ArgumentParser(description='Create parcel of multiple subjects')
    parser.add_argument('--Intersection-dir',type= str,help='Input intersection directory')
    parser.add_argument('--LVtk-file', type=str, help='Input file with the vtk')
    parser.add_argument('--RVtk-file', type=str, help='Input file with the vtk')
    parser.add_argument('--Lvlabels-file',type= str,help='Input file with the vertex labels')
    parser.add_argument('--Rvlabels-file',type= str,help='Input file with the vertex labels')
    parser.add_argument('--parcel-names',type= str, help='Input file with the names of the parcels')
    parser.add_argument('--output-dir', type=str, help='Output directory')
    parser.add_argument('--traceability',type= str, default='y', help='Write y, to obtain the traceability of the parcels')
    parser.add_argument('--size-thr', type=float, default='0.1',help='Size to delete small parcels')
    parser.add_argument('--dc-thr', type=float, default='0.25',help='Less probable triangles in a parcel (probability)')
    parser.add_argument('--idc', type=float, default='0.1',help='Percent of common triangles in the intersection of two density centers')
    parser.add_argument('--ero', type=int, default='1',help='Erosion threshold')
    parser.add_argument('--dil', type=int, default='4',help='Dilation threshold')
    args = parser.parse_args()
    start = time.time()
    atlas_path,trac_path,cc_path,hp_path,fp_path,ps_path = IO.create_atlas_dirs(args.output_dir,args.traceability)

    trac = args.traceability.lower()
    t1 = time.time()
    Lparcel_names, Lanatomic_parcels = IO.read_parcel_names(args.parcel_names)
    Rparcel_names, Ranatomic_parcels = IO.read_parcel_names(args.parcel_names)
    Lvertex_labels = IO.read_vertex_labels(args.Lvlabels_file) 
    Rvertex_labels = IO.read_vertex_labels(args.Rvlabels_file)

    Ltriangles = IO.read_mesh_obj(args.LVtk_file,Lvertex_labels)
    Ltriangles =  label_triangles(Ltriangles)  
    Rtriangles = IO.read_mesh_obj(args.RVtk_file, Rvertex_labels)
    Rtriangles = label_triangles(Rtriangles)

    print("Obtaining preliminary subparcels")
    Lanatomic_parcels, Lparcel_names = IO.preliminary_subparcels(args.Intersection_dir, Lparcel_names,Ltriangles,Lanatomic_parcels,'l') 
    Ranatomic_parcels, Rparcel_names = IO.preliminary_subparcels(args.Intersection_dir, Rparcel_names, Rtriangles,Ranatomic_parcels, 'r')
    t2 = time.time()
    print("Time for IO operations: ",t2-t1)
    idc = args.idc 
    dc_thr = args.dc_thr 
    size_thr = args.size_thr 
    if trac == "y":
        for aparcel in Lanatomic_parcels:
            recalc_probability(aparcel)
        for aparcel in Ranatomic_parcels:
            recalc_probability(aparcel)
        IO.write_atlas(Lanatomic_parcels, Lparcel_names,Ltriangles, trac_path,True,"L")
        IO.write_atlas(Ranatomic_parcels, Rparcel_names,Rtriangles, trac_path,True,"R")
    IO.write_preliminary_subparcels(Lanatomic_parcels,ps_path,"left",trac_path,Lparcel_names)
    IO.write_preliminary_subparcels(Ranatomic_parcels,ps_path,"right",trac_path,Rparcel_names)
    
    t1 = time.time()

    print("Processing parcels")
    Lanatomic_parcels, Lparcel_names = processing_parcels(Lanatomic_parcels,idc,dc_thr,size_thr,Lparcel_names,trac,trac_path,"L")
    Ranatomic_parcels, Rparcel_names = processing_parcels(Ranatomic_parcels,idc,dc_thr,size_thr,Rparcel_names,trac,trac_path,"R")
    IO.write_parcels_fussed(Lanatomic_parcels,fp_path,"left")
    IO.write_parcels_fussed(Ranatomic_parcels,fp_path,"right")

    print("Obtaining hard parcels")
    Lanatomic_parcels = get_hard_parcels(Lanatomic_parcels,Lparcel_names,trac,trac_path,"L")
    Ranatomic_parcels = get_hard_parcels(Ranatomic_parcels,Rparcel_names,trac,trac_path,"R")
    IO.write_hparcels_triangles(Lanatomic_parcels,hp_path,"left")
    IO.write_hparcels_triangles(Ranatomic_parcels,hp_path,"right")

    print("Obtaining principal connected components")
    Lparcel_cc = get_PCC(Lanatomic_parcels,Ltriangles)
    Rparcel_cc = get_PCC(Ranatomic_parcels,Rtriangles)
    IO.write_parcels_cc(Lparcel_cc,cc_path,"left")
    IO.write_parcels_cc(Rparcel_cc,cc_path,"right")

    t2 = time.time()
    print("Time to process the atlas: ", t2-t1)

    IO.write_atlas(Lanatomic_parcels,Lparcel_names,Ltriangles,atlas_path,False,"L")
    IO.write_atlas(Ranatomic_parcels,Rparcel_names,Rtriangles,atlas_path,False,"R")

    Lparcelcc_path = cc_path + "/left/"
    Rparcelcc_path = cc_path + "/right/"
    Lhemi_path = args.LVtk_file
    Rhemi_path = args.RVtk_file

    Lfinal_path = args.output_dir + '/final_parcels/left/'
    Rfinal_path = args.output_dir + '/final_parcels/right/'
    
    if not os.path.exists(Lfinal_path):
        os.makedirs(Lfinal_path)
    
    if not os.path.exists(Rfinal_path):
        os.makedirs(Rfinal_path)
    print("Morphologic opening")
    sproc.call(['make']);
    sproc.call(['./main', Lhemi_path, Rhemi_path, str(args.ero), str(args.dil),  Lparcelcc_path, Rparcelcc_path, Lfinal_path, Rfinal_path]);

    finish = time.time()
    print("Parcellation finished")
    print("Total time: ", finish-start)

if __name__ == '__main__':
    main()