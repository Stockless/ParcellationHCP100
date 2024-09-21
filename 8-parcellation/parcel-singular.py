import argparse
from time import time
import IO
from utils import *
from classes import *
import os

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

        parcel_names[sub_parcel1.label] = new_name
        sub_parcel_list.remove(sub_parcel_list[1])  
    return parcel_names

def remove_subparcel(aparcels,anatomic_parcel,subparcel):
    for k,sparcel in anatomic_parcel.sub_parcels.items():
        for tri in sparcel.triangles:
            tri.remove_labels([subparcel.label])
    anatomic_parcel.remove_subparcel(subparcel)

def unir_parcelas(aparcels,hemi):
    nombre_archivo = '/media/stocktan/DATA/documentos/universidad/TESIS/parcellation/parcellation-master/8-parcellation/output/trac_atlas_s/'+hemi+'fusion.txt'
    try:
        with open(nombre_archivo, 'r') as archivo:
            for linea in archivo:
                print(linea.strip())  #REALIZAR FUSION AQUI
    except FileNotFoundError:
        print(f"El archivo {nombre_archivo} no fue encontrado.")
    except Exception as e:
        print(f"Se produjo un error al leer el archivo: {e}")

def read_intersection(intersection_path,subj_dir,triangles,anatomic_parcels,parcel_names,selected_hemi):
    subject_name = subj_dir.split("-")[0]
    intersectionDir = os.listdir(intersection_path+"/"+subj_dir+"/"+subj_dir.split("_")[0]+"_"+selected_hemi)
    intersectionDir.sort()
    for file in intersectionDir:
        print(file)
        file_path = intersection_path+"/"+subj_dir+"/"+subj_dir.split("_")[0]+"_"+selected_hemi+"/"+file
        anatomic_parcels, aparcel_names = IO.assign_preliminary_subparcels(file_path,subject_name,file,parcel_names,triangles,anatomic_parcels)
    return anatomic_parcels, aparcel_names


def main():
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
    parser.add_argument('--dc-thr', type=float, default='0.15',help='Less probable triangles in a parcel (probability)')
    parser.add_argument('--idc', type=float, default='0.1',help='Percent of common triangles in the intersection of two density centers')
    parser.add_argument('--ero', type=int, default='1',help='Erosion threshold')
    parser.add_argument('--dil', type=int, default='6',help='Dilation threshold')
    args = parser.parse_args()
    #atlas_path,trac_path,cc_path,hp_path = IO.create_atlas_dirs(args.output_dir,args.traceability)

    #LEER ARCHIVOS
    trac = args.traceability.lower()
    Lparcel_names, Lanatomic_parcels = IO.read_parcel_names(args.parcel_names)
    Rparcel_names, Ranatomic_parcels = IO.read_parcel_names(args.parcel_names)
    Lvertex_labels = IO.read_vertex_labels(args.Lvlabels_file) 
    Rvertex_labels = IO.read_vertex_labels(args.Rvlabels_file)

    Ltriangles = IO.read_mesh_obj(args.LVtk_file,Lvertex_labels)
    Ltriangles = label_triangles(Ltriangles)  
    Rtriangles = IO.read_mesh_obj(args.RVtk_file, Rvertex_labels)
    Rtriangles = label_triangles(Rtriangles)

    idc = args.idc 
    dc_thr = args.dc_thr 
    size_thr = args.size_thr 
    print("Obtaining preliminary subparcels")
    subject_dirs = os.listdir(args.Intersection_dir)
    subject_dirs.sort()
    for subj_dir in subject_dirs:
        #Lanatomic_parcels,Lparcel_names = read_intersection(args.Intersection_dir,subj_dir,Ltriangles,Lanatomic_parcels,Lparcel_names,"left-hemi")
        #Ranatomic_parcels,Rparcel_names = read_intersection(args.Intersection_dir,subj_dir,Rtriangles,Ranatomic_parcels,Rparcel_names,"right-hemi")
        Lanatomic_parcels = unir_parcelas(Lanatomic_parcels,"L")
        Ranatomic_parcels = unir_parcelas(Ranatomic_parcels,"R")


if __name__ == '__main__':
    main()