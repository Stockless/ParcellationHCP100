import os
import numpy as np
import sys
import re
from dipy.segment.metric import mdf
from dipy.segment.metric import dist
from dipy.segment.metric import EuclideanMetric

#Imports bundleTools3 from parent directory
#sys.path.append('../bundleTools')
import bundleTools3 as bt

"""Alings fibers of a bundle """
def align_bundle_1(bundle, ref):
    ref_start, ref_end = ref[0], ref[-1]
    aligned_bundle = []
    
    for fiber in bundle:
        fiber_start, fiber_end = fiber[0], fiber[-1]
        
        dist_starts = np.linalg.norm(ref_start - fiber_start)
        dist_ends = np.linalg.norm(ref_end - fiber_end)
        dist_end_start = np.linalg.norm(ref_start - fiber_end)
        dist_start_end = np.linalg.norm(ref_end - fiber_start)
        invertir = False
        if dist_ends>dist_end_start:
            if dist_start_end < dist_starts:
                invertir = True
        if dist_starts<dist_start_end:
            if dist_end_start < dist_ends:
                invertir = True
        if invertir:
            aligned_bundle.append(fiber[::-1])
        else:
            aligned_bundle.append(fiber)
    
    return aligned_bundle

def euclidean_distance(p1, p2):
    return np.linalg.norm(np.array(p1) - np.array(p2))

def trajectory_distance(fiber1, fiber2):
    return sum(euclidean_distance(p1, p2) for p1, p2 in zip(fiber1, fiber2))

def is_fiber_in_correct_order(fiber, reference_fiber):
    direct_distance = trajectory_distance(fiber, reference_fiber)
    reversed_distance = trajectory_distance(fiber[::-1], reference_fiber)
    return direct_distance < reversed_distance

def align_bundle_2(fibers,reference_fiber):
    aligned_fibers = []
    for fiber in fibers:
        if is_fiber_in_correct_order(fiber, reference_fiber):
            aligned_fibers.append(fiber)
        else:
            aligned_fibers.append(fiber[::-1])
    return aligned_fibers

def reformat_string(input_string):
    # Buscar patrón con un número que se repite
    match = re.match(r'([a-z]+)_([A-Z]+)(\d+)_([A-Z]+)\3', input_string)
    
    if match:
        hemisferio = match.group(1)
        palabra1 = match.group(2)
        numero = match.group(3)
        palabra2 = match.group(4)
        
        # Formatear la nueva cadena
        new_string = f"{hemisferio}_{palabra1}-{palabra2}_{numero}"
        return new_string
    else:
        return input_string

def align_bundle(bundle, ref):
    aligned_bundle = []
    for fiber in bundle:
        d = dist(EuclideanMetric(), fiber, ref)/21
        f = mdf(fiber, ref)
        if f < d:
            aligned_bundle.append(fiber[::-1])
            
        else:
            aligned_bundle.append(fiber)

    return aligned_bundle

"""opens a directory and returns a list of all the files in it"""
def get_files(path):
    files = []
    for file in os.listdir(path):
        if file.endswith(".bundles"):
            files.append(file)
    return files

path = sys.argv[1] #bundles path
bundles = get_files(path)

#Creates aligned fibers folder
aligned_path = sys.argv[2]
if not os.path.exists(aligned_path):
    os.makedirs(aligned_path)
    os.makedirs(aligned_path+'/left')
    os.makedirs(aligned_path+'/right')

centroids_path = os.getcwd()+"/centroids"
centroides = set()
for centroid in os.listdir(centroids_path):
    centroides.add(centroid.split(".")[0])
n = 0
for bundle in bundles:
    hemi = '/right'
    if 'LEFT' in bundle or 'lh' in bundle:
        hemi = '/left'
    start = bundle.find("lh") #assumes that bundle name has lh or rh in its file name
    if start == -1:
        start = bundle.find("rh")
    name = "_".join(bundle[start:].split("_")) #saves the name of the regions of the bundle e.g "AR_ANT.bundles"
    bund = bt.read_bundle(path + bundle)
    #Get Atlas Bundle to get centroid reference
    atlas_path = (os.path.abspath(os.path.join(os.getcwd(),'../4-segmentation/atlasRo/bundles/'))).replace("\\","/")
    atlas_bundle = bt.read_bundle(atlas_path+"/MNI_atlas_"+name)
    name1 = reformat_string(name.replace(".bundles",""))
    if name1 in centroides:
        centroid = bt.read_bundle(centroids_path+"/"+name1+".bundles")[0]
        n += 1
    else:
        centroid = np.mean(atlas_bundle, axis = 0)
    aligned_bundles = align_bundle(bund, centroid)
    bt.write_bundle(aligned_path + hemi + '/aligned_' +  name, aligned_bundles)
print(n)