import os
import random
import numpy as np
import sys
#import utils.bundleTools as bt
import bundleTools3 as bt3
# import utils.visualizationTools as vt
from dipy.segment.metric import mdf
from dipy.segment.metric import dist
from dipy.segment.metric import EuclideanMetric

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

def calculate_centroid(bundle):
    return np.mean(bundle, axis = 0)

"""opens a directory and returns a list of all the .bundles files in it"""
def get_files(path):
    files = []
    for file in os.listdir(path):
        if file.endswith(".bundles"):
            files.append(file)
    return files

path = sys.argv[1]
bundles = get_files(path)

for bundle in bundles:
    name = "_".join(bundle.split("_")[2:])
    bundle = bt3.read_bundle(path + bundle)
    centroid = calculate_centroid(bundle)

    aligned_bundles = align_bundle(bundle, centroid)

    bt3.write_bundle(path + 'aligned_' + name, aligned_bundles)