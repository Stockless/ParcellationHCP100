import subprocess as sp
import sys
import os
import shutil

if len(sys.argv) <= 3:
    print("Usage: python multiple_segmentation.py ../subjects_dir/ folder/to/atlas/bundles/ atlas/folder/atlas_info.txt")
    sys.exit(0)
subjs_dir = sys.argv[1] #Data/subs
atlas_bundles = sys.argv[2] #atlas/bundles
atlas_info = sys.argv[3] #atlas/info.txt
clean = input("Clean old non segmented files? (y/n): ").capitalize()
sp.call(['g++', '-std=c++14', '-O3', 'segmentation.cpp', '-o', 'segmentation', '-fopenmp', '-ffast-math'])
for sub in os.listdir(subjs_dir):
    print("Segmenting subject "+sub)
    non_segmented = subjs_dir+"/"+sub
    subj_bundle = subjs_dir+"/"+sub+"/resampled/resampled_"+sub+".bundles"
    output_dir = subjs_dir+"/"+sub+"/segmented6"
    sp.call(['./segmentation', '21', subj_bundle, 'subject', atlas_bundles, atlas_info, output_dir])
    if clean[0] == 'Y':
        shutil.rmtree(non_segmented)