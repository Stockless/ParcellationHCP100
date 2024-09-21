from operator import itemgetter
import bundleTools3 as BT3
import sys
import random
import os
import re
"""Makes bundles from a fiber index list obtained from a segmented large fascicle (lh_AR, rh_CG2, etc.)
fiber list example:
[1 4 7 2 9 5]
"""
def get_bundle(idx,bundle):
    return list(map(bundle.__getitem__,idx))

fib_dir = sys.argv[1] #directory of fibers lists files
bundles_dir = sys.argv[2]
bundles_path = "bundles/"
bundles_files = os.listdir(bundles_dir)
if not os.path.exists(bundles_path):
    os.makedirs(bundles_path)
# print(bundles_files)
for folder in os.listdir(fib_dir): #iterates over the segmented fascicles folders e.g fibers/lh_AR/
    print(folder)
    if folder+'.bundles' not in bundles_files:
        continue
    start = folder.find("lh")
    if start == -1:
        start = folder.find("rh")
    hemi = folder[start:start+2]
    for filename in os.listdir(fib_dir+'/'+folder): #iterates over the fibers list on the segmented fascicle folder e.g fibers/lh_AR/fibs_CAC_Li.txt .. fibs_PrC_ANT.txt
        print(bundles_dir+folder)
        bundle = BT3.read_bundle(bundles_dir+folder+'.bundles') #reads the bundle of the fascicle to segment (must have the same name of the folder with the segmented fascicle)
        fibs = open(fib_dir+'/'+folder+'/'+filename,'r')
        fibs_idx = fibs.readline().split()
        idx = list(map(int,fibs_idx))
        name = "_".join(re.split('[_.-]',filename)[-3:-1])
        if len(idx) == 0: #skips lists with no fibers
            continue
        bun = get_bundle(idx,bundle)
        print(type(bun),len(bun))
        BT3.write_bundle(bundles_path+hemi+"_"+name+'.bundles',bun)
