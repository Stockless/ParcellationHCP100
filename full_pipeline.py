import subprocess as sp
import sys
subjs_dir = sys.argv[1]
meshes_dir = sys.argv[2]
atlas_bundles = 'atlas/bundles/'
atlas_info = 'atlas/atlas_info.txt'
parcel_inputs = 'inputs'
atlas_names = parcel_inputs + '/dk_names.txt'
lh_labels = parcel_inputs + '/lh_labels.txt'
rh_labels = parcel_inputs + '/rh_labels.txt'
lh_obj = parcel_inputs + '/lh.obj'
rh_obj = parcel_inputs + '/rh.obj'

sp.call('python3 multiple_resampling.py ' + subjs_dir,shell=True, cwd='2-resampling')

sp.call('python3 merge.py '+ subjs_dir,shell=True, cwd='1-merge')

sp.call('python3 multiple_transform.py ' + subjs_dir+ " resampled transformed_Tal T2_to_Tal_tr_tmp.trm",shell=True, cwd='3-transforms')

sp.call('python3 multiple_segmentation.py ' + subjs_dir + " " +atlas_bundles + " "+ atlas_info,shell=True, cwd='4-segmentation')

sp.call('python3 transform_T1.py ' + subjs_dir + " resampled transformed_Tal T2_to_Tal_tr_tmp.trm",shell=True, cwd='3-transforms')

sp.call('python3 multiple_align.py ' + subjs_dir,shell=True, cwd='5-alignment')

sp.call('python3 interx.py ' + subjs_dir + " " + meshes_dir,shell=True, cwd='6-intersection')
intersection_path = '../6-intersection/intersection/'
sp.call('python3 segment_large_fascicles.py ' + meshes_dir + ' ../6-intersection/intersection/ labels/ dk_names.txt',shell=True, cwd='7-filter')
sp.call('python3 filter_intersection.py ' + meshes_dir + ' ../6-intersection/intersection/ labels/ dk_names.txt',shell=True, cwd='7-filter')
intersection_path = '../7-filter/filtered_intersections/'
output = 'res_output/'
ero = 1
dil = 6
dc = 0.15
idc = 0.10
size_thr = 0.1
sp.call('python3 parcellation.py --output-dir '+ output+ ' --Lvlabels-file ' + lh_labels + ' --Rvlabels-file ' + rh_labels + ' --Intersection-dir ' + intersection_path +' --parcel-names ' + atlas_names + ' --LVtk-file ' + lh_obj + ' --RVtk-file '+ rh_obj,shell=True, cwd='8-parcellation')
