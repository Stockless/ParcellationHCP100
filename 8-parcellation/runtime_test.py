import os
import subprocess as sp

dc_list = [0.15]
idc_list = [0.10]
size_thr_list = [0.1]
ero_list = [1]
dil_list = [6]

output = 'Parcellation_test/'
if not os.path.exists(output):
    os.mkdir(output)
for i in range(30):
    for dc in dc_list:
        for idc in idc_list:
            for size_thr in size_thr_list:
                for dil in dil_list:
                    for ero in ero_list:
                        print('Running:\n'+'dc='+str(dc)+'\nidc='+str(idc)+'\nthr='+str(size_thr)+'\ndil='+str(dil)+'\nero='+str(ero))
                        output = 'Parcellation_test/output_dc'+str(dc)+'_idc'+str(idc)+'_thr'+str(size_thr)+'_dil'+str(dil)+'_ero'+str(ero)+'/'
                        if not os.path.exists(output):
                            os.mkdir(output)
                        sp.call(['python3', 'parcellation.py', '--output-dir', output, '--Lvlabels-file', 'inputs/lh_labels.txt', '--Rvlabels-file', 'inputs/rh_labels.txt', '--Intersection-dir', 'inputs/intersection', '--parcel-names', 'inputs/dk_names.txt', '--LVtk-file', 'inputs/lh.obj', '--RVtk-file', 'inputs/rh.obj', '--ero', str(ero), '--dil', str(dil), '--dc', str(dc), '--idc', str(idc), '--size-thr', str(size_thr)])


