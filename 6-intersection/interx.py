# -*- coding: utf-8 -*-

import os;
import numpy as np;
import subprocess as sp;
from time import time;
import sys

#========================== Cálculo de intersección ==========================================    

if len(sys.argv) < 3:
    print("Usage: python interx.py ../subs_dir/ ../subs_meshes_dir/")
    exit(1)
subs_path = sys.argv[1]
meshes_path = sys.argv[2]
sp.call(['make'])
for sub in os.listdir(subs_path):
    t1 = time()
    Lhemi_path = meshes_path+"/"+sub+'/lh.obj' # left mesh path
    Rhemi_path = meshes_path+"/"+sub+'/rh.obj' # right mesh path
    
    Lbundles_path = subs_path+"/"+sub+'/aligned4/left/' # bundles path
    Rbundles_path = subs_path+"/"+sub+'/aligned4/right/' # bundles path

    if not os.path.exists(Lbundles_path) or not os.path.exists(Rbundles_path):
        print("No aligned bundles found. Align the bundles with the bundle alignment step and try again.")
        sys.exit(1)
    
    Lintersection_path = os.getcwd() + '/intersection/' + sub + '/' + sub + '_left-hemi/' # left intersection path
    Rintersection_path = os.getcwd() + '/intersection/' + sub + '/' + sub + '_right-hemi/' # right intersection path
    
    if not os.path.exists(os.getcwd() + '/intersection/'):
        os.mkdir(os.getcwd() + '/intersection/')
    
    if not os.path.exists(os.getcwd() + '/intersection/'+ sub):
        os.mkdir(os.getcwd() + '/intersection/'+sub)
        
    sp.call(['./interx', Lhemi_path, Rhemi_path, Lbundles_path, Rbundles_path, Lintersection_path, Rintersection_path])
    
    print('Tiempo de ejecución Ldirect: ' + str(time()-t1) + '[s]')
    
    del Lbundles_path, Rbundles_path