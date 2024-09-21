import numpy as np
import random as rd
import math
import read_write_bundle as bt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from dipy.tracking.utils import length
from subprocess import call
import time
import os
from operator import itemgetter 
import seaborn as sns
from scipy.spatial import distance
import time
from joblib import Parallel, delayed
import multiprocessing
import matplotlib
import seaborn as sns
import numpy.matlib
import random
import pickle
import nibabel.streamlines.tck as TF
from nibabel.streamlines.tractogram import Tractogram
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import nibabel as nb
import seaborn as sns
import itertools
from numpy.linalg import inv
from scipy.spatial import cKDTree

num_cores = multiprocessing.cpu_count()
print(num_cores)

T= inv(np.load('tck2bundles.npy'))
def apply_aff_point(inPoint,T):
  Tfrm = T#N.array([[0.6, 0, 0, 0],[0, -0.6, 0, 5],[ 0, 0, -0.6, 0],[0, 0, 0, 1]])
  tmp = Tfrm * np.transpose(np.matrix(np.append(inPoint,1)))
  outpoint = np.squeeze(np.asarray(tmp))[0:3]
  return outpoint

def apply_aff_bundle(bunIn,bunOut):
  points=bt.read_bundle(bunIn)
  newPoints=[]
  for fib in points:
    newfib=[]
    for p in fib:
      pt=apply_aff_point(p,T)
      newfib.append((pt))
    newPoints.append(np.asarray(newfib,dtype=np.float32))
  #bt.write_bundle(bunOut, newPoints)
  return newPoints

def savetck(fibs,outname):
    centroids_tractogram_file = Tractogram(streamlines = fibs)
    centroids_tractogram_file.affine_to_rasmm = np.eye(4)
    centroids_tck = TF.TckFile(centroids_tractogram_file,header = {'timestamp':0})
    centroids_tck.save(outname)

def readtck(tractogram_file):
    print('Reading streamlines')
    fibra=TF.TckFile(tractogram_file)
    fibTCK=fibra.load(tractogram_file).streamlines
    #print(fibTCK)
    fibs=[]
    print('Iterating over streamlines')
    for i in range(len(fibTCK)):
        #print(i)
        fib=fibTCK[i]
        fibs.append(fib)
    return fibs

def compute_dice_voxel(density_1, density_2):
    """
    Compute the overlap (dice coefficient) between two
    density maps (or binary).
    Parameters
    ----------
    density_1: ndarray
        Density (or binary) map computed from the first bundle
    density_2: ndarray
        Density (or binary) map computed from the second bundle
    Returns
    -------
    A tuple containing
        float: Value between 0 and 1 that represent the spatial aggrement
            between both bundles.
        float: Value between 0 and 1 that represent the spatial aggrement
            between both bundles, weighted by streamlines density.
    """
    overlap_idx = np.nonzero(density_1 * density_2)
    numerator = 2 * len(overlap_idx[0])
    denominator = np.count_nonzero(density_1) + np.count_nonzero(density_2)

    if denominator > 0:
        dice = numerator / float(denominator)
    else:
        dice = np.nan

    overlap_1 = density_1[overlap_idx]
    overlap_2 = density_2[overlap_idx]
    w_dice = np.sum(overlap_1) + np.sum(overlap_2)
    denominator = np.sum(density_1) + np.sum(density_2)
    if denominator > 0:
        w_dice /= denominator
    else:
        w_dice = np.nan

    return dice


if not os.path.exists('AtlasRo_tck'):
    os.mkdir('AtlasRo_tck')
if not os.path.exists('AtlasRo_images'):
    os.mkdir('AtlasRo_images')

bundles = [f for f in os.listdir('AtlasRo') if f.endswith('.bundles')] #Se listan bundles del atlas
atlas_allfibers = [] #Aqui se iran guardando todas las fibras de los bundles del atlas
for idx,bun in enumerate(bundles): #Se itera en los bundles
    print(idx,bun)
    fibs_aff = apply_aff_bundle('/media/stocktan/DATA/documentos/universidad/TESIS/parcellation/parcellation-master/filtrarAtlasRo/centroids/'+bun,'') #Se aplica una transformacion afin para poder visualizar las fibras en MNI y en MRtrix (mrview)
    savetck(fibs_aff,'AtlasRo_tck/'+bun.split('.')[0]+'.tck') #Se guardan los bundles en formato .tck
    call(["tckresample AtlasRo_tck/"+bun.split('.')[0]+".tck AtlasRo_tck/"+bun.split('.')[0]+"_up1mm.tck -step_size 1 -force"],shell=True) #Se upsamplean las fibras para que tengan puntos cada 1mm
    call(["tckmap -contrast tdi AtlasRo_tck/"+bun.split('.')[0]+"_up1mm.tck AtlasRo_images/"+bun.split('.')[0]+"_up1mm.nii.gz -vox 1 -map_zero -force"],shell=True) #Se genera un mapa de densidad del bundle con voxel de 1x1x1 mm
    call(["mrthreshold AtlasRo_images/"+bun.split('.')[0]+"_up1mm.nii.gz AtlasRo_images/"+bun.split('.')[0]+"_mask.nii.gz -abs 1 -force"],shell=True) #Se genera una mascara binaria del bundle luego de aplicar umbralizacion al mapa de densidad
    atlas_allfibers+=fibs_aff #Se guardan las fibras 

savetck(atlas_allfibers,'atlas_allfibers.tck')  #Se guardan todas las fibras en .tck
call(["tckresample atlas_allfibers.tck atlas_allfibers_up1mm.tck -step_size 1 -force"],shell=True) #Se upsamplean las fibras con puntos cada 1mm
call(["tckmap -contrast tdi atlas_allfibers_up1mm.tck atlas_allfibers.nii.gz -vox 1 -map_zero -template MNI152_T1_1mm.nii.gz -force"],shell=True) #Se genera un mapa de densidad de todas las fibras
call(["mrthreshold atlas_allfibers.nii.gz atlas_allfibers_mask.nii.gz -abs 1 -force"],shell=True) #Se generan mascara binaria de todas las fibras


atlas_allfibers = readtck('atlas_allfibers.tck')
img_allfibers = nb.load('atlas_allfibers_mask.nii.gz').get_fdata()
vol_allfibers = np.where(img_allfibers>=1) 

volumes_bundles = []
bun_vol = {}  #Diccionario que guarda el volumen de cada bundle
for idx,bun in enumerate(bundles):
    #print(idx,bun)
    img_bun = nb.load("AtlasRo_images/"+bun.split('.')[0]+"_mask.nii.gz").get_fdata() #Se carga la mascara del bundle
    vol_bun = np.where(img_bun>=1) #Se determina el numero voxeles con intensidad 1
    volumes_bundles.append(len(vol_bun[0])) #Se guarda el volumen del bundle volumen=nvoxels x voxel_size, aqui voxel_size es 1mm^3
    bun_vol[bun] = len(vol_bun[0])


#Aqui se seleccionan los bundles con mayor volumen
bundles_filtered = [] #aqui se guarda la lista de bundles seleccionados
cutoff = np.percentile(volumes_bundles,50) #cutoff de volumen, se selecciona el 50% de bundles con mayor volumen
print("cutoff",cutoff)
for bun in bundles:
  if bun_vol[bun]>cutoff:
    bundles_filtered.append(bun.replace(".bundles",""))

print("len filt bund",len(bundles_filtered))

f = open("/media/stocktan/DATA/documentos/universidad/TESIS/parcellation/parcellation-master/atlasfiltro.txt",'a+')
for i in bundles_filtered:
    f.write(i+"\n")
f.close()