#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Agosto 10 20:09:15 2020

@author: fondecyt-1190701
"""
import numpy as np
import os
from collections import defaultdict

import bundleTools as bt
import visualizationTools as vt
import matplotlib.pyplot as plt
from scipy.stats import norm

import pickle


def load_restricted_triangles():
    #Existen ciertos triángulos que, por definición, no pueden representar parcelas. Estos corresponden a regiones posterior-inferior de la línea media...
    #...donde ambos hemisferios se unen, ya que en estricto rigor no corresponden a corteza, si no que a materia blanca/gris interna.
    
    #La parcelación obtenida por Lefranc indica cuales triángulos están "restringidos", los cuales están guardados en el archivo .pkl:
    #Lrestricted.pkl para el hemisferio izquierdo.
    #Rrestricted.pkl para el hemisferio derecho.
    
    with open('pickles/Lnope.pkl', 'rb') as f:
        Lrestricted = set(pickle.load(f))
    
    with open('pickles/Rnope.pkl', 'rb') as f:
        Rrestricted = set(pickle.load(f))
    
    return Lrestricted, Rrestricted

#Carga parcelas ('hard','cc','final')    
def load_parcels (parcels, parcels_path): #hard cc final
    Lparcels_path =  parcels_path + '/'+ parcels +'_parcels/left/'
    Rparcels_path = parcels_path + '/'+ parcels +'_parcels/right/'

    Lparcels = dict()
    Rparcels = dict()

    for Dir in os.listdir(Lparcels_path):
        Lparcels[Dir.split('.')[0]] = bt.read_parcels(Lparcels_path + Dir)
        
    for Dir in os.listdir(Rparcels_path):
        Rparcels[Dir.split('.')[0]] = bt.read_parcels(Rparcels_path + Dir)
            
    return Lparcels, Rparcels

#Formula del coeficiente DICE
def DSC(A,B):
    set_A = set(A)
    set_B = set(B)
    
    interx = set_A.intersection(set_B)
    if (len(set_A) + len(set_B)) == 0:
        return 0
    else:
        return (2*len(interx))/((len(set_A)+len(set_B)))

#Calculo de coeficiente DICE para comparar parcelas con atlas
def dise_comparision (atlas_comparision, parcels_path, dice_thr):
    # Lectura de parcelas propias.    
    lh_mine_dict, rh_mine_dict = load_parcels('final', parcels_path)
        
    # Cargar parcelación del atlas seleccionado.
    with open(atlas_comparision + '_Lparcels.pkl', 'rb') as f:   
        lh_atlas_dict = pickle.load(f)
    
    with open(atlas_comparision + '_Rparcels.pkl', 'rb') as f:   
        rh_atlas_dict = pickle.load(f)    
        
    # Cálculo de DSC entre parcelas.
    # Diccionario parcela_propia:parcela_atlas, para aquellos casos que superan el umbral de DSC.
    lh_dice_dict = defaultdict(int)
    rh_dice_dict = defaultdict(int)

    # Lista que almacena las parcelas ya utilizadas.
    used = []
    nl_parcels = len(lh_mine_dict.keys())
    nr_parcels = len(rh_mine_dict.keys())
    print('n° de parcelas izquierdo: %d' % nl_parcels)
    print('n° de parcelas derecho: %d' % nr_parcels)
    # Iteración sobre todas las parcelas propias obtenidas.
    suma_dice = 0
    len_dice = 0
    for parcel, tris in lh_mine_dict.items():
        #Se registra el mayor DSC obtenido entre todas las comparaciones, y el nombre de la parcela asociada.
        dice_max = 0 
        winner = ''
      
        #Se itera sobre las parcelas del atlas.
        for ref, rtris in lh_atlas_dict.items():
            
            #Si la parcela ya fue utilizada, ignorar.
            if ref in used:
                continue
            
            #Cálculo de DSC entre parcela propia y parcela del atlas.
            dice = DSC(tris,rtris)
            
            #Si es mayor al máximo obtenido, actualizar.
            if dice > dice_max:
                dice_max = dice
                winner = ref
        #print(dice_max)
        suma_dice += dice_max
        if dice_max != 0:
            len_dice += 1    
        #Si el mayor valor obtenido supera el umbral de similitud, se agrega al diccionario de parcelas similares, y a la lista de parcelas ya utilizadas.
        if dice_max >= dice_thr:
            lh_dice_dict[parcel] = winner
            used.append(winner)
            
        # Si no supera el umbral, continuar con la siguiente parcela.
        else:
            continue
    
    #Se repite lo mismo para el hemisferio derecho.
    used = []
    
    for parcel, tris in rh_mine_dict.items():
        dice_max = 0
        winner = ''
      
        for ref, rtris in rh_atlas_dict.items():
            if ref in used:
                continue
            dice = DSC(tris,rtris)
            if dice > dice_max:
                dice_max = dice
                winner = ref
        suma_dice += dice_max
        if suma_dice != 0:
            len_dice += 1    
        if dice_max >= dice_thr:
            rh_dice_dict[parcel] = winner
            used.append(winner)
            
        else:
            continue  
    
    #Diccionario de parcelas similares, considerando los triángulos asociados en el caso de la parcelación propia.
    lh_mine_common = {k:lh_mine_dict[k] for k in list(lh_dice_dict.keys()) if k in lh_mine_dict}
    rh_mine_common = {k:rh_mine_dict[k] for k in list(rh_dice_dict.keys()) if k in rh_mine_dict}
    
    #Diccionario de parcelas similares, considerando los triángulos asociados en el caso de la parcelación del atlas.
    lh_atlas_common = {k:lh_atlas_dict[k] for k in list(lh_dice_dict.values()) if k in lh_atlas_dict}
    rh_atlas_common = {k:rh_atlas_dict[k] for k in list(rh_dice_dict.values()) if k in rh_atlas_dict}
    
    if not os.path.exists(parcels_path+'Dise/'):
        os.makedirs(parcels_path+'Dise/')
        
    save_txt(parcels_path+'Dise/lh_list_'+atlas_comparision.replace("/",""),lh_mine_common, dice_thr)
    save_txt(parcels_path+'Dise/rh_list_'+atlas_comparision.replace("/",""), rh_mine_common, dice_thr)
    #print(suma_dice/len_dice)
    return lh_mine_common, rh_mine_common, lh_atlas_common, rh_atlas_common

def save_txt (name, diccionary, dice_thr):
    count=0
    with open(name+'_'+str(int(dice_thr*100))+'.txt', 'w') as f:
        for line in diccionary.keys():
            f.write(line)
            f.write('\n')
            count+=1
        f.write('\nTotal=%d' % count)
   #%%     

def visualize_parcellation(meshes_path, L_sp, R_sp, sub, seed = False):
    #Cargar trianguos restringidos
    Lrestricted, Rrestricted= load_restricted_triangles()
    
    #Semilla para colores aleatorios
    if seed != False:
        np.random.seed(seed)
    
    #Parcelas finales a graficar
    final_parcels = set()

    for k in L_sp.keys():
        final_parcels.add(k)
        
    for k in R_sp.keys():
        final_parcels.add(k)

    fp = list(final_parcels)    

    #Paleta de colores según cantidad de parcelas
    paleta = [(np.random.random(), np.random.random(), np.random.random()) for i in range(len(fp))]

    #Directorios de los mallados corticales
    Lhemi_path = meshes_path + sub + '/lh.obj'; # left hemisphere path
    Rhemi_path = meshes_path + sub + '/rh.obj'; # right hemisphere path

    #Lectura de mallados
    Lvertex, Lpolygons = bt.read_mesh_obj(Lhemi_path)
    Rvertex, Rpolygons = bt.read_mesh_obj(Rhemi_path)
    
    Lhemi = vt.Polygon(Lvertex, Lpolygons);
    Rhemi = vt.Polygon(Rvertex, Rpolygons);
    
    Lhemi.setOpacity(1);
    Rhemi.setOpacity(1);
    
    #Creación del render a visualizar
    render = vt.Render();
    
    #Se renderizan los mallados
    render.AddActor(Lhemi);
    render.AddActor(Rhemi);
    
    #Para cada parcela del hemisferio izquierdo...
    for k, v in L_sp.items():
        #Se selecciona un color de la paleta
        color = paleta[fp.index(k)]
        
        if len(v) == 0:
            continue
      
        #De todos los triángulos con parcelas, se eliminan aquellos que pertenecen al conjunto de triángulos restringidos.
        v_restricted = list(set(v).difference(Lrestricted))

        #Se renderizan los triángulos y polígonos.
        sp_tri = vt.Polygon(Lvertex, Lpolygons[v_restricted]);
        sp_tri.setColor((color[0], color[1], color[2]));
        render.AddActor(sp_tri);

    #Ídem para el derecho
    for k, v in R_sp.items():
        color = paleta[fp.index(k)]
    
        if len(v) == 0:
            continue
    
        v_restricted = list(set(v).difference(Rrestricted))
    
        sp_tri = vt.Polygon(Rvertex, Rpolygons[v_restricted]);
        sp_tri.setColor((color[0], color[1], color[2]));
        render.AddActor(sp_tri);
    
    render.Start();
    del render


def multiple_DSC_comp():

    #Comparacion Dice
    atlas_info = {}
    for atlas in atlases:
        print(''.join(atlas[6:]))
        atlas_info[''.join(atlas[6:])] = []
        for folder in os.listdir('Parcellation/'):
            params = folder.split('_')
            dc = ''.join(params[1][2:])
            idc = ''.join(params[2][3:])
            output_parcellation = 'Parcellation/'+folder+'/'
            lh_mine_common, rh_mine_common, lh_atlas_common, rh_atlas_common = dise_comparision(atlas, output_parcellation, dice_thr)
            atlas_info[''.join(atlas[6:])].append([dc, idc, 0.1, len(lh_mine_common.keys()), len(rh_mine_common.keys()), len(lh_mine_common.keys()) + len(rh_mine_common.keys())])
    
    """Para guardar los resultados en una tabla de latex"""
    file = open("comparisons.txt", "w")
    for k,v in atlas_info.items():
        file.write('\\multirow{'+k+'}')
        for info in v:
            i = 0
            file.write(' & '+info[i]+' & '+str(info[i+1])+' & '+str(info[i+2])+' & '+str(info[i+3])+' & '+str(info[i+4])+' & '+str(info[i+5])+' \\\\ \n')
        file.write('\\hline\n')
    file.close()

#Sujeto base. Puede ser cualquiera, ya que todos los mallados tienen triángulos correspondientes.
sub = '101410'
meshes_path= '../../../HCP100/Mallados/'
dice_thr=0.6

# Selección de atlas a comparar.
atlases = ['atlas/Lefranc','atlas/Brainnetome','atlas/Narciso','atlas/Richards']
semilla_visualizacion = 47

# Comparación Dice
final_parcels = 'Parcellation_individual/'
lh_atlas_dict, rh_atlas_dict = load_parcels('final', '../8-parcellation/output/')
#lh_atlas_dict, rh_atlas_dict = load_parcels('final', 'ParcellationHARDI/')
with open('atlas/Richards_Lparcels.pkl','wb') as handle:
    pickle.dump(lh_atlas_dict,handle)
with open('atlas/Richards_Rparcels.pkl','wb') as handle:
    pickle.dump(rh_atlas_dict,handle)
#lh_mine_common, rh_mine_common, lh_atlas_common, rh_atlas_common = dise_comparision(atlases[0], final_parcels, dice_thr)
#lh_mine_common, rh_mine_common, lh_atlas_common, rh_atlas_common = dise_comparision(atlases[1], final_parcels, dice_thr)
#lh_mine_common, rh_mine_common, lh_atlas_common, rh_atlas_common = dise_comparision(atlases[2], final_parcels, dice_thr)
lh_mine_common, rh_mine_common, lh_atlas_common, rh_atlas_common = dise_comparision(atlases[3], final_parcels, dice_thr)

# Visualizacion para Dice
#visualize_parcellation(meshes_path, lh_atlas_common, rh_atlas_common, sub, seed = 47)
#visualize_parcellation(meshes_path, lh_mine_common, rh_mine_common, sub, seed = 47)

#Lparcels_ps, Rparcels_ps= load_parcels('ps', final_parcels)
#Lparcels_fp, Rparcels_fp= load_parcels('fp', final_parcels)
#Lparcels_hard, Rparcels_hard= load_parcels('hard', final_parcels)
#Lparcels_cc, Rparcels_cc= load_parcels('cc', final_parcels)
Lparcels_final, Rparcels_final= load_parcels('final', final_parcels)
print(len(Lparcels_final),len(Rparcels_final))
Lbig, Rbig = 0,0
Lsmall, Rsmall = 9999999,999999
suma_parcelas = 0 
numero_parcelas = 0
sizes = []
for v in Lparcels_final:
    if Lbig < len(Lparcels_final[v]):
        Lbig = len(Lparcels_final[v])
    if Lsmall > len(Lparcels_final[v]):
        Lsmall = len(Lparcels_final[v])
    if len(Lparcels_final[v]):
        suma_parcelas += len(Lparcels_final[v])
        numero_parcelas += 1
        sizes.append(len(Lparcels_final[v]))
print(Lbig,Lsmall,suma_parcelas/numero_parcelas)
# Plot histogram of the vector to show the distribution
plt.hist(sizes, bins=30, density=True, alpha=0.6, color='b')

# Fit a normal distribution to the data
mu, std = norm.fit(sizes)

# Plot the normal distribution (bell curve)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)

# Add labels and title
plt.title('Distribución de tamaños de parcelas')
plt.xlabel('Numero de parcelas')
plt.ylabel('Densidad de parcelas de cada tamaño')

plt.show()
suma_parcelas = 0 
numero_parcelas = 0
sizes = []
for v in Rparcels_final:
    if Rbig < len(Rparcels_final[v]):
        Rbig = len(Rparcels_final[v])
    if Rsmall > len(Rparcels_final[v]):
        Rsmall = len(Rparcels_final[v])
    if len(Rparcels_final[v]):
        suma_parcelas += len(Rparcels_final[v])
        numero_parcelas += 1
        sizes.append(len(Rparcels_final[v]))
print(Rbig,Rsmall,suma_parcelas/numero_parcelas)
# Plot histogram of the vector to show the distribution
plt.hist(sizes, bins=30, density=True, alpha=0.6, color='b')

# Fit a normal distribution to the data
mu, std = norm.fit(sizes)

# Plot the normal distribution (bell curve)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)

# Add labels and title
plt.title('Distribución de tamaños de parcelas')
plt.xlabel('Numero de parcelas')
plt.ylabel('Densidad de parcelas de cada tamaño')

plt.show()
#
#visualize_parcellation(meshes_path, Lparcels_ps, Rparcels_ps, sub, seed = semilla_visualizacion)
#visualize_parcellation(meshes_path, Lparcels_fp, Rparcels_fp, '001', seed = semilla_visualizacion)
#visualize_parcellation(meshes_path, Lparcels_hard, Rparcels_hard, '001', seed = semilla_visualizacion)
#visualize_parcellation(meshes_path, Lparcels_cc, Rparcels_cc, '001', seed = semilla_visualizacion)
visualize_parcellation(meshes_path, Lparcels_final, Rparcels_final, sub, seed = semilla_visualizacion)
