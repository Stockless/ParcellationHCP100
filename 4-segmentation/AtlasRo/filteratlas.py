import os
import bundleTools as bt

# Leer los nombres del archivo de texto y guardarlos en un conjunto
def obtener_nombres(archivo_txt):
    nombres = {}
    with open(archivo_txt, 'r') as f:
        for linea in f:
            nombre = linea.split('\t')[0]
            fibras = linea.split('\t')[2].replace("\n","")
            nombres[nombre]=fibras
    return nombres

# Eliminar archivos en la carpeta que no estén en el conjunto de nombres
def eliminar_archivos_no_coincidentes(carpeta, nombres):
    n = 1
    archivos_en_carpeta = set(os.path.splitext(archivo)[0] for archivo in os.listdir(carpeta) if os.path.splitext(archivo)[1] == '.bundles')
    
    for archivo in os.listdir(carpeta):
        nombre_archivo, extension = os.path.splitext(archivo)
        if extension == '.bundles':
            bundle = bt.read_bundle(os.path.join(carpeta,archivo))
            len_bundle = len(bundle)
        if nombre_archivo not in nombres:
            print("No está en atlas_info:",n, nombre_archivo, nombres[nombre_archivo])#os.remove(os.path.join(carpeta, archivo))
            n += 1
        if nombre_archivo in nombres and str(len_bundle) != nombres[nombre_archivo]:
            print("Número distinto de fibras:",nombre_archivo,nombres[nombre_archivo],len_bundle)
    for nombre in nombres:
        if nombre not in archivos_en_carpeta:
            print("No está en la carpeta bundles: ", nombre, nombres[nombre])


# Especifica la ruta al archivo de texto y la carpeta
archivo_txt = os.getcwd() + '/atlas_info.txt'
carpeta = os.path.join(os.getcwd(),"bundles")

# Obtener los nombres y eliminar los archivos no coincidentes
nombres = obtener_nombres(archivo_txt)
eliminar_archivos_no_coincidentes(carpeta, nombres)

print("Proceso completado.")
