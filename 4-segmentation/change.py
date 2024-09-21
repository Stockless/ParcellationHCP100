
import os
import sys

# opens a directory and changes the name of all files in it
def change_name(dir_name):
    for file in os.listdir(dir_name):
        file_name = file.split('_')
        aux = file_name[-1]
        file_name.pop()
        file_name = file_name + aux.split('.')
        #print(len(file_name))
        if len(file_name) == 7:
            os.rename(dir_name + '/' + file, dir_name + '/' + file_name[0]+"_"+file_name[1]+"_"+file_name[2]+"_"+file_name[3]+file_name[5]+"_"+file_name[4]+file_name[5]+"."+file_name[6])
            #print(file_name[0]+"_"+file_name[1]+"_"+file_name[2]+"_"+file_name[3]+file_name[5]+"_"+file_name[4]+file_name[5]+"."+file_name[6])
        if len(file_name) == 6:
            #os.rename(dir_name + '/' + file, dir_name + '/' + file_name[0]+"_"+file_name[1]+file_name[2]+"."+file_name[3])
            print(file_name)

change_name("AtlasRo/bundles")