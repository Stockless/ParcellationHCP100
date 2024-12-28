
# Copyright (C) 2019  Andrea Vázquez Varela

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#Authors:
# Narciso López López
# Andrea Vázquez Varela
#Creation date: 31/01/2018
#Last update: 20/07/2019


import operator
import numpy as np
import collections
import re

def find_label(name,parcel_names):
    for i,n in enumerate(parcel_names):
        if n == name:
            return i
    return -1

def most_common(lst):
    return max(set(lst), key=lst.count)

def most_probable(map):
    if (len(map)!=0):
        return max(map.items(), key=operator.itemgetter(1))[0]
    else:
        return -1

def labels_by_prob(map,prob):
    labels = []
    for key,value in map.items():
        if (value < prob):
            labels.append(key)
    return labels

def fusion_names(n1,n2):
    index = n1.split("_")[1]
    name1_split = n1.split("_")[0].split("-")
    name2_split = n2.split("_")[0].split("-")
    new_name = ""
    for i,name1 in enumerate(name1_split):
        if i == 0:
            new_name += name1
        else:
            new_name+="-"+ name1
    for i,name2 in enumerate(name2_split):
        if i > 0:
            new_name+="-"+name2
    new_name+="_"+index
    return new_name

def fusion_names_2(n1,n2):
    n1_splitted = re.split("[_-]",n1)
    n2_splitted = re.split("[_-]",n2)
    new_name = ""
    for i,name1 in enumerate(n1_splitted):
        name1 = "" if name1.isdigit() else name1
        if name1 in new_name:
            continue

        if i == 0:
            new_name += name1
        else:
            new_name+="-"+ name1
    for i,name2 in enumerate(n2_splitted):
        name2 = "" if name2.isdigit() else name2
        if name2 in new_name:
            continue
        new_name+="-"+name2

    index = "_"+n1_splitted if n1_splitted[-1].isdigit() else ""
    new_name += index
    return new_name

def intersection(lst1, lst2):
    return set(lst1) & set(lst2)

def union(lst1, lst2):
    final_list = list(set(lst1) | set(lst2))
    return final_list

def get_inter_coords(coords,index):
    return np.asarray([coords[3*index], coords[3*index+1], coords[3*index+2]])

def sort_graphs(dic):
    dic = collections.OrderedDict(dic)
    return dict(sorted(dic.items(), key=lambda kv: len(kv[1])))

def map_triangles(triangles):
    map_tri = {}
    for tri in triangles:
        alabel = tri.label_parcel
        map_tri[alabel] = [tri]
    return map_tri

def get_bundle_names(bundle):
    start = bundle.find("lh") if "lh" in bundle else bundle.find("rh")
    if start == -1:
        return []  # No valid starting region found
    
    splitted = re.split("[._-]", bundle[start+3:])
    if splitted:
        splitted.pop(-1)  # Remove the last segment
    
    # Define and clean initial region
    dwm_labels = ["AR", "ANT", "POST", "CG", "2", "3", "IFO", "IL", "UN"]
    init_region = splitted[0]
    for label in dwm_labels:
        if init_region.startswith(label):
            init_region = init_region[len(label):]
    splitted[0] = init_region
    init_region = re.sub(r"\d+$", "", init_region)
    
    # Define and clean ending region if present
    end_region = re.sub(r"\d+$", "", splitted[1]) if len(splitted) > 1 else None

    return splitted, [init_region, end_region] if end_region else [init_region]