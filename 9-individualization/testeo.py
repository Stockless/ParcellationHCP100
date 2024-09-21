import json
from classes import *

file_path="D:/documentos/universidad/TESIS/parcellation/parcellation-master/8-parcellation/output/ps_parcels/right.json"
with open(file_path,'r') as f:
    json_string = f.read()

# Deserialize the JSON string
data = json.loads(json_string)

aparcels = []
for ap_data in data['aparcels']:
    anatomic_parcel = AnatomicParcel(ap_data['label'])
    
    for subparcel_data in ap_data['sub_parcels']:
        triangles = []
        for triangle_data in subparcel_data['triangles']:
            v1_data = triangle_data['v1']
            v1 = Vertex(v1_data['index'], v1_data['x'], v1_data['y'], v1_data['z'], v1_data['label_parcel'],[])
            
            v2_data = triangle_data['v2']
            v2 = Vertex(v2_data['index'], v2_data['x'], v2_data['y'], v2_data['z'], v2_data['label_parcel'],[])
            
            v3_data = triangle_data['v3']
            v3 = Vertex(v3_data['index'], v3_data['x'], v3_data['y'], v3_data['z'], v3_data['label_parcel'],[])
            
            triangle = Triangle(triangle_data['index'], v1, v2, v3, triangle_data['label_parcel'], triangle_data['label_subparcel'], triangle_data['fibers'])
            triangles.append(triangle)
        
        subparcel = SubParcel(subparcel_data['label'], subparcel_data['label_anatomic'], triangles, subparcel_data['inter_points'])
        
        anatomic_parcel.sub_parcels[subparcel.label] = subparcel
    
    aparcels.append(anatomic_parcel)
print(aparcels)