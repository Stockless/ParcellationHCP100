import bundleTools as BT
import os
import sys
import visual_tools as vt
from collections import defaultdict
import vtk
import numpy as np

#For loading more than one bundle
def load_bundles(path):
    dir = os.listdir(path);
    bundles = []
    for bundle in dir:
        if bundle.endswith('.bundles'):
            bundles.append(BT.read_bundle("bundles/"+bundle))
    return bundles

def read_intersection( infile ):

    f = open(infile, 'r');

    total_triangles = np.uint32(f.readline());
    InTri = list(map(np.uint32,f.readline().split()));
    print(total_triangles,len(InTri))
    FnTri = list(map(np.uint32,f.readline().split()));

    InPoints = list(map(np.float32,f.readline().split()))
    FnPoints = list(map(np.float32,f.readline().split()))

    fib_index = list(map(np.uint32,f.readline().split()))

    f.close();
    return InTri, FnTri, InPoints, FnPoints, fib_index;

#Loads the mesh
Lhemi_path = sys.argv[1] ; # 'lh.obj' or 'rh.obj'
Intersection_path = sys.argv[2]; # e.g: 'aligned_lh_AR_ANT.txt'
Lvertex, Lpolygons = BT.read_mesh_obj(Lhemi_path)
Lhemi = vt.Polygon(Lvertex, Lpolygons);
Lhemi.setOpacity(0.6)

InTri, FnTri, InPoints, FnPoints, fib_index = read_intersection(Intersection_path)

in_tri = vt.Polygon(Lvertex, Lpolygons[InTri])
in_tri.setColor((1.0,1.0/4,1.0/4))

fn_tri = vt.Polygon(Lvertex, Lpolygons[FnTri])
fn_tri.setColor((0.0,0.0,1.0))
# bundles = load_bundles('bundles/')

bundle_path = sys.argv[3]
bun = BT.read_bundle(bundle_path) #loads one bundle
vt.visual_allpoints(bun,Lhemi,in_tri,fn_tri)
