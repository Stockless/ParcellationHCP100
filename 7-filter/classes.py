import collections
class Fascicle:
    def __init__(self,r1,r2,tri_r1,tri_r2,prob,fiber,label_parcel):
        self.r1 = r1
        self.r2 = r2
        self.tri_r1 = tri_r1
        self.tri_r2 = tri_r2
        self.prob = prob
        self.fiber = fiber
        self.label_parcel = label_parcel

class Intersection:
    def __init__(self,numTri,InTri,FnTri,inter_in, inter_fn,id_fib):
        self.numTri = numTri
        self.InTri = InTri
        self.FnTri = FnTri
        self.inter_in = inter_in
        self.inter_fn = inter_fn
        self.id_fib = id_fib

class Vertex:
    def __init__(self,index,x,y,z,label_parcel,triangles):
        self.index = index
        self.x = x
        self.y = y
        self.z = z
        self.triangles = triangles
        self.label_parcel = label_parcel

class Triangle:
    def __init__(self,index,v1,v2,v3,label_parcel,labels_subparcel,fibers):
        self.index = index
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.label_parcel = label_parcel
        self.labels_subparcel = labels_subparcel
        self.neighbors = []
        self.prob = []
        self.prob_map =  collections.OrderedDict()
        self.fibers = fibers
        self.fibers_map = {}