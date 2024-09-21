import numpy as np
import bundleTools as bt
from dipy.data import get_fnames
from dipy.segment.clustering import QuickBundles
from dipy.viz import actor, colormap, window

input_path = "D:/documentos/universidad/TESIS/parcellation/parcellation-master/10-bundle-clustering/atlasMNI/bundles/MNI_atlas_rh_AR_ANT.bundles"
output_path = "D:/documentos/universidad/TESIS/parcellation/parcellation-master/10-bundle-clustering/segmented-bundles/MNI_atlas_rh_AR_ANT"
bundle = bt.read_bundle(input_path)
qb = QuickBundles(threshold=20)
clusters = qb.cluster(bundle)
interactive = True
print("Cluster sizes:", len(clusters[0]), len(clusters[1]))


colormap = colormap.create_colormap(np.arange(len(clusters)))

scene = window.Scene()
scene.SetBackground(1, 1, 1)
scene.add(actor.streamtube(clusters.centroids, colormap, linewidth=0.4))
window.record(scene, out_path="centroids.png", size=(600, 600))
if interactive:
    window.show(scene)
j=0
for i in clusters:
    bt.write_bundle(output_path+str(j)+".bundles",i)
    j += 1