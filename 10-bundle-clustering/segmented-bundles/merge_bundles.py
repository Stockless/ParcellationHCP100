import bundleTools as BT
import bundleTools3 as BT3
all_points=[]
fiber_bun_1 = 'D:/documentos/universidad/TESIS/parcellation/parcellation-master/10-bundle-clustering/segmented-bundles/MNI_atlas_rh_AR_ANT1.bundles'
fiber_points_1 = BT.read_bundle(fiber_bun_1)
for fiber in fiber_points_1:
    all_points.append(fiber)

fiber_bun_2 = 'D:/documentos/universidad/TESIS/parcellation/parcellation-master/10-bundle-clustering/segmented-bundles/MNI_atlas_rh_AR_ANT2.bundles'
fiber_points_2 = BT.read_bundle(fiber_bun_2)
for fiber in fiber_points_2:
    all_points.append(fiber)
BT3.write_bundle('D:/documentos/universidad/TESIS/parcellation/parcellation-master/10-bundle-clustering/segmented-bundles/MNI_atlas_rh_AR_ANT_nuevo.bundles', all_points)
        