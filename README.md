## From Coarse to Fine-Grained Parcellation Refactored
This repo contains a refactoring of the work developed in the paper [From Coarse to Fine-Grained Parcellation of the Cortical Surface Using a Fiber-Bundle Atlas](https://www.frontiersin.org/articles/10.3389/fninf.2020.00032/full).

For further details, please refer to the original paper.
![fs](/images/parcellation.png)


## Usage

The code is developed to work with tractography data from the ARCHI database. The data can be downloaded from [here](https://www.frontiersin.org/articles/10.3389/fninf.2020.00032/full).




#### ************************** Pipeline **********************************

1. Merge tractography files.
  - Go to 1-merge folder and run the following command:
    `python merge.py arg1`
    where arg1 is the path to the folder containing the subjects dir.
  - The tractography files are merged into one .bundles file and one .bundlesdata file.

2. Resample tractography files to 21 points.
  - Go to 2-resample folder and run the following command:
  `python multiple_resample.py arg1`
  where arg1 is the path to the folder containing the subjects dir.
  - The tractography file of each subject is resampled into a new bundle inside a new folder named 'resampled' inside each subject folder.
  - For further details of this process go to the readme inside 2-resample folder.
3. Transform tractography files to Talairach space.
  - Go to 3-transforms and run the following command:
  `python multiple_transform.py arg1 arg2 arg3 arg4`
  The args are the following:
    - arg1: subjects dir
    - arg2: tractography folder to be transformed
    - arg3: output folder for the new transformed tractography
    - arg4: input matrix file
  - The tractography of each subject is transformed into a new bundle inside at the output folder given in arg3.
  - For further details of this process go to the readme inside 3-transforms folder.
4. Segment the tractography files.
  - Go to 4-segmentation and run the following code:
  `python multiple_segmentation.py arg1 arg2 arg3`
  The args are the following:
    - arg1: subjects dir
    - arg2: folder containing the atlas bundles
    - arg3: file with the atlas information
  - The tractography file of each subject is segmented in several files on a new folder named `segmented_Tal`.
  - IMPORTANT: make sure that the tractography bundle of each subject is at the same space of the atlas (for this use case we are using the Talairach space for the atlas).
  - For further details of this process go to the readme inside 4-segmentation folder.

5. Obtain inverse affine transform matrix from Talairach to T1 space.
  - Go to 3-transforms and run the following code:
  `python inversaafintrm.py arg1 arg2 arg3`
  The args are the following:
    - arg1: subjects dir
    - arg2: input file name of the transform matrix
    - arg3: output file name of the inverse affine transform matrix
  - The new affine transformation matrix will be created inside the `TransformMatrices` folder inside each subject folder.
  - For this use case, arg2 should be `trmT1ToTal.trm` (transform matrix from T1 space to Talairach space) and arg3 could be named `Tal_to_T1.trm`.

6. Transform the segmented bundles to the cortical meshes space.
  - Go to 3-transforms and run the following command:
  `python transform_T1.py arg1`
  The args are the following:
    - arg1: subjects dir
  - The tractography of each subject is transformed into a new bundle inside at the output folder given in arg3.
  - For this use case we are transforming to T1 space. For other space transformation, we can use the multiple_transform.py script giving the required arguments for other type of trasnformation (further code/script modifications should be done to other pipeline codes if using other transformation).

7. Align fibers of the segmented tractographies.
  - Go to 5-alignment and run the following code:
  `python multiple_align.py arg1`
  The args are the following:
    - arg1: subjects dir
  - The segmented bundles of each subject will have its fibers aligned on new segmented bundles inside a new folder named `aligned_T1`.

8. Intersect the segmented tractographies with the cortical mesh.
  - Go to 6-intersection and run the following code:
  `python interx.py arg1 arg2`
  The args are the following:
    - arg1: subjects dir folder
    - arg2: meshes dir folder
  - Intersections files will be created containing the triangles of fibers intersecting the cortical mesh of each subject inside `intersection` folder inside the root of 6-intersection folder.
  - For further details of this process go to the readme inside 6-intersection folder.

9. Filter the intersection files.
  ### Filter DWM fascicles
  As this fascicles aren't included in Desikan-Killiany regions, we have to segment them with its triangles intersection files.
  - Go to 7-filter and run the following code:
  `python segment_large_fascicles.py arg1 arg2`
  The args are the following:
    - arg1: meshes dir folder
    - arg2: intersections folder
    - arg2: folder containing the labels of each hemisphere
    - arg3: file containing the name of the regions of the atlas (regions of Desikan-Killiany for this use case)
  - This will delete the DWM intersections files (like AR, CG, UN, etc.) and create new intersection files with the connecting regions of the atlas for this fascicles.

  ### Filter incorrect segmented fibers
  Fibers of segmented fascicles can be pointing at wrong regions inside the fascicle. For this we delete the fibers that don't belong to the connecting regions of the fascicles (given by the name of the fascicle).
  - Go to 7-filter and run the following code:
  `python filter_intersection.py arg1 arg2`
  The args are the following:
    - arg1: meshes dir folder
    - arg2: intersections folder
    - arg2: folder containing the labels of each hemisphere
    - arg3: file containing the name of the regions of the atlas (regions of Desikan-Killiany for this use case)
  - This will delete the wrong segmented fibers of fascicles and create new filtered intersection files inside `filtered_intersections` at the root of 7-filter folder.

9. Parcellation:
  - Go to 8-parcellation and run the following code:
  python parcellation.py --output-dir arg1 --Lvlabels-file arg2 --Rvlabels-file arg3 --Intersection-dir arg4 --parcel-names arg5 --LVtk-file arg6 --RVtk-file arg7

  Additional arguments could be added:
  - --traceability: For traceability of the process (like save info of preliminar subparcels information).
  Values:
    - y
    - n
  - --size_thr: threshold for deleting small parcels. Value should be a float number bigger than 0
  - 0.1 to 0.4 (could be more but higher number can delete bigger parcels)
  - --dc_thr: threshold for density center. Value should be a float number bigger than 0
  - --idc: threshold for density center intersection. Value should be a float number bigger than 0
  - --ero: erotion factor. Value should be a integer number starting from 0 (for no erotion) or 1 (for a 1 factor of erotion)
  - --dil: dilation factor. Value should be a integer number starting from 0 (for no erotion) or 1 (for a 1 factor of erotion)


