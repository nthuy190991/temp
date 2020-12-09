# Thesis: A Versatile and Efficient Data Fusion Methodology for Heterogeneous Airborne LiDAR and Optical Imagery Data Acquired Under Unconstrained Conditions

Matlab codes developed by TH Nguyen, Université Laval & IMT Atlantique

## Main folders:
* **phd_dev:** common folder contains the coarse registration development, and other miscellaneous codes
* **Snakes:** dedicated to the early development of snake models
* **SR:** dedicated to the super-resolution and fine registration
* **SRSM:** dedicated to the super-resolution-based snake models used in Remote Sensing MDPI article
* **R2Sonic:** dedicated to the R2Sonic suspended sediment challenges

## Registration dev: 

### Building Segment Extraction from LiDAR
- *building_region_seg:* extraction of building regions from LiDAR data, resulting a building mask
- *buildingBoundaryExtraction:* extracting boundaries given a mask (as seeding regions)

### Building Segment Extraction from Image
- *meanshift_seg_final:* run meanshift segmentation and filtering, call MeanShiftCluster
- *MeanShiftCluster:* meanshift algorithm
- *segmentFiltering*: filtering segments based on their area
- *building_img_seg_meanshift:* refinement of meanshift resulting segments with the MBR filling percentage

### Segment Matching
- *compare_segs:* compute segment MBR filling percentage (in preparation for the segment matching)
- *compare_segs_b4regis:* running GTM/RANSAC for the segment matching
- *minBoundingBox:* MBR detection algorithm
- *GTM_algo_original:* GTM algorithm. Sub-functions for GTM in the **GTM** folder

### Transformation Model Estimation
- *TME_GS_algo:* transformation model estimation, calling gold_standard_algo
- *gold_standard_algo:* Gold Standard algorithm for the estimation of the general camera matrix
- *gold_standard_algo_affine:* Gold Standard algorithm for the estimation of the affine camera matrix
- *normalize_pts:* centering of a set of points, by their center, i.e. x<sub>i</sub> = x<sub>i</sub>-mean(x<sub>i</sub>)
- *runGoldStandard:* Gold Standard algorithm resulting in different camera components K, C, R

### Fine Registration
- *regis_6params_9patches:* fine registration on patches, involving the maximization of MI/NCMI
- *calc_mi_with_P:* calculate MI given a camera matrix P
- *calc_mi_with_theta_6params:* calculate MI given 6 exterior parameters
- *calc_ncmi_with_P:* calculate NCMI given a camera matrix P
- *SR_for_MI:* super-resolution to generate i-image (later used for the maximization of MI), calling SR_by_grad.m
- *SR_for_NCMI:* super-resolution to generate i-image and z-image (later used for the maximization of NCMI), calling SR_by_grad.m
- *SR_by_grad:* run the FISTA algorithm for the super-resolution
- *mi:* measure Mutual Information between two images
- *ncmi:* measure Normalized Combined Mutual Information between three images
- *patch_based_interpolation:* IDW-based interpolation for the smoothing between patches

## Super-resolution (SR)
- *depth_map_SR* and *run_SR_by_grad:* scripts for the super-resolution of a depth-map
- *cost_func_Dxy:* cost function defined by the sum squared directional gradients (SSGD)
- *gradient_Dxy:* gradient of the cost function SSGD
- *gradient_Dxy_fast:* gradient of the cost function SSGD (faster version)
- *soft_thres:* shrinkage operator

## Snake Model
*(The implementation of original snake models was developed by D. Kroon University of Twente. To download the original codes, please visit [link](https://www.mathworks.com/matlabcentral/fileexchange/28149-snake-active-contour?s_tid=prof_contriblnk))*
- *building_seg_isprs:* preliminary extraction of building points on ISPRS Vaihingen dataset
- *BE_by_snakes_ISPRS:* call snake model for building extraction
- *Snake2D:* run snake model algorithm
- *ExternalForceImage2D:* compute external image-based energy function
- *GetContourNormals2D:* compute the normals of the contour points using the neighbouring points (used in traditional balloon force)
- *GVFOptimizeImageForces2D:* gradient vector flow (GVF)
- *ImageDerivatives2D:* gaussian based image derivatives
- *InterpolateContourPoints2D:* resamples a given contour to a smooth contour of uniform sampled points
- *MakeContourClockwise2D:* make contour clockwise
- *SnakeInternalForceMatrix2D:* compute internal energy function
- *SnakeMoveIteration2D:* compute the snake movement at iteration
- *snake_polygonization:* apply a polygonization (Douglas-Peucker or Dutter) on snake model results
- *polygonization_Dutter:* polygonization algorithm proposed by (Dutter 2007)

## SRSM
*(Many of the codes are the same with the ones in folder **Snakes** and **SR**)*
- *SRSM_prototype:* run SRSM to extract buildings 
- *SRSM_prototype_ISPRS:* run SRSM to extract buildings on ISPRS Vaihingen dataset 
- *SR_z*: call super-resolution algorithm to generate z-image
- *SR_by_grad:* super-resolution algorithm
- *SnakeMoveIteration2D_pixel_wise:* compute the snake movement at iteration with option for modified balloon force
- *las2txt:* convert LAS files to TXT

## Misc.
- *calc_IoU:* calculate Intersection-over-Union metric
- *calcPolisMetric*: calculate PoLIS metric
- *convertGeoCoorMatCoor:* function to convert georeferenced coordinates into matrix coordinates
- *convertMatCoorGeoCoor:* function to convert matrix coordinates into georeferenced coordinates
- *countObject*: function to count an object (for object-based accuracy) if it overlaps 50% of its reference boundary. 
- *display_TP_FP_FN:* display true positives, false positives, and false negatives, for building extraction results evaluated with a ground truth
- *dist_euclid:* compute Euclidean distance from a point to a line
- *dpsimplify:* Douglas-Peucker line simplification algorithm
- *findNearestPoints:* find the nearest points within a cloud to a given location (x,y)
- *findPoints:* find the points within a cloud that have horizontal coordinates within a rectangle [xmin xmax ymin ymax]
- *HausdorffDist:* measure Hausdorff distance between two sets of points
- *load_datasets:* load Quebec datasets (image and LiDAR)
- *load_datasets_isprs:* load Vaihingen datasets (image and LiDAR)
- *maxk_myfunc:* find k largest elements of an array (if maxk is not avaible)
- *mink_myfunc:* find k smallest elements of an array (if mink is not avaible)
- *project_pts_given_TM:* project points given a transformation model (i.e. camera matrix P)

## R2Sonic Challenge
- *readXTFFiles_R2Sonic_Truepix:* process XTF files to read and extract TruePix data, calling R2SONIC_TRUEPIX
- *R2SONIC_TRUEPIX:* function to process XTF files and read TruePix data header
- *R2SONIC_TRUEPIX_D:* function to process XTF files and read TruePix data packets
- *XTFPINGHEADER:* function to process XTF files and read the sonar header (metadata, config, etc.)
- *XTFATTITUDEDATA:* function to process XTF files and read attitude data
- *XTFPOSRAWNAVIGATION:* function to process XTF files and read raw navigation data

- *createTruePix_sediments:* create bathymetric point cloud from TruePix data
- *createWCI_from_Truepix:* create water column fan-shaped images from TruePix data
- *truepix_sediments_Sv:* attempt to compute backscattering strength per unit volume of ensonified water, inspired by [Simmons et al 2010]
- *localisation_gage:* try to locate the gage and extract the TruePix data in proximity
- *physical_samples_analysis:* read and plot physical samples, in attempt to analyze them and construct historical data