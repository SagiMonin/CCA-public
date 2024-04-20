## Continuous Cost Aggregation (CCA) for Dual-Pixel Disparity Extraction ##
Link to paper: https://arxiv.org/abs/2306.07921


# Description folders
1. CSGM : folder containing all implementations divided based on data used (DP-DSLR/PHONE, Stereo).
2. demo_DP: folder with different demo files for DP data: DSLR + Google phone
3. demo_Stereo: folder with demo for traditional stereo - comparison between C-SGM and SGM. Filter is the same for both SGM and C-SGM (names are changed).
4. DP_data_example: example data of DP-DSLR/PHONE
5. Stereo_data_example: example data of Stereo middlebury at Quarter resultion. 
6. evaluate_funcs: functions used for evaluation
7. openCvMatlab : this is complied for windows, might need to mex-install in different eviorment 

# Demos: 
1. demo_csgm_DSLR_quantitive - test on DSLR data set (1 in data)
2. demo_csgm_phone+filter_test_phone.py: test on phone data (2 in data).
3. demo_csgm_stereo+bilaterl_filter_csgm.py: test on middlebury data (3 in data).
4. demo_sgm_stereo+bilaterl_filter_sgm.py: Adapted SGM on middlebury data (3 in data).

# How to use
1. Add all files to path.
2. Run one of the above demos.
3. OpenCV for Matlab might need to be installed to load images.

# Data:

1. DP-data:
There is an example for one image from ICCP-2020 paper (DSLR): DP_data_example\ICCP2020\GT_data.
Rest of data can be found in: 
https://github.com/abhijithpunnappurath/dual-pixel-defocus-disparity/blob/master/README.md

2. Google-phone data: 
There is one example from google data from ICCV 2019 paper: DP_data_example\google2019\test
Rest of data can be found in:
https://github.com/google-research/google-research/blob/master/dual_pixels/README.md

3. Stereo middlebury:
https://vision.middlebury.edu/stereo/submit3/

# Different algorithms (compared in paper):

1. DLP: Learning Single Camera Depth Estimation using Dual-Pixels.
Implementation not available, predictions on google data set found in: https://github.com/google-research/google-research/blob/master/dual_pixels/README.md

2. SDoF: Synthetic Depth-of-Field with a Single-Camera Mobile Phone
Implementation not available. 

3. DPE: Dual Pixel Exploration: Simultaneous Depth Estimation and Image Restoration
Implementation in: https://github.com/panpanfei/Dual-Pixel-Exploration-Simultaneous-Depth-Estimation-and-Image-Restoration

4. DPdisp: Modeling Defocus-Disparity in Dual-Pixel Sensors
Implementation found in: https://github.com/abhijithpunnappurath/dual-pixel-defocus-disparity/blob/master/README.md

5. SGM: 
SGM implementation, customed to use BT score instead of SAD. 
Implementation adapted from: https://github.com/kobybibas/semi_global_matching
