# VENuS_cosine_correction

## Description
### VENuS_TopCos_LINUX_20200813.py
 This code is used to process a cosine correction on VENµS (Vegetation and Environment New Micro Satellite) SRE (Surface REflectance) product using the 8-m HMA (High Mountain Asia) DEM (Digital Elevation Model); improve the CNES CLM_XS cloud mask (cloud shadows are not included in this mask) and 
 compute SCA (Snow Cover Area) maps in the Khumbu valley (Nepal).
### Haralick_Energy_LINUX_20200813.py
 This code is used to retrieve only the Energy matrix using HaralickTextureExtraction orfeo toolbox tool (https://www.orfeo-toolbox.org/CookBook/Applications/app_HaralickTextureExtraction.html)

## Installation
 These codes has been developed on Python 3 using Linux. A Windows version will follow.
 
## Data
VENµS tiles can be dowloaded on https://theia.cnes.fr/atdistrib/rocket/#/search?collection=VENUS. To automatically download Theia products, please use the code developed by @olivierhagolle here: https://github.com/olivierhagolle/theia_download

## Citation
If using this script, please cite: Z. Bessin, J.P. Dedieu, Y. Arnaud, P. Wagnon, F. Brun, M. Esteves, L. B. Perry, T. Matthews, "Processing of VENµS Images in High Mountain: A Case Study for Cryospheric and Hydro-Climatic applications in the Everest Region (Nepal)", Remote Sensing, under review.
