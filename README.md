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

## Quote
If using this script, please cite: 
Bessin, Z.; Dedieu, J.-P.; Arnaud, Y.; Wagnon, P.; Brun, F.; Esteves, M.; Perry, B.; Matthews, T. Processing of VENµS Images of High Mountains: A Case Study for Cryospheric and Hydro-Climatic Applications in the Everest Region (Nepal). Remote Sens. 2022, 14, 1098. https://doi.org/10.3390/rs14051098
