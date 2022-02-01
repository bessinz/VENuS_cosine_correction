#!/usr/bin/env python
# -*- coding: utf-8 -*-

#@author : Zoe BESSIN
#LINUX VERSION

import os
import numpy as np
import otbApplication



#-------------------------------------------------------------------------------------------------------------
#---------------------------------------TO GET ALL KIND OF GRANULE PATHS--------------------------------------
#-------------------------------------------------------------------------------------------------------------

def sensor_data_path_provider(sensor,path_images):    
    path_granule = []
    
    ## At this step, path to the sensor directory is provided, it will now go deeper until path of each files are gathered    
    if sensor =='Venus':
        images_all=os.listdir(path_images)
                
        for ii in range(len(images_all)):
            if images_all[ii] != "images_zip":
                path_inter=os.listdir(path_images+images_all[ii]) 
                idx2 = np.asarray([cpt2 for cpt2, s2 in enumerate(path_inter) if "SRE" in s2 and "B11" in s2 and ".aux.xml" not in s2])
                path_granule.append(path_images+str(images_all[ii])+'/'+str(path_inter[idx2[-1]]))
#                idx3 = np.asarray([cpt2 for cpt2, s2 in enumerate(path_inter) if "FRE" in s2 and "B11" in s2 and ".aux.xml" not in s2])
#                path_granule_FRE.append(path_images+str(images_all[ii])+'/'+str(path_inter[idx3[-1]]))
#                idx4 = np.asarray([cpt2 for cpt2, s2 in enumerate(path_inter) if "FRE" in s2 and "B7" in s2 and ".aux.xml" not in s2])
#                path_granule_B7.append(path_images+str(images_all[ii])+'/'+str(path_inter[idx4[-1]]))
#                path_inter2=os.listdir(path_images+str(images_all[ii])+'/DATA/')
#                path_pv = os.listdir(path_images+str(images_all[ii])+'/DATA/'+str(path_inter2[1]))
#                idx_pv = np.asarray([cpt2 for cpt2, s2 in enumerate(path_pv) if "CLD" in s2 and '.xml' not in s2])
#                path_granule_cloud.append(path_images+str(images_all[ii])+'/DATA/'+str(path_inter2[1])+'/'+str(path_pv[idx_pv[-1]]))
#                idx_pv = []
        print(sensor+' archive successfully found')
    return path_granule#, path_granule_FRE, path_granule_B7, path_granule_cloud





def haralick(path_granule):
    img_in = path_granule
    img_out = path_res+os.path.basename(path_granule)[:-4]+"_Haralick.tif"
    
    app = otbApplication.Registry.CreateApplication("HaralickTextureExtraction")
    
    app.SetParameterString("in", img_in)
    app.SetParameterInt("channel", 1)
    app.SetParameterInt("step", 1)
    app.SetParameterInt("parameters.xrad", 2)
    app.SetParameterInt("parameters.yrad", 2)
    app.SetParameterInt("parameters.xoff", 1)
    app.SetParameterInt("parameters.yoff", 1)
    app.SetParameterInt("parameters.nbbin", 8)
    app.SetParameterString("texture","simple")
    app.SetParameterString("out", img_out)
    app.ExecuteAndWriteOutput()
    return str(img_out)
    
def keep_energy_band(image_in, image_out, band_to_keep):
    
    keep_cast_shadow = os.path.join(path_gdal+'gdal_translate -b '+str(band_to_keep)+' '+image_in+' '+image_out)
    print(keep_cast_shadow)
    os.system(keep_cast_shadow)
    

## Provides path for the rest of the code
sensor = 'Venus'
path_general = '/home/bessinz/'
path_images  = path_general+'IMG2/' 
path_gdal = '/home/bessinz/anaconda2/envs/orfeo/bin/'


path_granule = sensor_data_path_provider(sensor,path_images)
nb_images = len(path_granule)
print(nb_images)
for i in range(nb_images):
    pg = path_granule[i]
    print(pg)
    
    
    #path_res = path_general+'Results_haralick/'
    path_res = path_images+os.path.basename(os.path.dirname(pg))+'/'
    #Make the directory to store the results 
    if not os.path.exists(path_res):
        os.makedirs(path_res)
    
    #img_in = pg
    #img_out = path_res+os.path.basename(pg)[:-4]+"_Haralick.tif"
    
    img_out_har = haralick(pg)
    img_out = path_res+os.path.basename(img_out_har)[:-4]+"_Energy.tif"
    keep_energy_band(img_out_har, img_out, 1)
    os.remove(img_out_har)
