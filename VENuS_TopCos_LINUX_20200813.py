#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on Tue Mar 10 14:22:13 2020

#@author: Zoe BESSIN
#LINUX VERSION


#=======================================================================================================================
#===================================================LIBRARIES===========================================================
#=======================================================================================================================

import os
import ogr
import gdal, osr
import osgeo.gdal
osgeo.gdal.GetDriverByName
import numpy as np
import xml.etree.ElementTree as ET
import numpy


#---------------summury of the fonctions defined in this program---------------
#sensor_data_path_provider(sensor,path_images)
#findGranuleMetadata(sensor,path_granule_SRE_B11)
#getprojection(raster)
#slope_aspect(path_input_DEM)
#createBuffer(inputfn, outputBufferfn, bufferDist)
#create_shp_tile(path_granule, srcnodata)
#normalisation(raster_in, raster_out, dstmin, dstmax)
#clip_to_tile_and_normalized(path_granule, path_granule_SRE_B11)
#clip_to_tile(path_granule, path_granule_SRE_B11, srcnodata, dstnodata)
#create_matrix_of_ones(path_granule_clip)
#extract_cloud_mask(path_granule_cloud, path_granule_SRE_B11)
#calcul_NDVI(path_granule_B11, path_granule_B7)
#Haralick_cloud_masking_Venus(Haralick, CLM_MASK, NDVI, FRE_B11_clip, path_granule_SRE_B11, path_granule_SRE_B11_clip)
#lakes_mask(rasterin, path_granule_SRE_B11)
#TopCor_Cosine(path_granule, path_granule_FRE, Zs,As, path_input_DEM, path_granule_SRE_B11)
#apply_cld_mask(rasterin, CLM_MASK_VENUS, path_granule_SRE_B11)
#SCA_maker(NDVI, Energy, CLM_MASK_VENUS, path_granule_SRE_B11)
#clip_SCA_to_watershed(SCA, path_granule_SRE_B11, srcnodata, dstnodata)

#=======================================================================================================================
#===================================================FONCTIONS===========================================================
#=======================================================================================================================

#-------------------------------------------------------------------------------------------------------------
#---------------------------------------TO GET ALL KIND OF GRANULE PATHS--------------------------------------
#-------------------------------------------------------------------------------------------------------------

def sensor_data_path_provider(sensor,path_images):    
    path_granule_SRE_B11 = []
    path_granule_FRE_B11 = []
    path_granule_SRE_B7 = []
    path_granule_FRE_B7 = []
    path_granule_CLM = []
    path_granule_Energy = []
    path_granule_CLDCOV = []
    ordre_images = {}
    ## At this step, path to the sensor directory is provided, it will now go deeper until path of each files are gathered    
    if sensor =='Venus':
        images_all=os.listdir(path_images)
                
        for ii in range(len(images_all)):
            if images_all[ii] != "zip" and images_all[ii] != "search.json" and images_all[ii] != "Results_Energy" :
                path_inter=os.listdir(path_images+images_all[ii])
                idx2 = np.asarray([cpt2 for cpt2, s2 in enumerate(path_inter) if "SRE" in s2 and "B11" in s2 and "CLOUD_MASK_VENUS_8bits_m" not in s2 and ".aux.xml" not in s2 and "_Haralick_Energy" not in s2])
                path_granule_SRE_B11.append(path_images+str(images_all[ii])+'/'+str(path_inter[idx2[-1]]))
                
                
                
                idx3 = np.asarray([cpt2 for cpt2, s2 in enumerate(path_inter) if "FRE" in s2 and "B11" in s2 and ".aux.xml" not in s2])
                path_granule_FRE_B11.append(path_images+str(images_all[ii])+'/'+str(path_inter[idx3[-1]]))
                
                idx5 = np.asarray([cpt2 for cpt2, s2 in enumerate(path_inter) if "SRE" in s2 and "B7" in s2 and ".aux.xml" not in s2])
                path_granule_SRE_B7.append(path_images+str(images_all[ii])+'/'+str(path_inter[idx5[-1]]))
                
                idx4 = np.asarray([cpt2 for cpt2, s2 in enumerate(path_inter) if "FRE" in s2 and "B7" in s2 and ".aux.xml" not in s2])
                path_granule_FRE_B7.append(path_images+str(images_all[ii])+'/'+str(path_inter[idx4[-1]]))
                
                idx8 = np.asarray([cpt2 for cpt2, s2 in enumerate(path_inter) if "CLOUD_MASK_VENUS_8bits_m.tif" in s2 or "CLOUD_MASK_VENUS_8bits_m.TIF" in s2 and ".pox" not in s2 and ".aux.xml" not in s2])
		print(idx8)
		path_granule_CLDCOV.append(path_images+str(images_all[ii])+"/"+str(path_inter[idx8[-1]]))  
               
                path_inter2=os.listdir(path_images+str(images_all[ii])+'/MASKS/')
                
                
                idx_pv2 = np.asarray([cpt2 for cpt2, s2 in enumerate(path_inter2) if "CLM_XS" in s2 and '.xml' not in s2])
                path_granule_CLM.append(path_images+str(images_all[ii])+'/MASKS/'+str(path_inter2[idx_pv2[-1]]))
                
                idx6 = np.asarray([cpt2 for cpt2, s2 in enumerate(path_inter) if "SRE" in s2 and "B11_Haralick_Energy" in s2 and ".aux.xml" not in s2])
                path_granule_Energy.append(path_images+str(images_all[ii])+'/'+str(path_inter[idx6[-1]]))
        
                #path_inter3=os.listdir(path_images+str(images_all[ii])+'/Results/')
                
                
                #idx_pv3 = np.asarray([cpt2 for cpt2, s2 in enumerate(path_inter3) if "CLOUD_MASK_VENUS_8bits_m.tif" in s2 or "CLOUD_MASK_VENUS_8bits_m.TIF" in s2 and '.pox' not in s2 and '.aux.xml' not in s2])
                #path_granule_CLDCOV.append(path_images+str(images_all[ii])+'/Results/'+str(path_inter3[idx_pv3[-1]]))
              
        f= open("ordre_images.txt","a+")
        for e in path_granule_SRE_B11:
            ordre_images[path_granule_SRE_B11.index(e)]=e
            f.write(str(path_granule_SRE_B11.index(e))+' '+str(e)+'\n')
        f.close()
        #print(ordre_images)
	print(path_granule_CLDCOV)
        print('\n----'+sensor+' archive successfully found. To see image order, have a look to ordre_images.txt')
    return path_granule_SRE_B11, path_granule_FRE_B11, path_granule_SRE_B7, path_granule_FRE_B7, path_granule_CLM, path_granule_Energy, path_granule_CLDCOV


#-------------------------------------------------------------------------------------------------------------
#--------------------------------------TO GET ZENITH AND AZIMUT ANGLE-----------------------------------------
#-------------------------------------------------------------------------------------------------------------

def findGranuleMetadata(sensor,path_granule_SRE_B11):
    # This routine goes into the index of all existing granules for a given sensor and return the solar zenith and azimuth angle
    ##
    ##
    ## sensor = name of sensor according to previous used typo. (Sentinel2, Landsat-WRS_1, Landsat-WRS_2, ASTER etc...)
    ## path_granule_SRE_B11 = list of granule path from root.
    
    print("\n2-Searching for metadata of the image in catalog...")
        
    if sensor == 'Venus':
        idx1=np.asarray([cpt for cpt, s2 in enumerate(os.listdir(os.path.dirname(path_granule_SRE_B11))) if '.xml' in s2 and '.aux.xml' not in s2])
        if (len(idx1))>0:
            collection_file = (os.path.dirname(path_granule_SRE_B11))+'/'+(os.listdir(os.path.dirname(path_granule_SRE_B11)))[idx1[-1]]
            tree = ET.parse(collection_file)
            root = tree.getroot()
            for i,test in enumerate(root.iter('ZENITH_ANGLE')):
                if i == 0:
                    Zs = float(test.text)
                    
            for i,test in enumerate(root.iter('AZIMUTH_ANGLE')):
                if i == 0:
                    As = test.text
                    
        elif (len(idx1)) == 0:
            print('Could not grab metadata of this granule')
            Zs = np.nan
            As = np.nan
    print('----As = '+str(As)+'°')
    print('----Zs = '+str(Zs)+'°')
    return float(Zs), float(As)



#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------TO GET LAYER PROJECTION-------------------------------------------
#-------------------------------------------------------------------------------------------------------------

def getprojection(raster):
    try:
        print('----Get the projection of '+os.path.basename(raster))
        ds = gdal.Open(raster) #open the granule
        prj = ds.GetProjection() #get its projection
        geotransform = ds.GetGeoTransform() #
        ds = None
        pixelwidth = geotransform[1] #give us the spatial resolution
        pixelheight = geotransform[-1]
        srs=osr.SpatialReference(wkt=prj)
        #print('\nThe layer you load has :')
        #print('-Geographic Coordinate System (GCS) = EPSG:'+str(srs.GetAttrValue("GEOGCS|AUTHORITY", 1))) #give us the GCS (Geographic Coordinate System of the data)
        #print('-Projected Coordinate Systems (PROJCS) EPSG:'+str(srs.GetAttrValue("PROJCS|AUTHORITY", 1))+'('+str(srs.GetAttrValue("PROJCS|PROJECTION", 1))+')') #give us the Projected Coordinate System of the data
        SRS='EPSG:'+srs.GetAttrValue("PROJCS|AUTHORITY", 1)
        #print('SRS '+os.path.basename(raster)+'\n'+SRS)
        return pixelwidth, pixelheight, SRS
    except:
        print('Could not grab SRS for'+raster)



def slope_aspect(path_input_DEM):
    #Calculate slope and aspect then cropped to tiles so you don't loose accuracy at the border
    slope_tmp = TopCor_material_Fn+'/'+os.path.basename(path_input_DEM[:-4])+'_slope_sensor.TIF'    
    aspect_tmp = TopCor_material_Fn+'/'+os.path.basename(path_input_DEM[:-4])+'_aspect_sensor.TIF'
    
    if not os.path.isfile(slope_tmp):
        print("----Creation of slope")      
        slope = os.path.join(path_gdal+'gdaldem slope -compute_edges '+path_input_DEM+' '+slope_tmp) #work for Windows
        os.system(slope)
        
    if not os.path.isfile(aspect_tmp):            
        print("----Creation of aspect")    
        aspect = os.path.join(path_gdal+'gdaldem aspect -zero_for_flat -compute_edges '+path_input_DEM+' '+aspect_tmp) #work for Windows
        os.system(aspect)
    return str(slope_tmp), str(aspect_tmp)    

 
        
#-------------------------------------------------------------------------------------------------------------
#---------------------------------------------TO CREATE A BUFFER----------------------------------------------
#-------------------------------------------------------------------------------------------------------------
        
def createBuffer(inputfn, outputBufferfn, bufferDist):
    inputds = ogr.Open(inputfn) #open the file
    inputlyr = inputds.GetLayer()#get its layers

    shpdriver = ogr.GetDriverByName('ESRI Shapefile') #define the format
    if os.path.exists(outputBufferfn):
        shpdriver.DeleteDataSource(outputBufferfn) #remove the existing file to replace it
    outputBufferds = shpdriver.CreateDataSource(outputBufferfn)
    bufferlyr = outputBufferds.CreateLayer(outputBufferfn, geom_type=ogr.wkbPolygon) #create the polygon geometry
    featureDefn = bufferlyr.GetLayerDefn()

    for feature in inputlyr:
        ingeom = feature.GetGeometryRef()
        geomBuffer = ingeom.Buffer(bufferDist, 8)
        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(geomBuffer)
        bufferlyr.CreateFeature(outFeature)
        outFeature = None
      
        
#-------------------------------------------------------------------------------------------------------------
#---------------------------------------TO CREATE SHAPEFILE USED TO CROP--------------------------------------
#-------------------------------------------------------------------------------------------------------------
    
def create_shp_tile(path_granule, srcnodata):
    print('\n1-Create the shapefile used to crop')
    # set name of outfile shapefile
    granule_extent_shapefile_tmp = path_res+os.path.basename(path_granule[:-4])+'_tmp.'
    granule_extent_shapefile = path_res+os.path.basename(path_granule[:-4])+'.'
    granule_extent_shapefile_25m_tmp = path_res+os.path.basename(path_granule[:-4])+'_25m_tmp.'
    granule_extent_shapefile_25m = path_res+os.path.basename(path_granule[:-4])+'_25m.'
    if not os.path.exists(granule_extent_shapefile_25m+'shp'):        
        # set the name of temporary layers
        granule_only_valid_fn = path_res+os.path.basename(path_granule[:-4])+'_only_valid.TIF' # Set tmp name
    
        # First remove nodata from the granule
        # Use gdalwarp fonction to delete nodata in source data and introduce an alpha band to replace nodata. 
        # You get granule_only_valid_fn.TIF in output
        print("----Remove nodata from SRE")
        remove_nodata = os.path.join(path_gdal+'gdalwarp -overwrite -srcnodata '+srcnodata+' -dstalpha '+path_granule+' '+granule_only_valid_fn) #work for Windows
        os.system(remove_nodata)
      
        # then grab the real extent of valid pixels by using gdal_polygonize.py witch give us a shapefile delimiting the tile.
        # Use the raster and the band number choosen and the output is in granule_extent_shapefile.shp
        # Band 2 is the band alpha witch contains only 2 different features
        print("----Create outline shapefile")
        grab_extent_valid_pix = os.path.join('python '+path_gdal+'gdal_polygonize.py '+granule_only_valid_fn+' -b 2 -f "ESRI Shapefile" '+granule_extent_shapefile_tmp+'shp') #work for Windows
        os.system(grab_extent_valid_pix)
    
        print("----Select the right feature to crop the rasters")
        command2 = os.path.join(path_gdal+'ogr2ogr -overwrite -where DN!="0" '+granule_extent_shapefile+'shp '+granule_extent_shapefile_tmp+'shp')
        os.system(command2)
    
        print("----Create a buffer of 25 m from the shapefile, in order to suppress aberrant values introduced by MAJA on the bounderies of the tile")
        createBuffer(granule_extent_shapefile+'shp', granule_extent_shapefile_25m_tmp+'shp', -25.0)
        
        attribute_proj = os.path.join(path_gdal+'ogr2ogr -f "ESRI Shapefile" -t_srs EPSG:32645 -s_srs EPSG:32645 '+granule_extent_shapefile_25m+'shp '+granule_extent_shapefile_25m_tmp+'shp ')
        os.system(attribute_proj)
        
        print("----Remove the temporary shapefile")
        os.remove(granule_extent_shapefile_tmp+'shp')
        os.remove(granule_extent_shapefile_tmp+'shx')
        os.remove(granule_extent_shapefile_tmp+'dbf')
        os.remove(granule_extent_shapefile_tmp+'prj') 
        os.remove(granule_extent_shapefile_25m_tmp+'shp')
        os.remove(granule_extent_shapefile_25m_tmp+'shx')
        os.remove(granule_extent_shapefile_25m_tmp+'dbf')
        os.remove(granule_extent_shapefile+'shp')
        os.remove(granule_extent_shapefile+'shx')
        os.remove(granule_extent_shapefile+'dbf')
        os.remove(granule_extent_shapefile+'prj')     
        os.remove(granule_only_valid_fn)
    return str(granule_extent_shapefile_25m)
    
        

#-------------------------------------------------------------------------------------------------------------
#-----------------------------------------TO NORMALIZE RASTER VALUES------------------------------------------
#-------------------------------------------------------------------------------------------------------------

def normalisation(raster_in, raster_out, dstmin, dstmax):
    geotif = gdal.Open(raster_in)
    srcband = geotif.GetRasterBand(1)
    stats = srcband.GetStatistics(True, True)
    print("[ STATS %s] =  Minimum=%.3f, Maximum=%.3f, Mean=%.3f, StdDev=%.3f" % (os.path.basename(raster_in), stats[0], stats[1], stats[2], stats[3]))
    normalize = os.path.join(path_gdal+'gdal_translate -a_srs EPSG:32645 -scale '+str(stats[0])+' '+str(stats[1])+' '+str(dstmin)+' '+str(dstmax)+' -a_nodata 20000 '+raster_in+' '+raster_out)
    os.system(normalize)


#-------------------------------------------------------------------------------------------------------------
#----------------------------------TO CLIP AND NORMALIZED A RASTER----------------------------------
#-------------------------------------------------------------------------------------------------------------
    
    
def clip_to_tile_and_normalized(path_granule, path_granule_SRE_B11):
    path_granule_clip = path_res+os.path.basename(path_granule[:-4])+'_clip.TIF'
    path_granule_clip_norm = path_res+os.path.basename(path_granule[:-4])+'_clip_norm.TIF'
    
    if not os.path.exists(path_granule_clip):
    
        # First, open the granule and grab projection, Xres ...
        S_pixelwidth, S_pixelheight, SRSsensor = getprojection(path_granule_SRE_B11)
        
        print("----Clip of "+os.path.basename(path_granule))
        clip_FRE_to_tile = os.path.join(path_gdal+'gdalwarp -overwrite --config GDALWARP_IGNORE_BAD_CUTLINE YES -crop_to_cutline -cutline '
                     +granule_extent_shapefile_25m+'shp '+'-r '+Rmethod+' -tr '+str(S_pixelwidth)+' '+str(S_pixelheight)+' -s_srs '+SRSinput_DEM+' -t_srs '
                                        +SRSsensor+' -srcnodata -10000 -dstnodata 20000 '+path_granule+' '+path_granule_clip)
        os.system(clip_FRE_to_tile)
    
#    if not os.path.exists(path_granule_clip_norm):
#
#        print("Raster normalisation of "+os.path.basename(path_granule))
#        normalisation(path_granule_clip, path_granule_clip_norm, 0.0, 10000.0)
        
    return str(path_granule_clip)#, str(path_granule_clip_norm)

def clip_to_tile(path_granule, path_granule_SRE_B11, srcnodata, dstnodata):
    path_granule_clip = path_res+os.path.basename(path_granule[:-4])+'_clip.TIF'
    
    if not os.path.exists(path_granule_clip):
    
        # First, open the granule and grab projection, Xres ...
        S_pixelwidth, S_pixelheight, SRSsensor = getprojection(path_granule_SRE_B11)
        
        print("Clip of "+os.path.basename(path_granule))
        clip_FRE_to_tile = os.path.join(path_gdal+'gdalwarp -overwrite --config GDALWARP_IGNORE_BAD_CUTLINE YES -crop_to_cutline -cutline '
                     +granule_extent_shapefile_25m+'shp '+'-r '+Rmethod+' -tr '+str(S_pixelwidth)+' '+str(S_pixelheight)+' -s_srs '+SRSinput_DEM+' -t_srs '
                                        +SRSsensor+' -srcnodata '+str(srcnodata)+' -dstnodata '+str(dstnodata)+' '+path_granule+' '+path_granule_clip)
        print(clip_FRE_to_tile)
        os.system(clip_FRE_to_tile)
        
    return str(path_granule_clip)
    
#-------------------------------------------------------------------------------------------------------------
#------------------------------------------TO CREATE A MATRIX OF ONES-----------------------------------------
#-------------------------------------------------------------------------------------------------------------

def create_matrix_of_ones(path_granule_clip): #path of the file you want to keep the size for your matrix
    # Create an array of ones
    print("----Creation of an array of ones")
    # Define the filename of the output array
    path_granule_sensor_array_copy2 = path_res+os.path.basename(path_granule_clip[:-4])+'_array2.tif'
    if not os.path.exists(path_granule_sensor_array_copy2):
        ds5 = gdal.Open(path_granule_clip) 
        path_granule_sensor_array2 = ds5.GetRasterBand(1).ReadAsArray().astype(numpy.int16)
        # Fill the array with ones
        path_granule_sensor_array2 = numpy.ones((numpy.size(path_granule_sensor_array2,0),numpy.size(path_granule_sensor_array2,1)), dtype = numpy.float32)
        #set the format of output data
        driver5 = gdal.GetDriverByName('GTiff')
        #create a copy of the raster
        ds5_out = driver5.CreateCopy(path_granule_sensor_array_copy2, ds5) 
        # write the modified array to the raster
        ds5_out.GetRasterBand(1).WriteArray(path_granule_sensor_array2)
        # set the NoData metadata flag
        ds5_out.GetRasterBand(1).SetNoDataValue(20000)
        # clear the buffer, and ensure file is written
        ds5_out.FlushCache()
    return str(path_granule_sensor_array_copy2)

#-------------------------------------------------------------------------------------------------------------
#------------------------------------------TO EXTRACT THE CLOUD BIT-------------------------------------------
#-------------------------------------------------------------------------------------------------------------



def extract_cloud_mask(path_granule_cloud, path_granule_SRE_B11):
    print('Extract the bit 7 from CLM')
    #path_cloud_mask = path_res+os.path.basename(path_granule_SRE_B11[:-4])+'_'+os.path.basename(path_granule_cloud[:-4])+'_mask.TIF'    
    path_CLM_mask = path_res+os.path.basename(path_granule_SRE_B11[:-4])+'_CLM_mask.TIF' 
    
#    path_cloud_mask_clip = path_cloud_mask_clip = path_cloud_mask[:-4]+'_clip.TIF'
#    if not os.path.exists(path_cloud_mask_clip):
#        path_cloud_mask_clip = clip_to_tile(path_granule_cloud, path_granule_SRE_B11, 0, 0)
        
    if not os.path.exists(path_CLM_mask):
        extract_bit = os.path.join('python '+path_gdal_calc_beta+'gdal_calc_modified.py --overwrite --type=Int16 -A '+path_granule_cloud+' --outfile='
                                 +path_CLM_mask+' --calc="numpy.bitwise_and(A, (10000000))" --CalcAsDT --NoDataValue=0')
        
        os.system(extract_bit)
              
    return str(path_CLM_mask)

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------TO CALCULATE NDVI-------------------------------------------
#-------------------------------------------------------------------------------------------------------------
            
def calcul_NDVI(path_granule_B11, path_granule_B7):
    print('Create the NDVI raster')
    #filename of the output : B7 cropped
    path_granule_B7_clip = path_res+os.path.basename(path_granule_B7[:-4])+'_clip.TIF'
    if not os.path.exists(path_granule_B7_clip):
        #crop the B7 band
        clip_B7_to_tile = os.path.join(path_gdal+'gdalwarp -overwrite --config GDALWARP_IGNORE_BAD_CUTLINE YES -crop_to_cutline -cutline '
                     +granule_extent_shapefile_25m+'shp '+'-r '+Rmethod+' -tr '+str(SRE_B11_pixelwidth)+' '+str(SRE_B11_pixelheight)+' -s_srs '+SRSinput_DEM+' -t_srs '+SRS_SRE_B11+
                                              ' -srcnodata -10000 -dstnodata 20000 '+path_granule_B7+' '+path_granule_B7_clip)
        os.system(clip_B7_to_tile)
    
    #filename of the output : B11 cropped
    path_granule_B11_clip = path_res+os.path.basename(path_granule_B11[:-4])+'_clip.TIF'
    if not os.path.exists(path_granule_B11_clip):
        #crop the B11 band
        clip_B11_to_tile = os.path.join(path_gdal+'gdalwarp -overwrite --config GDALWARP_IGNORE_BAD_CUTLINE YES -crop_to_cutline -cutline '
                     +granule_extent_shapefile_25m+'shp '+'-r '+Rmethod+' -tr '+str(SRE_B11_pixelwidth)+' '+str(SRE_B11_pixelheight)+' -s_srs '+SRSinput_DEM+' -t_srs '+SRS_SRE_B11+
                                              ' -srcnodata -10000 -dstnodata 20000 '+path_granule_B11+' '+path_granule_B11_clip)
        os.system(clip_B11_to_tile)

    #filename of the output : NDVI
    NDVI = path_res+os.path.basename(path_granule_B11)[:-4]+'_NDVI.TIF'
    if not os.path.exists(NDVI):
        #calculate the NDVI
        ndvi = os.path.join('python '+path_gdal_calc_beta+'gdal_calc_modified.py --overwrite --type=float32 -A '+path_granule_B11_clip+' -B '+path_granule_B7_clip+' --outfile='
                                 +NDVI+' --calc="(numpy.divide((A-B),(A+B), out=numpy.zeros(A.shape), where= ((A+B)<>0)))" --CalcAsDT --NoDataValue=0')
        os.system(ndvi)
    os.remove(path_granule_B7_clip)
    os.remove(path_granule_B11_clip)
    return str(NDVI)#, str(path_granule_B11_clip), str(path_granule_B7_clip)
    

#-------------------------------------------------------------------------------------------------------------
#---------------------------------------------TO CREATE THE CLOUD MASK----------------------------------------
#-------------------------------------------------------------------------------------------------------------


def Haralick_cloud_masking_Venus(Haralick, CLM_MASK, NDVI, FRE_B11_clip, path_granule_SRE_B11, path_granule_SRE_B11_clip):
    print('Create the cloud mask')        
    path_granule_sensor_array_copy2 = path_res+os.path.basename(path_granule_SRE_B11[:-4])+'_array2.tif'
        
    #check if matrix of ones exist or not
    if not os.path.isfile(path_granule_sensor_array_copy2):        
        path_granule_sensor_array_copy2 = create_matrix_of_ones(path_granule_SRE_B11)
    
    #path_granule_sensor_array_copy2_clip = clip_to_tile(path_granule_sensor_array_copy2, path_granule_SRE_B11)

    #filename of the output
    path_cloud_mask_clip = CLM_MASK[:-4]+'_clip.TIF'
    print(path_cloud_mask_clip)
    #clip path_cloud_mask to the tile
    path_cloud_mask_clip = clip_to_tile(CLM_MASK, path_granule_SRE_B11, 0, 0)
    print(path_cloud_mask_clip)
    
    Haralick_clip = path_res+os.path.basename(Haralick[:-4])+'_clip.TIF'
    print(Haralick_clip)
    Haralick_clip = clip_to_tile(Haralick, path_granule_SRE_B11, -10000, 0)
    print(Haralick_clip)
    
    #filename of the output
    CLD_MASK_VENUS = path_res+os.path.basename(path_granule_SRE_B11[:-4])+'_CLOUD_MASK_VENUS.TIF'
    if not os.path.isfile(CLD_MASK_VENUS): 
        #create the cloud mask with CLOUD CNES and NDVI
        mask = os.path.join('python '+path_gdal_calc_beta+'gdal_calc_modified.py --overwrite --type=Float32 -A '+NDVI+' -B '+path_cloud_mask_clip+
                        ' -C '+Haralick_clip+' --outfile='+CLD_MASK_VENUS+
                        ' --calc="numpy.logical_and(numpy.logical_and(numpy.logical_and(A>=-0.06, A<=0.05),B==128),C>0.8)" --CalcAsDT --NoDataValue=0')
        #print(mask)
        os.system(mask)
    
    CLD_MASK_VENUS_8bits = path_res+os.path.basename(path_granule_SRE_B11[:-4])+'_CLOUD_MASK_VENUS_8bits.TIF'
    if not os.path.isfile(CLD_MASK_VENUS_8bits):
        to_byte = os.path.join(path_gdal+'gdal_translate -ot Byte -a_nodata 0 '+CLD_MASK_VENUS+' '+CLD_MASK_VENUS_8bits)
        #print(to_byte)
        os.system(to_byte)
        
        
    os.remove(path_granule_sensor_array_copy2)
    os.remove(Haralick_clip)
    
    return str(CLD_MASK_VENUS)
 
#-------------------------------------------------------------------------------------------------------------
#---------------------------------------------TO CREATE THE LAKE MASK----------------------------------------
#------------------------------------------------------------------------------------------------------------ 

def lakes_mask(rasterin, path_granule_SRE_B11):
    print('Apply the lake mask to '+os.path.basename(rasterin))
    lakes = path_general+'QGIS/shp/2010_Glacial_lake_utm45_3'
    #filename of the vector layer used to crop
    croplake = path_res+os.path.basename(lakes)+'_crop'
    if not os.path.isfile(croplake):
        crop_with_lakes = os.path.join(path_gdal+'ogr2ogr -overwrite -f "ESRI Shapefile" -clipsrc '+granule_extent_shapefile_25m+'shp '+croplake+'.shp '+lakes+'.shp')
        #print(crop_with_lakes)
        os.system(crop_with_lakes)
    
    lakes_raster = path_res+os.path.basename(croplake)+'_raster.TIF'
    if not os.path.isfile(lakes_raster):
        data = gdal.Open(rasterin)
        geoTransform = data.GetGeoTransform()
        minx = geoTransform[0]
        maxy = geoTransform[3]
        maxx = minx + geoTransform[1] * data.RasterXSize
        miny = maxy + geoTransform[5] * data.RasterYSize
        #print('------------------------------------------------------------------------------------')
        #print [minx, miny, maxx, maxy]
        data = None
                
        rasterize_lakes = os.path.join(path_gdal+'gdal_rasterize -i -burn 1.0 -tr 5.0 5.0 -a_nodata 20000 -te '+str(minx)+' '+str(miny)+' '+str(maxx)+' '+str(maxy)+' -ot Float32 -of GTiff '
                                       +croplake+'.shp '+lakes_raster)   
        #print(rasterize_lakes)
        os.system(rasterize_lakes)
    
    lakes_raster_clip = clip_to_tile(lakes_raster, path_granule_SRE_B11, 20000,0)
    
    array_ones = create_matrix_of_ones(rasterin)
    
    raster_without_lakes = rasterin[:-4]+'_no_lakes.TIF'
    if not os.path.isfile(raster_without_lakes):
        suppress_lakes = os.path.join('python '+path_gdal_calc_beta+'gdal_calc_modified.py --overwrite --type=float32 -A '+rasterin+' -B '+lakes_raster_clip+
                        ' -C '+array_ones+' --outfile='+raster_without_lakes+
                        ' --calc="numpy.multiply(A,C, out = numpy.zeros(A.shape), where= B==1)" --CalcAsDT --NoDataValue=0')
        #print suppress_lakes
        os.system(suppress_lakes)
        
    os.remove(lakes_raster)
    os.remove(croplake+'.shp')
    os.remove(croplake+'.prj')
    os.remove(croplake+'.dbf')
    os.remove(croplake+'.shx')
 
    os.remove(lakes_raster_clip)
    os.remove(array_ones)
        
    return str(raster_without_lakes)
    
        
#-------------------------------------------------------------------------------------------------------------
#-----------------------------------------TO PROCESS THE COSINE CORRECTION------------------------------------
#-------------------------------------------------------------------------------------------------------------

def TopCor_Cosine(path_granule, path_granule_FRE, Zs,As, path_input_DEM, path_granule_SRE_B11):
    print("\n4-Start processing the cosine correction for "+os.path.basename(path_granule))
    
    # Make the directory to store the results 
    if not os.path.exists(path_res):
        os.makedirs(path_res)
        
    # Set the name of the final product (band_TopCor)
    TopCor_fn = path_res+os.path.basename(path_granule[:-4])+'_TopCor.TIF'   
    if os.path.exists(TopCor_fn):
        print('----TopCor corrections already exist for this image')
        return
        
    # Set names for slope and aspect of the input DEM
    TopCor_material_Fn = path_res+'TopCor_Material'
    slope_tmp_fn_sensor = TopCor_material_Fn+'/'+os.path.basename(path_input_DEM[:-4])+'_slope_sensor'+'.TIF'
    aspect_tmp_fn_sensor = TopCor_material_Fn+'/'+os.path.basename(path_input_DEM[:-4])+'_aspect_sensor'+'.TIF'
    
    # Make the directory to store the TopCor materials (slope, aspect, DEM ...)
    if not os.path.exists(TopCor_material_Fn):
        os.makedirs(TopCor_material_Fn)
        
    # set the fn of the futur clipped DEM to the tile extent of the desired sensor
    path_DEM_sensor = TopCor_material_Fn+'/'+os.path.basename(path_input_DEM[:-4])+'.TIF'
    
    #check if the slope and the aspect data exist or not    
    if not os.path.isfile(slope_tmp_fn_sensor) and not os.path.isfile(aspect_tmp_fn_sensor):         
        print('----Clip DEM Slope and Aspect to sensor extent')
                         
        print("----Reproj and clip DEM")
        # reproj and clip the input DEM to the extent of the tile (do that once per tile) using the granule_extent_shapefile.shp        
        reproj_clip_input_to_tile = os.path.join('gdalwarp -overwrite --config GDALWARP_IGNORE_BAD_CUTLINE YES -crop_to_cutline -cutline '
                 +granule_extent_shapefile_25m+'shp '+'-r '+Rmethod_DEM+' -tr '+str(SRE_B11_pixelwidth)+' '+str(SRE_B11_pixelheight)+' -s_srs '+SRSinput_DEM+' -t_srs '+SRS_SRE_B11+
                 ' '+path_input_DEM+' '+path_DEM_sensor)
        os.system(reproj_clip_input_to_tile)
        
        slope_tmp_fn_sensor, aspect_tmp_fn_sensor = slope_aspect(path_DEM_sensor)
        
                 
    # Define outputs filenames    

    path_granule_SRE_clip = path_res+os.path.basename(path_granule[:-4])+'_clip.TIF'
    path_granule_SRE_clip_norm = path_res+os.path.basename(path_granule[:-4])+'_clip_norm.TIF'
    
    #check if the SRE sensor data exist or not    
    if not os.path.isfile(path_granule_SRE_clip) and not os.path.isfile(path_granule_SRE_clip_norm): 
        path_granule_SRE_clip = clip_to_tile_and_normalized(path_granule, path_granule)
        
    path_granule_sensor_array_copy2 = path_res+os.path.basename(path_granule_SRE_clip[:-4])+'_array2.tif'  
    #check if matrix of ones exist or not
    if not os.path.isfile(path_granule_sensor_array_copy2):        
        path_granule_sensor_array_copy2 = create_matrix_of_ones(path_granule_SRE_clip)
        
            
    granule_neg_a_zero_sensor = path_res+os.path.basename(path_granule[:-4])+'_neg_a_zero_sensor.tif'
    #check if the SRE sensor with negative values to zero data exist or not    
    if not os.path.isfile(granule_neg_a_zero_sensor): 
        #Put SRE negative values to zero
        print("----Put SRE negative values to zero")
        calc_corr1 = os.path.join('python '+path_gdal_calc_beta+'gdal_calc_modified.py --overwrite --type=float32 -A '+path_granule_SRE_clip+' -B'+path_granule_sensor_array_copy2+
                          ' --outfile='+granule_neg_a_zero_sensor+' --calc="numpy.multiply(A,B, out=numpy.zeros(A.shape), where = A>=0)" --CalcAsDT --NoDataValue=20000')
        os.system(calc_corr1)
        
    path_granule_FRE_clip = path_res+os.path.basename(path_granule_FRE[:-4])+'_clip.TIF'
    path_granule_FRE_clip_norm = path_res+os.path.basename(path_granule_FRE[:-4])+'_clip_norm.TIF'
        
    #check if the FRE sensor data exist or not    
    if not os.path.isfile(path_granule_FRE_clip) and not os.path.isfile(path_granule_FRE_clip_norm) :
        path_granule_FRE_clip = clip_to_tile_and_normalized(path_granule_FRE, path_granule_FRE)

    # Define outputs filenames
    cos_Zs = path_res+os.path.basename(path_granule[:-4])+'_cos_Zs.TIF'
    cos_IL = path_res+os.path.basename(path_granule[:-4])+'_cos_IL.TIF'
    test_cos_IL = path_res+os.path.basename(path_granule[:-4])+'_cos_IL_02.TIF'
    
    if not os.path.isfile(cos_Zs):
        print("----Create cos(Zs) array")   
        calc_corr1 = os.path.join('python '+path_gdal_calc_beta+'gdal_calc_modified.py --overwrite --type=float32 -A '+path_granule_sensor_array_copy2+
                          ' --outfile='+cos_Zs+' --calc="A * cos(radians('+str(Zs)+'))" --CalcAsDT --NoDataValue=20000')
        os.system(calc_corr1)
        
    if not os.path.isfile(cos_IL):
        print("----Create cos(IL) array")
        calc_corr2 = os.path.join('python '+path_gdal_calc_beta+'gdal_calc_modified.py --overwrite --type=float32 -A '+slope_tmp_fn_sensor+
                          ' -B '+aspect_tmp_fn_sensor+' --outfile='+cos_IL+' --calc="cos(radians('+str(Zs)+'))*cos(radians(A)) + sin(radians('+str(Zs)+'))*sin(radians(A))*cos(radians('+str(As)+
                            ' - B))" --CalcAsDT --NoDataValue=20000')
        #print(calc_corr2)
        os.system(calc_corr2)
            
    if not os.path.isfile(test_cos_IL):
        print("----Create test cos(IL)>0.2 -- denominator")
        calc_corr5 = os.path.join('python '+path_gdal_calc_beta+'gdal_calc_modified.py --overwrite --type=float32 -A '+cos_IL+' --outfile='
                                 +test_cos_IL+' --calc="A*(A>=0.2)" --CalcAsDT --NoDataValue=0')
        os.system(calc_corr5)

    print("----Create TopCor")        
    calc_corr3 = os.path.join('python '+path_gdal_calc_beta+'gdal_calc_modified.py --overwrite --type=float32 -A '+granule_neg_a_zero_sensor+
            ' -B '+cos_Zs+' -C '+test_cos_IL+' -D '+slope_tmp_fn_sensor+' --outfile='+TopCor_fn+
            ' --calc="numpy.divide(numpy.divide(numpy.multiply(A,B, out = numpy.zeros(A.shape, dtype=float), where= D<>20000),C, out = numpy.zeros(A.shape, dtype = float), where=numpy.logical_and(C!=0,D<>20000)),1000)" --CalcAsDT --NoDataValue=20000')    
    #print(calc_corr3)
    os.system(calc_corr3)
    
    

    # Define filenames of the FRE and SRE products without the nodata of the TopCor
    SRE_mask_nodata = path_res+os.path.basename(path_granule[:-4])+'_mask_nodata.tif'
    FRE_mask_nodata = path_res+os.path.basename(path_granule_FRE[:-4])+'_mask_nodata.tif'
    SRE_mask_nodata_norm = path_res+os.path.basename(path_granule[:-4])+'_mask_nodata_norm.tif'
    FRE_mask_nodata_norm = path_res+os.path.basename(path_granule_FRE[:-4])+'_mask_nodata_norm.tif'
    
    print("----Create SRE with a mask of nodata")
    calc_corr3 = os.path.join('python '+path_gdal_calc_beta+'gdal_calc_modified.py --overwrite --type=float32 -A '+TopCor_fn+' -B '+path_granule_SRE_clip+
                              ' --outfile='+SRE_mask_nodata+' --calc="numpy.divide(numpy.multiply(B,numpy.ones(B.shape), out = numpy.zeros(B.shape, dtype=float), where = A<>20000), 1000)" --CalcAsDT --NoDataValue=20000')
    os.system(calc_corr3)

    print("----Create FRE with a mask of nodata")
    calc_corr3 = os.path.join('python '+path_gdal_calc_beta+'gdal_calc_modified.py --overwrite --type=float32 -A '+TopCor_fn+' -B '+path_granule_FRE_clip+
                              ' --outfile='+FRE_mask_nodata+' --calc="numpy.divide(numpy.multiply(B,numpy.ones(B.shape), out = numpy.zeros(B.shape, dtype=float), where = A<>20000), 1000)" --CalcAsDT --NoDataValue=20000')
    os.system(calc_corr3)
    
#    print("Create SRE with a mask of nodata and the SRE normalized")
#    calc_corr3 = os.path.join('python '+path_gdal_calc_beta+'gdal_calc_modified.py --overwrite --type=float32 -A '+TopCor_fn+' -B '+path_granule_SRE_clip_norm+
#                              ' --outfile='+SRE_mask_nodata_norm+' --calc="numpy.multiply(B,numpy.ones(B.shape), out = numpy.zeros(B.shape, dtype=float), where = A<>20000)" --CalcAsDT --NoDataValue=20000')
#    os.system(calc_corr3)
#
#    print("Create FRE with a mask of nodata and the FRE normalized")
#    calc_corr3 = os.path.join('python '+path_gdal_calc_beta+'gdal_calc_modified.py --overwrite --type=float32 -A '+TopCor_fn+' -B '+path_granule_FRE_clip_norm+
#                              ' --outfile='+FRE_mask_nodata_norm+' --calc="numpy.multiply(B,numpy.ones(B.shape), out = numpy.zeros(B.shape, dtype=float), where = A<>20000)" --CalcAsDT --NoDataValue=20000')
#    os.system(calc_corr3)
    

    # Test if TopCor has been created
    if os.path.isfile(TopCor_fn):                    
        #TopCor_fn_norm = path_res+os.path.basename(path_granule[:-4])+'_TopCor_norm.TIF'
        #normalisation(TopCor_fn, TopCor_fn_norm, 0.0, 10000.0)
        print('====TopCor correction finished for : '+os.path.basename(path_granule))
        os.remove(test_cos_IL)
        os.remove(cos_Zs)
        os.remove(path_granule_sensor_array_copy2)
        os.remove(granule_neg_a_zero_sensor)
        os.remove(path_granule_FRE_clip)
        os.remove(path_granule_SRE_clip)
        
        os.remove(cos_IL)
        
        
        return str(TopCor_fn)        
    else:
        print("!!!! We got a problem, try to look for dimension issues !!!!")
        return
    


def apply_cld_mask_m(rasterin, path_granule_CLDCOV, path_granule_SRE_B11):
    print('Apply cloud mask to '+os.path.basename(rasterin))
    cld_raster_clip = clip_to_tile(path_granule_CLDCOV, path_granule_SRE_B11, 0,0)

    CLDCOV_alpha = path_granule_CLDCOV[:-4]+'_alpha.tif'
    if not os.path.isfile(CLDCOV_alpha):
        remove_nodata = os.path.join(path_gdal+'gdalwarp -overwrite -srcnodata 0 -dstalpha '+cld_raster_clip+' '+CLDCOV_alpha) #work for Windows
        print(remove_nodata)
        os.system(remove_nodata)
    
    CLDCOV_alpha_b = path_granule_CLDCOV[:-4]+'_alpha_b.tif'
    if not os.path.isfile(CLDCOV_alpha_b):    
        trslt = os.path.join(path_gdal+'gdal_translate -a_srs EPSG:32645 '+CLDCOV_alpha+' '+CLDCOV_alpha_b+' -a_nodata 255 -b 2')
        print(trslt)
        os.system(trslt)
        
    edit = os.path.join('python '+path_gdal+'gdal_edit.py -unsetnodata '+CLDCOV_alpha_b)
    os.system(edit)

    raster_without_cld = rasterin[:-4]+'_no_cld.TIF'   
    if not os.path.isfile(raster_without_cld):        
        suppress_cld = os.path.join('python '+path_gdal_calc_beta+'gdal_calc_modified.py --overwrite --type=float32 -A '+rasterin+' -B '+CLDCOV_alpha_b+
                        ' --outfile='+raster_without_cld+
                        ' --calc="(A*(B<1))+((B>0)*128)" --CalcAsDT --NoDataValue=20000')
        print(suppress_cld)
        os.system(suppress_cld)
    
    os.remove(CLDCOV_alpha)
    os.remove(cld_raster_clip)
        
    return raster_without_cld



def SCA_maker(NDVI, Energy, path_granule_CLDCOV, path_granule_SRE_B11):
    
    print('Create the SCA map')
    Energy_clip= clip_to_tile(Energy, path_granule_SRE_B11, -10000, 0)
    
    CLDCOV_alpha_b = path_granule_CLDCOV[:-4]+'_alpha_b.tif'
        
    SCA_no_cld = path_res+os.path.basename(path_granule_SRE_B11[:-4])+'_SCA_no_cld.TIF'
    if not os.path.isfile(SCA_no_cld):
        mask = os.path.join('python '+path_gdal_calc_beta+'gdal_calc_modified.py --overwrite --type=float32 -A '+NDVI+' -B '+Energy_clip+' -C '+CLDCOV_alpha_b+' --outfile='+SCA_no_cld+
                        ' --calc="numpy.logical_and(numpy.logical_and(numpy.logical_and(A>=-0.16, A<=-0.02),B>0.8), C==0)" --CalcAsDT --NoDataValue=20000')
        print(mask)
        os.system(mask)
    
    SCA_8bits_no_cld = path_res+os.path.basename(path_granule_SRE_B11[:-4])+'_SCA_8bits_no_cld.TIF'
    if not os.path.isfile(SCA_8bits_no_cld):
        to_byte = os.path.join(path_gdal+'gdal_translate -ot Byte -a_nodata 255 '+SCA_no_cld+' '+SCA_8bits_no_cld)
        print(to_byte)
        os.system(to_byte)
        
    edit = os.path.join('python '+path_gdal+'gdal_edit.py -unsetnodata '+SCA_8bits_no_cld)
    print(edit)
    os.system(edit)
    edit = os.path.join('python '+path_gdal+'gdal_edit.py -unsetnodata '+CLDCOV_alpha_b)
    print(edit)
    os.system(edit)
    
    SCA_no_cld_fn = path_res+os.path.basename(path_granule_SRE_B11[:-4])+'_SCA_no_cld_fn.TIF'
    if not os.path.isfile(SCA_no_cld_fn):
        mask = os.path.join('python '+path_gdal_calc_beta+'gdal_calc_modified.py --overwrite --type=float32 -A '+SCA_8bits_no_cld+' -B '+CLDCOV_alpha_b+' --outfile='+SCA_no_cld_fn+
                        ' --calc="(A*(B<1))+((B>0)*128)" --CalcAsDT --NoDataValue=255')
        print(mask)
        os.system(mask)
    
    
    os.remove(CLDCOV_alpha_b)
    os.remove(SCA_8bits_no_cld)
    os.remove(SCA_no_cld)
    os.remove(Energy_clip)    
    return str(SCA_no_cld_fn)

def clip_SCA_to_watershed(SCA, path_granule_SRE_B11, srcnodata, dstnodata):
    dingboche = path_general+'QGIS/shp/ws_Dingboche.'
    pheriche = path_general+'QGIS/shp/ws_Pheriche.'    
    SCA_Dingboche_clip = path_res+os.path.basename(SCA[:-4])+'_Dingboche_clip.TIF'
    SCA_Pheriche_clip = path_res+os.path.basename(SCA[:-4])+'_Pheriche_clip.TIF'
    
    if not os.path.exists(SCA_Dingboche_clip):
    
        # First, open the granule and grab projection, Xres ...
        S_pixelwidth, S_pixelheight, SRSsensor = getprojection(path_granule_SRE_B11)
        
        print("Clip of "+os.path.basename(SCA))
        clip_SCA = os.path.join(path_gdal+'gdalwarp -overwrite --config GDALWARP_IGNORE_BAD_CUTLINE YES -crop_to_cutline -cutline '
                     +dingboche+'shp '+'-r '+Rmethod+' -tr '+str(S_pixelwidth)+' '+str(S_pixelheight)+' -s_srs '+SRSinput_DEM+' -t_srs '
                                        +SRSsensor+' -srcnodata '+str(srcnodata)+' -dstnodata '+str(dstnodata)+' '+SCA+' '+SCA_Dingboche_clip)
        #print(clip_SCA)
        os.system(clip_SCA)
        
    if not os.path.exists(SCA_Pheriche_clip):
    
        # First, open the granule and grab projection, Xres ...
        S_pixelwidth, S_pixelheight, SRSsensor = getprojection(path_granule_SRE_B11)
        
        print("Clip of "+os.path.basename(SCA))
        clip_SCA = os.path.join(path_gdal+'gdalwarp -overwrite --config GDALWARP_IGNORE_BAD_CUTLINE YES -crop_to_cutline -cutline '
                     +pheriche+'shp '+'-r '+Rmethod+' -tr '+str(S_pixelwidth)+' '+str(S_pixelheight)+' -s_srs '+SRSinput_DEM+' -t_srs '
                                        +SRSsensor+' -srcnodata '+str(srcnodata)+' -dstnodata '+str(dstnodata)+' '+SCA+' '+SCA_Pheriche_clip)
        #print(clip_SCA)
        os.system(clip_SCA)
        
    return str(SCA_Dingboche_clip), str(SCA_Pheriche_clip)






        
  

      
#=======================================================================================================================
#====================================================PROCESS============================================================
#=======================================================================================================================
    
print("================================================================================")
print("===================================MAIN=========================================")
print("================================================================================\n")


## Provides path for the rest of the code
sensor = 'Venus'
path_general = '/home/bessinz/'
path_images  = path_general+'IMG/' 
path_input_DEM = path_general+"MNT/MNT_HMA_clip_Venus_UTM45_5m.tif"
path_script = "/home/bessinz/path_script/"
path_gdal = '/home/bessinz/anaconda2/envs/zoe/bin/'
Rmethod_DEM = 'cubic'
Rmethod ='near'

# Path to the modified gdal_calc that add the option CalcAsDt. 
#You can find and download it on https://trac.osgeo.org/gdal/ticket/4064
path_gdal_calc_beta = path_script


path_granule_SRE, path_granule_FRE_B11, path_granule_SRE_B7, path_granule_FRE_B7, path_granule_CLM, path_granule_Energy, path_granule_CLDCOV = sensor_data_path_provider(sensor, path_images)
nb_images = len(path_granule_SRE)



DEM_pixelwidth, DEM_pixelheight, SRSinput_DEM = getprojection(path_input_DEM)    

#for i in range(nb_images):


for i in range(142,238): 
    
#---Take an image from the list of Venus images and define all the dependent paths
    pg_SRE_B11 = path_granule_SRE[i] #SRE_B11
    pg_FRE_B11 = path_granule_FRE_B11[i] #FRE_B11
    pg_SRE_B7 =path_granule_SRE_B7[i] #SRE_B7
    pg_FRE_B7 = path_granule_FRE_B7[i] #FRE_B7
    pg_CLM = path_granule_CLM[i] #CLM
    pg_energy = path_granule_Energy[i] #Energy
    pg_CLDCOV = path_granule_CLDCOV[i] #CLOUDCOV manually corrected


#---get projection of the choosen image
    SRE_B11_pixelwidth, SRE_B11_pixelheight, SRS_SRE_B11 = getprojection(pg_SRE_B11)
    
#---Create the path to save all the results    
    path_res = path_images+os.path.basename(os.path.dirname(pg_SRE_B11))+'/Results_T2/'
    # Make the directory to store the results 
    if not os.path.exists(path_res):
        os.makedirs(path_res)

#---Create the path to save slope and aspect clipped
    TopCor_material_Fn = path_res+'TopCor_Material'        
    # Make the directory to store the TopCor materials (slope, aspect, DEM ...)
    if not os.path.exists(TopCor_material_Fn):
        os.makedirs(TopCor_material_Fn)

#---Create the shapefile using the nodata. This layer will be used to crop all the rasters depending of the image    
    granule_extent_shapefile_25m = create_shp_tile(pg_SRE_B11, '-10000')
    
#---Get the sun angles from the metadata file
    Zs, As = findGranuleMetadata(sensor, pg_SRE_B11)  

#================================================================================

#---Clip the FRE and SRE images to the tile with the shapefile created before
    print('\n3-Clip SRE an FRE to the useful extent of the tile for the NIR band')
    pg_FRE_B11_clip = clip_to_tile_and_normalized(pg_FRE_B11, pg_FRE_B11)  
    pg_SRE_B11_clip = clip_to_tile_and_normalized(pg_SRE_B11, pg_SRE_B11)
    
#---Process the cosine correction on the NIR band (11)            
    TopCor_fn_B11 = TopCor_Cosine(pg_SRE_B11, pg_FRE_B11, Zs,As, path_input_DEM, pg_SRE_B11)
    
#================================================================================

#---Clip the FRE and SRE images to the tile with the shapefile created before   
    print('\n5-Clip SRE an FRE to the useful extent of the tile for the Red band') 
    pg_FRE_B7_clip = clip_to_tile_and_normalized(pg_FRE_B7, pg_FRE_B11)     
    pg_SRE_B7_clip = clip_to_tile_and_normalized(pg_SRE_B7, pg_SRE_B11) 

#---Process the cosine correction on the red band (7)               
    TopCor_fn_B7 = TopCor_Cosine(pg_SRE_B7, pg_FRE_B7, Zs,As, path_input_DEM, pg_SRE_B11)

#================================================================================

    TopCor_fn_B11 = path_res+os.path.basename(pg_SRE_B11[:-4])+'_TopCor.TIF'
    TopCor_fn_B7 = path_res+os.path.basename(pg_SRE_B7[:-4])+'_TopCor.TIF'
   
#---Calculate the NDVI with NIR and red bands corrected
    print('\n7-Calculate NDVI on Red and NIR bands corrected')
    NDVI_TopCor = calcul_NDVI(TopCor_fn_B11, TopCor_fn_B7)
    
#---Put the lakes to nodata in the NDVI layer
    print('\n8-Apply lakes mask on the NDVI raster')
    NDVI_TopCor_no_lakes = lakes_mask(NDVI_TopCor, pg_SRE_B11)
        
#---Put the clouds areas to nodata in NIR TopCor
    print('\n10-Apply cloud mask to the NIR band corrected')
    TopCor_B11_no_cld = apply_cld_mask_m(TopCor_fn_B11, pg_CLDCOV, pg_SRE_B11)

#---Create SCA map using NDVI, Energy matrix and the cloud mask we produced    
    print('\n11-Create the SCA map')
    SCA_8bits_no_cld_fn = SCA_maker(NDVI_TopCor, pg_energy, pg_CLDCOV, pg_SRE_B11)
    
    print('\n12-Clip SCA map to the two watersheds')
    SCA_Dingboche_no_cld_fn, SCA_Pheriche_no_cld_fn = clip_SCA_to_watershed(SCA_8bits_no_cld_fn, pg_SRE_B11, 255, 255)
    
    print('\n13-Remove all useless data')
    os.remove(NDVI_TopCor)
    os.remove(NDVI_TopCor_no_lakes)
    os.remove(TopCor_fn_B7)
    
    
    slope_tmp_fn_sensor = TopCor_material_Fn+'/'+os.path.basename(path_input_DEM[:-4])+'_slope_sensor.TIF'
    aspect_tmp_fn_sensor = TopCor_material_Fn+'/'+os.path.basename(path_input_DEM[:-4])+'_aspect_sensor'+'.TIF'
    path_DEM_sensor = TopCor_material_Fn+'/'+os.path.basename(path_input_DEM[:-4])+'.TIF'
     
    os.remove(slope_tmp_fn_sensor)
    os.remove(aspect_tmp_fn_sensor)
    os.remove(path_DEM_sensor)
    
    cos_IL_B7 = path_res+os.path.basename(pg_SRE_B7[:-4])+'_cos_IL.TIF'
    SRE_mask_nodata_B7 = path_res+os.path.basename(pg_SRE_B7[:-4])+'_mask_nodata.tif'
    FRE_mask_nodata_B7 = path_res+os.path.basename(pg_FRE_B7[:-4])+'_mask_nodata.tif'
    SRE_mask_nodata_B11 = path_res+os.path.basename(pg_SRE_B11[:-4])+'_mask_nodata.tif'
    FRE_mask_nodata_B11 = path_res+os.path.basename(pg_FRE_B11[:-4])+'_mask_nodata.tif'

    os.remove(SRE_mask_nodata_B7)
    os.remove(FRE_mask_nodata_B7)
    os.remove(granule_extent_shapefile_25m+"dbf")
    os.remove(granule_extent_shapefile_25m+"prj")
    os.remove(granule_extent_shapefile_25m+"shp")
    os.remove(granule_extent_shapefile_25m+"shx")
    os.remove(TopCor_fn_B11)
    os.remove(SRE_mask_nodata_B11)
    os.remove(FRE_mask_nodata_B11)

 


    







    

    







