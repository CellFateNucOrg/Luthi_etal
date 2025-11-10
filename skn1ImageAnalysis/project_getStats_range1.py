#!/usr/bin/env python
# coding: utf-8

# In[9]:


#from skimage.io import imread, imsave
from microfilm.microplot import microshow
from nd2reader import ND2Reader
import numpy as np
import re
import os
import glob
import nd2
import random
from fnmatch import fnmatch
import matplotlib.pyplot as plt
import pyclesperanto_prototype as cle
#from stardist.models import StarDist2D
#from stardist import random_label_cmap, _draw_polygons, export_imagej_rois
from csbdeep.utils import normalize
from stardist import random_label_cmap
import stackview
from skimage import measure, morphology, filters, exposure
import pandas as pd 
from napari_segment_blobs_and_things_with_membranes import threshold_otsu
import pyclesperanto_prototype as cle
from tifffile import imread, imwrite
np.random.seed(6)
from PIL import Image
lbl_cmap = random_label_cmap()


# 

# In[10]:


def filter_labels_by_size(label_image,min_size,max_size):
    out = np.copy(label_image)
    component_sizes = np.bincount(label_image.ravel())
    too_small = component_sizes < min_size
    too_small_mask = too_small[label_image]
    out[too_small_mask] = 0
    too_big = component_sizes > max_size
    too_big_mask = too_big[label_image]
    out[too_big_mask] = 0
    return(out)


def get_label_props(label_image,intensity_image):
    properties = measure.regionprops_table(label_image=label_image, intensity_image=np.round(intensity_image), 
                            properties=('label','area', 'intensity_mean', 'intensity_max', 'intensity_min'))
    data = pd.DataFrame(properties)
    return(data)


def rescaleImageIntensity(image,min_val,max_val):
    # rescale to 0-1 range
    normalized_image = ((image.astype(np.float32) - min_val)/(max_val-min_val))
    # stretch scale
    rescaled_image = exposure.rescale_intensity(normalized_image, in_range=(0, 1))
    rescaled_image_uint16 = (rescaled_image * 65535).astype(np.uint16)
    return(rescaled_image_uint16)

def makeComposite(green,brightfield):
    # Normalize to the range 0-1
    green_rescaled = (green - green.min()) / (green.ptp() / 255.0)
    brightfield_rescaled = (brightfield - brightfield.min()) / (brightfield.ptp() / 255.0)
    allGreen = green_rescaled + brightfield_rescaled
    allGreen[allGreen>1] = 1
    # convert to 8-bit
    green_uint8 = green_rescaled.astype(np.uint8)
    #green_uint8 = allGreen.astype(np.uint8)
    brightfield_uint8 = brightfield_rescaled.astype(np.uint8)
    # Convert the numpy array to a PIL image
    green_uint8 = Image.fromarray(green_uint8)
    brightfield_uint8 = Image.fromarray(brightfield_uint8)
    composite = Image.merge('RGB',(brightfield_uint8, green_uint8, brightfield_uint8))
    return(composite)

def plot_labels_and_image(intensity_image,label_image,bf_image,title=""):
    binary_borders = cle.detect_label_edges(label_image)
    labeled_borders = binary_borders * label_image
    axis_norm = (0,1)
    plt.figure(figsize=(15,8))
    plt.subplot(131); 
    plt.imshow(intensity_image, cmap='gray')
    plt.axis('off')
    plt.title("GFP")
    plt.subplot(132); 
    plt.imshow(label_image, cmap=lbl_cmap)
    plt.axis('off')
    plt.title("ASI labels")
    plt.subplot(133); 
    plt.imshow(bf_image,vmin=np.min(bf_image),vmax=np.max(bf_image)) 
    plt.axis('off')
    plt.title("composite")
    plt.tight_layout()
    figure = plt.gcf()
    return figure


# In[11]:


work_dir = "/mnt/external.data/MeisterLab/pmeister/skn-1_deletions/240411_skn-1_deletions"
#work_dir = '/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/microscopy/20240411_skn-1_deletions'

lowerQ=0.5
upperQ=99.95
lowerVal = 142
upperVal = 434

useQuantiles = True

os.makedirs(work_dir+'/imageStats', exist_ok=True)  
os.makedirs(work_dir+'/maxProj', exist_ok=True)
os.makedirs(work_dir+'/rescaled_bf', exist_ok=True)

os.makedirs(work_dir+'/rescaled_gfp_'+str(lowerQ)+'_'+str(upperQ), exist_ok=True)
os.makedirs(work_dir+'/composite_'+str(lowerQ)+'_'+str(upperQ), exist_ok=True)
os.makedirs(work_dir+'/blind_gfp_'+str(lowerQ)+'_'+str(upperQ), exist_ok=True)
os.makedirs(work_dir+'/labels_'+str(lowerQ)+'_'+str(upperQ), exist_ok=True)

if not useQuantiles:
    os.makedirs(work_dir+'/rescaled_gfp_'+str(lowerVal)+'_'+str(upperVal), exist_ok=True)
    os.makedirs(work_dir+'/composite_'+str(lowerVal)+'_'+str(upperVal), exist_ok=True)
    os.makedirs(work_dir+'/blind_gfp_'+str(lowerVal)+'_'+str(upperVal), exist_ok=True)
    os.makedirs(work_dir+'/labels_'+str(lowerVal)+'_'+str(upperVal), exist_ok=True)


# ## Get intensity values from control images

# In[12]:


pattern=work_dir+'/**/1116*.nd2'
file_list=glob.glob(pattern,recursive=True)
#print(len(file_list))
def getControlImageProps(file_list,bf_loQ,bf_hiQ,gfp_loQ,gfp_hiQ):
    df = pd.DataFrame(index=range(0,len(file_list)),
    columns=['base_name','max_slice','max_slice_final','gfp_mean','gfp_median','gfp_loQ','gfp_hiQ',
    'otsu_threshold','mean_blob','min_blob','max_blob','bf_loQ','bf_hiQ'])

    #file_name='/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/microscopy/20240411_skn-1_deletions/1116_HS/1116_HS_001.nd2'
    #file_name='/mnt/external.data/MeisterLab/pmeister/skn-1_deletions/240411_skn-1_deletions/1116_HS/1116_HS_001.nd2'
    i=0
    for file_name in file_list:
        base_name=str(re.sub('.nd2','', os.path.basename(file_name)))
        df.loc[i,"base_name"] = base_name
        
        # read in raw image and make max projection
        image = nd2.imread(file_name)
        channel1 = image[:,0,:,:]
        channel2 = image[:,1,:,:]
        c1 = np.max(channel1,axis=0)

        # choose brightfield slice to use
        highSigPerSlice = [np.percentile(channel1[slice,:,:],99.9) for slice in range(0,channel1.shape[0])]
        maxSlice = np.argmax(highSigPerSlice)
        df.loc[i,"max_slice"] = maxSlice
        
        if(maxSlice<6 or maxSlice>15):
            maxSlice = 10
        df.loc[i,"max_slice_final"] = maxSlice
        c2 = channel2[maxSlice,:,:]
        df.loc[i,"bf_loQ"] = np.percentile(c2,bf_loQ)
        df.loc[i,"bf_hiQ"] = np.percentile(c2,bf_hiQ)

        #get image stats
        denoised_median = filters.median(c1, morphology.disk(1))
        df.loc[i,"gfp_mean"] = np.mean(denoised_median)
        df.loc[i,"gfp_median"] = np.median(denoised_median)
        df.loc[i,"gfp_loQ"] = np.percentile(denoised_median,gfp_loQ)
        df.loc[i,"gfp_hiQ"] = np.percentile(denoised_median,gfp_hiQ)

        # get asi stats
        asi_props = pd.DataFrame()
        threshold = filters.threshold_otsu(denoised_median)
        df.loc[i,"otsu_threshold"] = threshold
        binary_c1 = denoised_median >= threshold
        labeled = measure.label(binary_c1, connectivity=1)
        filtered_labels = filter_labels_by_size(labeled,min_size=100,max_size=3000)
        tmp = get_label_props(filtered_labels,c1)
        if(asi_props.empty):
            asi_props = tmp
        else:
            asi_props = pd.concat([asi_props,tmp])
        col_means=tmp.mean()
        df.loc[i,"mean_blob"] = col_means["intensity_mean"]
        df.loc[i,"min_blob"] = col_means["intensity_min"]
        df.loc[i,"max_blob"] = col_means["intensity_max"]
        i = i+1

    return(df)

df =  getControlImageProps(file_list, bf_loQ=0.01, bf_hiQ=99.99, gfp_loQ=lowerQ, gfp_hiQ=upperQ)


if (useQuantiles): 
    min_val = df.loc[:, 'gfp_loQ'].mean()
    max_val = df.loc[:, 'gfp_hiQ'].mean()
else:
    min_val = lowerVal
    max_val = upperVal

min_val_bf = df.loc[:, 'bf_loQ'].mean()
max_val_bf = df.loc[:, 'bf_hiQ'].mean()

df.to_csv(path_or_buf=work_dir+'/imageStats/controlImageStats_'+str(lowerQ)+'_'+str(upperQ)+'_'+str(int(min_val))+'_'+str(int(max_val))+'.csv')

threshold_df = pd.DataFrame.from_dict({'min_val': np.round([min_val]),'max_val': np.round([max_val]), 'min_val_bf': np.round([min_val_bf]), 'max_val_bf': np.round([max_val_bf])})
threshold_df.to_csv(path_or_buf=work_dir+'/imageStats/thresholds.csv')
threshold_df


# In[13]:


threshold_df = pd.read_csv(work_dir+'/imageStats/thresholds.csv')
min_val = threshold_df.loc[0,'min_val']
max_val = threshold_df.loc[0,'max_val']
min_val_bf = threshold_df.loc[0,'min_val_bf']
max_val_bf = threshold_df.loc[0,'max_val_bf']


# ## Process all images

# In[14]:


threshold_df = pd.read_csv(work_dir+'/imageStats/thresholds.csv')
min_val = threshold_df.loc[0,'min_val']
max_val = threshold_df.loc[0,'max_val']
min_val_bf = threshold_df.loc[0,'min_val_bf']
max_val_bf = threshold_df.loc[0,'max_val_bf']

pattern=work_dir+'/**/*.nd2'
file_list=glob.glob(pattern,recursive=True)
min_blob_size = 500
max_blob_size = 6000

df1 = pd.DataFrame(index=range(0,len(file_list)),
   columns=['blind_key','base_name','max_slice','max_slice_final','gfp_mean','gfp_median','gfp_loQ','gfp_hiQ',
   'otsu_threshold','mean_blob','min_blob','max_blob','bf_loQ','bf_hiQ','worm_num','moved','remove_slices'])

# default values
df1.loc[:,"worm_num"] = 1
df1.loc[:,"moved"] = 0
df1.loc[:,"remove_slices"] = -1

asi_props = pd.DataFrame()

# shuffle indeces to get blind key for scoring
random.seed(33)
blind_keys = list(range(len(file_list)))
random.shuffle(blind_keys)
df1.loc[:,"blind_key"] = ['image'+str(x).zfill(3) for x in blind_keys]

#file_name='/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/microscopy/20240411_skn-1_deletions/1116_HS/1116_HS_002.nd2'
file_name='/mnt/external.data/MeisterLab/pmeister/skn-1_deletions/240411_skn-1_deletions/1116_HS/1116_HS_001.nd2'
i=0
for file_name in file_list:
    base_name=str(re.sub('.nd2','', os.path.basename(file_name)))
    df1.loc[i,"base_name"] = base_name

    # read in denoised and raw images
    image = nd2.imread(file_name)
    channel1 = image[:,0,:,:]
    channel2 = image[:,1,:,:]
    channel1.shape

    # max projection of gfp
    c1 = np.max(channel1,axis=0)
    imwrite(work_dir+'/maxProj/'+base_name+'_maxProj.tif', c1)

    # choose gfp slice with highest signal to choose brightfield slice
    highSigPerSlice = [np.percentile(channel1[slice,:,:],99.9) for slice in range(0,channel1.shape[0])]
    maxSlice = np.argmax(highSigPerSlice)
    df1.loc[i,"max_slice"] = maxSlice

    # take the middle slice if a very topmost or bottommost slices were chosen
    if(maxSlice<6 or maxSlice>15):
        maxSlice = 10
    df1.loc[i,"max_slice_final"] = maxSlice
    c2 = channel2[maxSlice,:,:]

    # get brightfield image stats
    df1.loc[i,"bf_loQ"] = np.percentile(c2,0.01)
    df1.loc[i,"bf_hiQ"] = np.percentile(c2,99.99)
    imwrite(work_dir+'/maxProj/'+base_name+'_bf.tif', c2)

    #get gfp image stats
    df1.loc[i,"gfp_mean"] = np.mean(c1)
    df1.loc[i,"gfp_median"] = np.median(c1)
    df1.loc[i,"gfp_loQ"] = np.percentile(c1,lowerQ)
    df1.loc[i,"gfp_hiQ"] = np.percentile(c1,upperQ)

    # detect ASI 
    denoised_median = filters.median(c1, morphology.disk(1))
    threshold = filters.threshold_otsu(denoised_median)
    df1.loc[i,"otsu_threshold"] = threshold
    binary_c1 = denoised_median >= threshold
    labeled = measure.label(binary_c1, connectivity=1)
    filtered_labels = filter_labels_by_size(labeled, min_size=min_blob_size, max_size=max_blob_size)
    tmp = get_label_props(filtered_labels,c1)
    tmp['base_name'] = base_name

    if(asi_props.empty):
        asi_props = tmp
    else:
        asi_props = pd.concat([asi_props,tmp])

    # save average values from ASIs to main table
    col_means=tmp.iloc[:,2:5].mean()
    df1.loc[i,"mean_blob"] = col_means["intensity_mean"]
    df1.loc[i,"min_blob"] = col_means["intensity_min"]
    df1.loc[i,"max_blob"] = col_means["intensity_max"]

    if not useQuantiles:
        lq = lowerQ
        uq = upperQ
        lowerQ = min_val
        upperQ = max_val

    # save stats to csv
    df1.to_csv(path_or_buf=work_dir+'/imageStats/allImageStats_'+str(lowerQ)+'_'+str(upperQ)+'_'+str(int(min_val))+'_'+str(int(max_val))+'.csv')
    asi_props.to_csv(path_or_buf=work_dir+'/imageStats/blobStats_'+str(lowerQ)+'_'+str(upperQ)+'_'+str(int(min_val))+'_'+str(int(max_val))+'.csv')

    # rescale gfp image
    rescaled = rescaleImageIntensity(c1,min_val,max_val)
    imwrite(work_dir+'/rescaled_gfp_'+str(lowerQ)+'_'+str(upperQ)+'/'+base_name+'.tif', rescaled)

    # save gfp image with random number for blind scoring
    imwrite(work_dir+'/blind_gfp_'+str(lowerQ)+'_'+str(upperQ)+'/'+df1.loc[i,"blind_key"]+'.tif', rescaled)

    # rescale brightfield image
    rescaled_bf = rescaleImageIntensity(c2,min_val_bf,max_val_bf)
    imwrite(work_dir+'/rescaled_bf/'+base_name+'_bf.tif', rescaled_bf)

    # make composite image
    composite = makeComposite(rescaled,rescaled_bf)
    composite.save(work_dir+'/composite_'+str(lowerQ)+'_'+str(upperQ)+'/'+base_name+'_composite.tif')

    plot_labels_and_image(rescaled,filtered_labels,composite,title=base_name)
    plt.savefig(work_dir+'/labels_'+str(lowerQ)+'_'+str(upperQ)+'/'+base_name+'_labels.png')
    plt.close()

    if not useQuantiles:
        lowerQ = lq
        upperQ = uq

    i = i+1

 


# 

# #plot_labels_on_image(labeled,c1,title=base_name)
# #microshow(labeled)
# cle.detect_label_edges(labeled)
