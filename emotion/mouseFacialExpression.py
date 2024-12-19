# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 19:45:34 2021

@author: pc
"""

# %%

import pandas as pd
import seaborn as sns
import numpy as np  
import matplotlib
import matplotlib.pyplot as plt 
from skimage import io
from skimage.feature import hog
from joblib import Parallel, delayed
# import cv2
import glob
import matplotlib.image as mpimg
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import scipy.io as sio

# %%
###example image - full frame
image = mpimg.imread("C:/Users/yihudu/Desktop/Yihui/data/m1271/frames_downsample_3/lick (10001).jpg")

plt.imshow(image, cmap="gray")

# %%
###example image - cropped
image_crop = image[80:350, 120:630]

width = 510
height = 270
dpi = 300
plt.figure(figsize=(width/dpi, height/dpi), dpi=dpi)
    
plt.subplots_adjust(top=1,bottom=0,left=0,right=1,hspace =0, wspace =0)
plt.margins(0,0)
    
plt.imshow(image_crop, cmap="gray")  

plt.savefig("C:/Users/yihudu/Desktop/Yihui/data/m1271/frames_downsample_3/lick_crop (10001).png")
# plt.close()     

# %%
###Data loading and preprocessing
#the imagesToHogsCellCrop() function will take all image files in a given folder, 
#crop the area as given in the cropCoords 
#and convert them into their HOG (histogram of oriented gradients) descriptors.
#pixels_per_cell argument defines the sliding window size for HOG creation.
# n_jobs argument defines number of threads used

def imagesToHogsCellCrop(imgFolderAndFileType, pixelsPerCell):
    from skimage import io
    from skimage.feature import hog
    from joblib import Parallel, delayed
    
    coll = io.ImageCollection(imgFolderAndFileType)
    
    kwargs = dict(orientations=8, pixels_per_cell=(pixelsPerCell, pixelsPerCell), cells_per_block=(1, 1), transform_sqrt=True, channel_axis=-1)
    # kwargs = dict(orientations=9, pixels_per_cell=(pixelsPerCell, pixelsPerCell), cells_per_block=(2, 2), transform_sqrt=True, multichannel=True)
    #kwargs = dict(orientations=8, pixels_per_cell=(pixelsPerCell, pixelsPerCell), cells_per_block=(1, 1), visualize=True, transform_sqrt=True)
    
    return Parallel(n_jobs=32)(delayed(hog)(image[80:350, 120:630], **kwargs) for image in coll)
    # return Parallel(n_jobs=32)(delayed(hog)(image, **kwargs) for image in coll)
    
# %%    
#replace the filepath below for own use
fear = imagesToHogsCellCrop("C:/Users/yihudu/Desktop/Yihui/data/m1271/frames_downsample_3/*.jpg", 8)
fear = pd.DataFrame(data = fear)

# %%
# ###Paiwrise similarity

# #calculate similarities (correlation coefficients) between all possible pairs of frames
corr_fear = fear.iloc[np.r_[0:5245]].T.corr()   

# %%
writer = pd.ExcelWriter('C:/Users/yihudu/Desktop/Yihui/data/m1271/frames_downsample_3/lick_8.xlsx')
corr_fear.to_excel(writer)
# writer.book.use_zip64()
writer.close()
# writer.save()

# #clear repetitive pattern in the pairwise correlation matrix can be observed (follow the diagonal)
sns.heatmap(corr_fear, robust = True, rasterized = True)

# %%
# ###Paiwrise similarity - hiearachical clustering

# #hierarchically clustered pairwise correlation matrix
g = sns.clustermap(corr_fear, robust=True, rasterized = True)

# %%
#post-hoc labelling the frames based on their temporal origin (gray = before each stimulus, red = during each stimulus) 
# colorVec1 = (["very light green"] * 500 + ["kiwi"] * 625 + ["light green"] * 1375 + ["green"] * 375 + ["darkgreen"] * 813
#              + ["light cyan"] * 375 + ["aqua"] * 625 + ["sky blue"] * 500 + ["dodger blue"] * 312 + ["blue"] * 500 + ["dark blue"] * 438)
# colorVec1 = (["yellow"] * 894 + ["orange"] * 894 + ["red"] * 894 + ["dark coral"] * 894 + ["reddy brown"] * 894)
# colorVec1 = (["light khaki"] * 2000 + ["green"] * 1937 + ["darkgreen"] * 2562)
# colorVec1 = (["yellow"] * 1031 + ["kiwi"] * 1000 + ["orange"] * 750 + ["green"] * 969 + ["red"] * 656 + ["darkgreen"] * 1281)
# colorVec1 = (["yellow"] * 625 + ["light khaki"] * 1125 + ["orange"] * 459 + ["kiwi"] * 1084 + ["red"] * 500 + ["green"] * 1834 + ["reddy brown"] * 417 + ["darkgreen"] * 938)
# colorVec1 = (["light cyan"] * 894 + ["aqua"] * 894 + ["sky blue"] * 894 + ["dodger blue"] * 894 + ["blue"] * 894 + ["dark blue"] * 894)
# colorVec1 = (["egg shell"] * 480 + ["very light green"] * 250 
             # + ["yellow"] * 292 + ["kiwi"] * 313 
             # + ["orange"] * 688 + ["light green"] * 688 
             # + ["red"] * 167 + ["green"] * 188 
             # + ["dark coral"] * 209 + ["darkgreen"] * 407 
             # + ["egg shell"] * 292 + ["light cyan"] * 188
             # + ["yellow"] * 167 + ["aqua"] * 313 
             # + ["orange"] * 188 + ["sky blue"] * 250 
             # + ["red"] * 375 + ["dodger blue"] * 156 
             # + ["dark coral"] * 292 + ["blue"] * 250 
             # + ["reddy brown"] * 334 + ["dark blue"] * 219)
# colorVec1 = (["kiwi"] * 4100 + ["green"] * 2227)
# colorVec1 = (["yellow"] * 1500 + ["kiwi"] * 3750 + ["orange"] * 1139 + ["green"] * 1114)
# colorVec1 = sns.xkcd_palette(colorVec1)
# sns.palplot(colorVec1)
# colorVec1 = (["yellow"] * 500 + ["blue"] * 750 + ["yellow"] * 1000 + ["green"] * 3000 + ["blue"] * 1375 + ["yellow"] * 875)
# colorVec1 = (["yellow"] * 812 + ["green"] * 2750 + ["blue"] * 438)
colorVec1 = (["lightgray"] * 37 + ["firebrick"] * 456 + ["lightgray"] * 287 + ["blue"] * 195 + ["firebrick"] * 812 + ["lightgray"] * 253 + ["blue"] * 288
             + ["firebrick"] * 223 + ["lightgray"] * 559 + ["blue"] * 249 + ["firebrick"] * 547 + ["lightgray"] * 278 + ["blue"] * 298 + ["firebrick"] * 509 + ["lightgray"] * 254)
g = sns.clustermap(corr_fear, robust=True, col_colors = colorVec1, row_colors = colorVec1, rasterized = True)


# %%
Folder = 'C:/Users/yihudu/Desktop/Yihui/data/m1271/'

col_linkage = g.dendrogram_col.linkage
# row_linkage = g.dendrogram_row.linkage
col_indices = g.dendrogram_col.reordered_ind
# row_indices = g.dendrogram_row.reordered_ind


# plt.figure()
width = 900
height = 300
dpi = 300
plt.figure(figsize=(width/dpi, height/dpi), dpi=dpi)
matplotlib.rcParams['lines.linewidth'] = 0.5

# cluster_indices = fcluster(col_linkage, t=0.7*max(col_linkage[:,2]), criterion='distance')

# dendrogram(col_linkage, no_labels=True, color_threshold=None)
# plt.axhline(y=0.7*max(col_linkage[:,2]), c="k", linestyle='--') 

cluster_indices = fcluster(col_linkage, t=2.1, criterion='distance')

dendrogram(col_linkage, no_labels=True, color_threshold=2.1)
plt.axhline(y=2.1, c="k", linestyle='--') 

plt.yticks(np.arange(0, max(col_linkage[:,2]), step=1))
plt.yticks(fontsize=6)
plt.savefig(Folder + 'cluster_map_8.png')
plt.close()  
  

matfn = Folder + 'cluster_indices_8.mat'
sio.savemat(matfn,{'cluster_indices': cluster_indices})

matfn = Folder + 'cluster_all_8.mat'
sio.savemat(matfn,{'cluster_indices': cluster_indices,'col_linkage': col_linkage,'col_indices': col_indices})

# %%
###Create neutral prototypical face

neutralFolder = "D:/duyh/video/20210628 optogenetics_video/CT (Baseline)/sleep_trial5/nrem/*.jpg"

protoAll_neutral = imagesToHogsCellCrop(neutralFolder, 8)
protoAll_neutral = pd.DataFrame(data = protoAll_neutral)
proto_neutral = protoAll_neutral.mean(axis = 0)

# %%
####Creating an emotion prototype - fear

from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
t1 = fear.corrwith(proto_neutral, axis = 1) #comparing this single tail shock presentation to the neutral prototypical face we created above
plt.plot(scaler.fit_transform(t1.to_numpy().reshape(-1, 1)), color = "gray")

# %%
proto_fear = fear.iloc[t1.nsmallest(10).index].mean(axis = 0) #10 frames most dissimilar to the neutral face are averaged and used to form the prototypical pain face

t2 = fear.corrwith(proto_fear, axis = 1)
plt.plot(scaler.fit_transform(t2.to_numpy().reshape(-1, 1)), color = "firebrick")

# %%
plt.plot(scaler.fit_transform(t1.to_numpy().reshape(-1, 1)) - 1, color = "gray") #neutral face similarity
plt.plot(scaler.fit_transform(t2.to_numpy().reshape(-1, 1)), color = "firebrick") #pain face similarity. They follow almost, but not completely, inverted paths

