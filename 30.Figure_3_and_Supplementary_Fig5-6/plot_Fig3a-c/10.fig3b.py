import numpy as np
import os
import pandas as pd
from skimage.io import imsave
from skimage import exposure, draw
import matplotlib.pyplot as plt

well = "E10"
objId = 12200 
imFiles = os.listdir('../../30.mosaics') 
imFiles = [x for x in imFiles if x.startswith(well) & x.endswith('.trimmed.npy') & \
        (('cyc1' in x) | ('cyc2' in x))]

im = {}
for fl in imFiles:
    im[fl[:11]] = np.load('../../30.mosaics/' + fl)
    im[fl[:11]] = np.nan_to_num(im[fl[:11]])

ftrFiles = os.listdir('../../40.segmentation/iter1.step3.features/')
ftrFiles = [x for x in ftrFiles if x.startswith(well) & x.endswith('.HOECHST_33342.txt')]

ftr = {}
for fl in ftrFiles:
    ftr[fl[:8]] = pd.read_csv('../../40.segmentation/iter1.step3.features/' + fl, sep = '\t')

r = int(ftr[well + '-cyc1'][ftr[well + '-cyc1'].label == objId]['centroid-0'].iloc[0])
c = int(ftr[well + '-cyc1'][ftr[well + '-cyc1'].label == objId]['centroid-1'].iloc[0])

saturated_pixels = 0.3 #Percent of saturated pixels in the image
height = 940
width = 600

#----------------------
#cycle1
im11 = im[well + '-C1-cyc1'][(r - int(height/2)):(r + int(height/2)), \
        (c - int(width/2)):(c + int(width/2))]
im12 = im[well + '-C2-cyc1'][(r - int(height/2)):(r + int(height/2)), \
        (c - int(width/2)):(c + int(width/2))]
im13 = im[well + '-C3-cyc1'][(r - int(height/2)):(r + int(height/2)), \
        (c - int(width/2)):(c + int(width/2))]
im1 = np.dstack((im11, im12, im13))
enhanced_im1 = np.empty_like(im1)
for channel in range(im1.shape[2]):
    v_min, v_max = np.percentile(im1[:, :, channel], \
            (saturated_pixels / 2, 100 - saturated_pixels / 2))
    enhanced_im1[:, :, channel] = exposure.rescale_intensity(im1[:, :, channel], \
            in_range=(v_min, v_max))

x_start, y_start, x_end = 500, 900, 560
bar_length, bar_color, bar_thickness = 60, 1, 10   

for i in range(-bar_thickness//2 +1, bar_thickness//2+1):
    rr1, cc1 = draw.line(y_start+i, x_start, y_start+i, x_end)
    enhanced_im1[rr1, cc1] =  np.array(bar_color)

imsave(well + '-cyc1_' + str(r) + '_' + str(c) + '.tiff', enhanced_im1)
imsave(well + '-cyc1_' + str(r) + '_' + str(c) + '.png', \
        (enhanced_im1*255/np.max(enhanced_im1)).astype(np.uint8))

#----------------------
#cycle2
im21 = im[well + '-C1-cyc2'][(r - int(height/2)):(r + int(height/2)), \
        (c - int(width/2)):(c + int(width/2))]
im22 = im[well + '-C2-cyc2'][(r - int(height/2)):(r + int(height/2)), \
        (c - int(width/2)):(c + int(width/2))]
im23 = im[well + '-C3-cyc2'][(r - int(height/2)):(r + int(height/2)), \
        (c - int(width/2)):(c + int(width/2))]
im2 = np.dstack((im21, im22, im23))
enhanced_im2 = np.empty_like(im2)
for channel in range(im2.shape[2]):
    v_min, v_max = np.percentile(im2[:, :, channel], \
            (saturated_pixels / 2, 100 - saturated_pixels / 2))
    enhanced_im2[:, :, channel] = exposure.rescale_intensity(im2[:, :, channel], \
            in_range=(v_min, v_max))

#enhanced_im2 = exposure.adjust_gamma(enhanced_im2, 0.7)
imsave(well + '-cyc2_' + str(r) + '_' + str(c) + '.tiff', enhanced_im2)
imsave(well + '-cyc2_' + str(r) + '_' + str(c) + '.png', \
        (enhanced_im2*255/np.max(enhanced_im2)).astype(np.uint8))

#----------------------
#cycles3-5
nucInROI = ftr[well + '-cyc1'].loc[ \
        ((r - int(height/2)) < ftr[well + '-cyc1']['centroid-0']) & \
        (ftr[well + '-cyc1']['centroid-0'] < (r + int(height/2))) & \
        ((c - int(width/2)) < ftr[well + '-cyc1']['centroid-1']) & \
        (ftr[well + '-cyc1']['centroid-1'] < (c + int(width/2))), ]
intDat = pd.read_csv('../../60.barcodes/integrated_intensity_table.txt.gz', sep = '\t') 
intROI = intDat.loc[intDat.objId.isin([well + 'n' + str(x) for x in nucInROI.label])].copy()
intROI2 = intROI.filter(like = "cycle").apply(lambda x: x*100/np.median(x))
intROI.loc[:, intROI2.columns] = intROI2

sizePerBcd = 30*30
sizePerElem = int(np.sqrt(sizePerBcd)/3)
bcd = {}
bcdCh = [['cycle3.ch1', 'cycle3.ch2', 'cycle3.ch3'], \
        ['cycle4.ch1', 'cycle4.ch2', 'cycle4.ch3'], \
        ['cycle5.ch1', 'cycle5.ch2', 'cycle5.ch3']]
outIm = np.empty(im[well + '-C3-cyc1'].shape[:2])
for obj in intROI.objId:
    bcd[obj] = np.empty((int(np.sqrt(sizePerBcd)), int(np.sqrt(sizePerBcd))))
    thisDat = nucInROI.loc[nucInROI.label == int(obj.split('n')[1])]
    rstart, rend = [int(thisDat['centroid-0'].iloc[0]) - int(np.sqrt(sizePerBcd)/2), \
            int(thisDat['centroid-0'].iloc[0]) + int(np.sqrt(sizePerBcd)/2)]
    cstart, cend = [int(thisDat['centroid-1'].iloc[0]) - int(np.sqrt(sizePerBcd)/2), \
            int(thisDat['centroid-1'].iloc[0]) + int(np.sqrt(sizePerBcd)/2)]
    for ix_r in range(3):
        for ix_c in range(3):
            bcd[obj][int(ix_r*np.sqrt(sizePerBcd)/3):int((ix_r + 1)*np.sqrt(sizePerBcd)/3), \
                    int(ix_c*np.sqrt(sizePerBcd)/3):int((ix_c + 1)*np.sqrt(sizePerBcd)/3)] = \
                    intROI[bcdCh[ix_r][ix_c]].loc[intROI.objId == obj].iloc[0]
    #outIm[int(thisDat['centroid-1'].iloc[0]), int(thisDat['centroid-0'].iloc[0])] = 1000
    outIm[rstart:rend, cstart:cend] = bcd[obj]

outIm = outIm[(r - int(height/2)):(r + int(height/2)), (c - int(width/2)):(c + int(width/2))]

enhanced_outIm = np.empty_like(outIm)
v_min, v_max = np.percentile(outIm, \
        (saturated_pixels / 2, 100 - saturated_pixels / 2))
enhanced_outIm = exposure.rescale_intensity(outIm, \
        in_range=(v_min, v_max))

imsave(well + '-cyc3-5_' + str(r) + '_' + str(c) + '.tiff', 1 - enhanced_outIm)
imsave(well + '-cyc3-5_' + str(r) + '_' + str(c) + '.png', \
        ((1 - enhanced_outIm)*255/np.max((1 - enhanced_outIm))).astype(np.uint8))










