import pandas as pd
import numpy as np
import os
from skimage.io import imread, imsave
from skimage.exposure import rescale_intensity

well = 'A08'
im1 = np.load('../30.mosaics/' + well + '-C2-cyc1.mosaic.trimmed.npy')
im2 = np.load('../30.mosaics/' + well + '-C4-cyc2.mosaic.trimmed.npy')

xstart = 4000
xend = 4400
ystart = 5000
yend = 5400

im1x = im1[ystart:yend, xstart:xend].astype(np.uint16)
im1x = rescale_intensity(im1x)
imsave("2e_1.png", im1x)
im2x = im2[ystart:yend, xstart:xend].astype(np.uint16)
im2x = rescale_intensity(im2x)
imsave("2e_2.png", im2x)











