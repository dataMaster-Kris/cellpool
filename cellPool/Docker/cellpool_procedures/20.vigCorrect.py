import skimage.io
import os
import json
import numpy as np
import sys
from basicpy import BaSiC

files = sys.argv[1:]
files = [x.strip() for x in files]
uniqueChannels = np.unique([os.path.basename(x)[7:] for x in files])

for ch in uniqueChannels:
	thisFiles = [x for x in files if ch in os.path.basename(x) and x.endswith('.tiff')]
	images = []
	for x in thisFiles:
		images.append(skimage.io.imread(x))

	images = np.array(images)
	
	b = BaSiC(get_darkfield = True, fitting_mode = 'ladmap', reweighting_tol = 0.001)
	b.fit(images)
	flatfield, darkfield = b.flatfield, b.darkfield
	imCorr = [(x - darkfield)/flatfield for x in images]
	
	for ix in range(len(thisFiles)):
		skimage.io.imsave(fname = os.path.basename(thisFiles[ix])[:-5] + \
						'.vigCorr.tiff', \
						arr = imCorr[ix], check_contrast = False)
	
	skimage.io.imsave(fname = os.path.basename(thisFiles[ix])[4:-5] + \
			'.flat-field.tiff', arr = flatfield, check_contrast = False)
	skimage.io.imsave(fname = os.path.basename(thisFiles[ix])[4:-5] + \
			'.dark-field.tiff', arr = darkfield, check_contrast = False)





