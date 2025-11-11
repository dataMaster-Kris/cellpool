import os
import numpy as np
import pandas as pd
import sys
from PIL import Image
Image.MAX_IMAGE_PIXELS = 1000000000 #Turn off DecompressionBombWarning for large mosaics
import imageio
from skimage.io import imread
from skimage.measure import regionprops_table
import json
import tifffile

wellOmeTiff = sys.argv[1]
paramsFile = sys.argv[2]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.utils import distribSummarize

maskDir = os.path.join(params['analysisDir'], params['ashlarAnalysisDir'], \
	params['segmentationDir'], \
	'iter' + str(params['segmentationParams']['segIter']) + \
	params['segmentationParams']['maskDirSuffix'])

properties = ('label', 'area', 'centroid', 'eccentricity', 'feret_diameter_max', 'orientation')
im = imread(wellOmeTiff)
im = np.transpose(im, axes = [1, 2, 0])
mask = imageio.imread(os.path.join(maskDir, wellOmeTiff.replace('.ome.tiff', '') + \
	'-cyc1_' + params['segmentationParams']['compartment'] + '_cp_masks.png')) 
props = regionprops_table(mask, intensity_image = im, \
			properties = properties, \
			extra_properties = (distribSummarize,))
keysOfProps = [x for x in props if x != 'distribSummarize']
for ch in range(im.shape[2]):
	outTbl = pd.concat([pd.DataFrame({key: value for key, value \
				in props.items() \
				if key in keysOfProps}), \
			pd.DataFrame([x[ch] for x in props['distribSummarize'].tolist()])], \
			axis = 1)
	outTbl.to_csv(path_or_buf = wellOmeTiff.replace('.ome.tiff', '') + \
			'.ch' + str(ch) +  '.csv', index = None)

