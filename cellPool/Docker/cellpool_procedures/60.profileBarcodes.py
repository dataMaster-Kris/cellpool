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

well = sys.argv[1]
paramsFile = sys.argv[2]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.utils import distribSummarize

mosaicDir = os.path.join(params['analysisDir'], params['mosaicDir'])
maskDir = os.path.join(params['analysisDir'], params['segmentationDir'], \
	'iter' + str(params['segmentationParams']['segIter']) + \
	params['segmentationParams']['maskDirSuffix'])

properties = ('label', 'area', 'centroid')

mosaicSuffix = '.mosaic.trimmed.npy'
cycles = ['cycle' + str(x) for x in params['cycles']]
for cyc in params['cycles']:
	thisMaskExt = os.path.splitext( \
				[x for x in os.listdir(maskDir) if x.startswith(well) and \
					('-cyc' + str(cyc) in x) and \
					('_cp_masks' in x)][0])[1]
	mask = imageio.imread(os.path.join(maskDir, well + '-cyc' + str(cyc) + '_' + \
		params['segmentationParams']['compartment'] + '_cp_masks' + thisMaskExt))
	chImFiles = [x for x in os.listdir(mosaicDir) if x.endswith(mosaicSuffix) \
			and x.startswith(well) and ('-cyc' + str(cyc) in x)]
	chIm = {}
	for nextFile in chImFiles:
		chIm[nextFile.replace(mosaicSuffix, '')] = \
			np.load(os.path.join(mosaicDir, nextFile))
	#------------------------------------
	#Get features
	for cycleCh in chIm:
		props = regionprops_table(mask, intensity_image = chIm[cycleCh], \
					properties = properties, \
					extra_properties = (distribSummarize,))
		keysOfProps = [x for x in props if x != 'distribSummarize']
		outTbl = pd.concat([pd.DataFrame({key: value for key, value \
					in props.items() \
					if key in keysOfProps}), \
				pd.DataFrame(props['distribSummarize'].tolist())], axis = 1)
		outTbl['gridXStart'] = [int(x/params["confluencyQcGridSize"]) * \
				params["confluencyQcGridSize"] \
				for x in outTbl['centroid-1']]
		outTbl['gridXEnd'] = [min(mask.shape[1], \
				(1 + int(x/params["confluencyQcGridSize"])) * \
				params["confluencyQcGridSize"]) \
				for x in outTbl['centroid-1']]
		outTbl['gridYStart'] = [int(x/params["confluencyQcGridSize"]) * \
				params["confluencyQcGridSize"] \
				for x in outTbl['centroid-0']]
		outTbl['gridYEnd'] = [min(mask.shape[0], \
				(1 + int(x/params["confluencyQcGridSize"])) * \
				params["confluencyQcGridSize"]) \
				for x in outTbl['centroid-0']]
		outTbl.to_csv(path_or_buf = cycleCh +  '.csv', index = None)



