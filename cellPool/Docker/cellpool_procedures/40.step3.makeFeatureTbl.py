import os
import numpy as np
import pandas as pd
import sys
from PIL import Image
Image.MAX_IMAGE_PIXELS = 1000000000 #Turn off DecompressionBombWarning for large mosaics
import imageio
from skimage.segmentation import expand_labels
from skimage.measure import regionprops_table
import json

cyc = sys.argv[1]
well = sys.argv[2]
paramsFile = sys.argv[3]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.utils import *
from cellpool.mosaic import *

wellBorderWidth = params['featureExtractionParams']['wellBorderWidth']
mosaicDir = os.path.join(params['analysisDir'], params['mosaicDir'])
maskDir = os.path.join(params['analysisDir'], params['segmentationDir'], \
	'iter' + str(params['segmentationParams']['segIter']) + \
	params['segmentationParams']['maskDirSuffix'])
properties = tuple(params['featureExtractionParams']['regionProperties'])

thisMaskExt = os.path.splitext( \
	[x for x in os.listdir(maskDir) if x.startswith(well + \
		'-cyc' + str(cyc) + '_' + params['segmentationParams']['compartment'] + \
		'_cp_masks')][0])[1]
im = imageio.imread(os.path.join(maskDir, \
	well + '-cyc' + cyc + '_' + params['segmentationParams']['compartment'] + \
	'_cp_masks' + thisMaskExt))
im = im > 0

breakX = [params['confluencyQcGridSize']*x for x in \
	range(np.floor(im.shape[1]/params['confluencyQcGridSize']).astype(int) + 1)]
breakX.append(im.shape[1])

breakY = [params['confluencyQcGridSize']*x for x in \
	range(np.floor(im.shape[0]/params['confluencyQcGridSize']).astype(int) + 1)]
breakY.append(im.shape[0])

confluencyEsts = {"gridXStart" : [], "gridXEnd" : [], \
	"gridYStart" : [], "gridYEnd" : [], "confluency" : []}
for x in range(len(breakX) - 1):
	for y in range(len(breakY) - 1):
		confluencyEsts["gridXStart"].append(breakX[x])
		confluencyEsts["gridXEnd"].append(breakX[x + 1])
		confluencyEsts["gridYStart"].append(breakY[y])
		confluencyEsts["gridYEnd"].append(breakY[y + 1])
		confluencyEsts["confluency"].append(round(np.sum(im[breakY[y]:breakY[y + 1], \
			breakX[x]:breakX[x + 1]])/(params['confluencyQcGridSize'] ** 2), 2))

confluencyEsts = pd.DataFrame.from_dict(confluencyEsts)
confluencyEsts.to_csv(well + '-cyc' + cyc + '.confluencyEsts.txt', \
	sep = '\t', index = None)

#----------------------------------------------------

def borderToucher(regionmask, intensity):
	this = intensity[regionmask]
	return any(np.isnan(this))

#----------------------------------------------------
chInfo = pd.read_csv(os.path.join(params['analysisDir'], params['chToStainId']), sep = '\t')
nucCh = chInfo[chInfo.cycle == 'cycle' + str(cyc)].T
nucCh = nucCh[nucCh == \
	params['segmentationParams']['cellposeParams']['useStain']].dropna().index[0][-1] + \
	'-cyc' + str(cyc)

chIm = np.load(os.path.join(mosaicDir, well + '-C' + nucCh + '.mosaic.trimmed.npy'))
thisMaskExt = os.path.splitext( \
	[x for x in os.listdir(maskDir) if x.startswith(well) and \
		('-cyc' + str(cyc) in x) and \
		('_cp_masks' in x)][0])[1]
mask = imageio.imread(os.path.join(maskDir, well + '-cyc' + cyc + '_' + \
			params['segmentationParams']['compartment'] + '_cp_masks' + thisMaskExt))

borderIm = imageio.imread(os.path.join(mosaicDir, well + '-cyc' + str(cyc) + \
	'.trimmed.tileOverlapZone.png'))

#------------------------------------
#Get features
props = regionprops_table(mask, intensity_image = chIm, \
			properties = properties)
props2 = regionprops_table(mask, intensity_image = borderIm, \
	properties = ('label', 'intensity_max'))
props3 = regionprops_table(expand_labels(mask, \
				distance = wellBorderWidth).astype(int), \
		intensity_image = np.where(np.isnan(chIm), np.nan, mask), \
		properties = ('label', ), \
		extra_properties = (borderToucher, ))

props['isTouchingTileBorder'] = props2['intensity_max'] > 0
props['isTouchingWellBorder'] = props3['borderToucher']
keysOfProps = [x for x in props]
outTbl = pd.DataFrame({key: value for key, value \
			in props.items() \
			if key in keysOfProps})

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

outTbl = outTbl.merge(confluencyEsts)

outTbl.to_csv(os.path.join(well + '-cyc' + str(cyc) + '.' + \
	params['segmentationParams']['cellposeParams']['useStain'].replace(' ', '_') + '.txt'), \
	sep = '\t', index = False)






