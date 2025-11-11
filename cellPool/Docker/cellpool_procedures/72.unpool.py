import os
import numpy as np
import sys
import json
from skimage.io import imsave
from collections import OrderedDict
import pandas as pd

objList = sys.argv[1]
paramsFile = sys.argv[2]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.unpool import *

prefix = os.path.basename(objList).removesuffix('.toUnpool.txt')
dpWidth = params['unpoolClustParams']['dpWidth'] #Width of region to extract around each object
iteration = params['segmentationParams']['segIter']
cycle = params['unpoolClustParams']['imagingCycle']
#Font for labeling images.
fontFile = '/app/UbuntuMono-R.ttf' 

mosaicDir = os.path.join(params['analysisDir'], params['mosaicDir'])
maskDir = os.path.join(params['analysisDir'], params['segmentationDir'], \
	'iter' + str(iteration) + params['segmentationParams']['maskDirSuffix'])
vals = pd.read_csv(objList, sep = '\t', header = None)
#-------------------------------------------
#-------------------------------------------
#image showing all the nuclei within the mask and neighborhood with boundary on top
ims = unpoolObjects( \
        fileWithListOfObjIds = objList, dirWithImages = mosaicDir, keepImIfWellPrefix = True, \
        keepImIfContains = '-cyc' + str(cycle), keepImIfSuffix = '.mosaic.trimmed.npy', \
        objMaskFilePath = maskDir, objMaskIdContains = '-cyc' + str(cycle) + '_' + \
                        params['segmentationParams']['compartment'] + '_cp_masks', \
	objMaskIdWellPrefix = True, objMaskIdFindExt = True, \
        returnAllMasksInDir = False, returnIndObjMask = True, showOnlyTgtObjs = False, \
        dirWithMasks = maskDir, dpWidth = dpWidth)

objIds = [x for x in ims['images']['C1-cyc' + str(cycle) + '.mosaic.trimmed.npy']]
idsAndAnnot = [x + ' ' + str(round(vals[1][vals[0] == x].iloc[0], 3)) for \
		x in ims['images']['C1-cyc' + str(cycle) + '.mosaic.trimmed.npy']]

chInfo = pd.read_csv(os.path.join(params['analysisDir'], \
		params['chToElemIdFile']), sep = '\t')
rgbCh = chInfo[chInfo.cycle == 'cycle' + str(cycle)].T
rgbCh = ['C' + rgbCh[rgbCh == x].dropna().index[0][-1] + \
	'-cyc' + str(cycle) for x in params['unpoolClustParams']['channels']]

reds = normalizeForViewing(dict(zip(idsAndAnnot, \
	[ims['images'][rgbCh[0] + '.mosaic.trimmed.npy'][x] for \
		x in objIds])), \
	method = params['unpoolClustParams']['normMethod'], \
	limits = params['unpoolClustParams']['rescaleRangeLimits'])
greens = normalizeForViewing(dict(zip(idsAndAnnot, \
	[ims['images'][rgbCh[1] + '.mosaic.trimmed.npy'][x] for \
		x in objIds])), \
	method = params['unpoolClustParams']['normMethod'], \
	limits = params['unpoolClustParams']['rescaleRangeLimits'])
blues = normalizeForViewing(dict(zip(idsAndAnnot, \
	[ims['images'][rgbCh[2] + '.mosaic.trimmed.npy'][x] for \
		x in objIds])), \
	method = params['unpoolClustParams']['normMethod'], \
	limits = params['unpoolClustParams']['rescaleRangeLimits'])
maskImgs = dict(zip(idsAndAnnot, \
	[ims['mainMasks'][x] for x in objIds]))
outIm = showUnpooled(gray = None, \
	red = reds, \
	green = greens, \
	blue = blues, maskImgs = maskImgs, markBoundaries = True, \
	ncols = params['unpoolClustParams']['nPerFile'][1], \
	showObjIds = True, dpWidth = dpWidth)

outIm.save(prefix + '.png')

