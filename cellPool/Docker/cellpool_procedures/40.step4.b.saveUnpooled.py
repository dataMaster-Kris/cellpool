import os
import numpy as np
import sys
import json
import pandas as pd

cycle = sys.argv[1]
well = sys.argv[2]
paramsFile = sys.argv[3]
with open(paramsFile, 'r') as inFile:
        params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.unpool import *

chInfo = pd.read_csv(os.path.join(params['analysisDir'], params['chToStainId']), sep = '\t')

#Font for labeling images.
fontFile = '/app/UbuntuMono-R.ttf' 
iteration = params['segmentationParams']['segIter']
dpWidth = params['unpoolSegQcParams']['dpWidth'] #Width of region to extract around each object
mosaicDir = os.path.join(params['analysisDir'], params['mosaicDir'])
maskDir = os.path.join(params['analysisDir'], params['segmentationDir'], \
	'iter' + str(iteration) + params['segmentationParams']['maskDirSuffix'])
objList = os.path.join(params['analysisDir'], params['segmentationDir'], \
	'iter' + str(iteration) + params['segmentationParams']['unpoolDirSuffix'], \
	params['segmentationParams']['unpooledNominationsDir'], \
	well + '-cyc' + cycle + '.nominations.txt')
#-------------------------------------------
#-------------------------------------------
#image showing all the nuclei within the mask and neighborhood with boundary on top
thisMaskExt = os.path.splitext( \
	[x for x in os.listdir(maskDir) if x.startswith(well) and \
		('-cyc' + cycle in x) and ('_cp_masks' in x)][0])[1]
ims = unpoolObjects( \
	fileWithListOfObjIds = objList, dirWithImages = mosaicDir, keepImIfPrefix = well, \
	keepImIfContains = '-cyc' + cycle, keepImIfSuffix = '.mosaic.trimmed.npy', \
	objMaskFilePath = maskDir, objMaskId = well + '-cyc' + cycle + '_' + \
			params['segmentationParams']['compartment'] + '_cp_masks' + thisMaskExt, \
	returnAllMasksInDir = False, returnIndObjMask = True, showOnlyTgtObjs = False, \
	dirWithMasks = maskDir, keepMaskIfPrefix = cycle, dpWidth = dpWidth, \
	keepMaskIfSuffix = thisMaskExt)

mainCh = chInfo[chInfo.cycle == 'cycle' + cycle].T
mainCh = mainCh[mainCh == \
	params['segmentationParams']['cellposeParams']['useStain']].dropna().index[0][-1] + \
	'-cyc' + cycle

chOtherThanMain = [x for x in ims['images'].keys() if 'C' + mainCh not in x]
imCh = {'red': [x for x in ims['images'].keys() if 'C' + mainCh in x][0]}

if len(chOtherThanMain) > 0:
	imCh['green'] = chOtherThanMain[0]
else:
	imCh['green'] = imCh['red']

if len(chOtherThanMain) > 1:
	imCh['blue'] = chOtherThanMain[1]
else:
	imCh['blue'] = imCh['red']

objIds = [x for x in ims['images'][imCh['blue']]]
dpIds = [x.replace('n', '-cyc' + cycle + '_n') for x in ims['images'][imCh['blue']]]
reds = normalizeForViewing(dict(zip(dpIds, \
	[ims['images'][imCh['red']][x] for x in objIds])))
greens = normalizeForViewing(dict(zip(dpIds, \
	[ims['images'][imCh['green']][x] for x in objIds])))
blues = normalizeForViewing(dict(zip(dpIds, \
	[ims['images'][imCh['blue']][x] for x in objIds])))

#Need to update keys associated with image masks
for (objLabel, objMask), objLabelNCycle in zip(list(ims['mainMasks'].items()), dpIds):
    del ims['mainMasks'][objLabel]
    ims['mainMasks'][objLabelNCycle] = objMask

saveUnpooled(gray = None, \
	red = reds, \
	green = greens, \
	blue = blues, 
	markBoundaries = True, \
	maskImgs = ims['mainMasks'],
	ncols = 1, \
	showObjIds = True, \
	font = fontFile, dpWidth = dpWidth)

