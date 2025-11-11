import os
import numpy as np
import sys
import json
import pandas as pd
from collections import OrderedDict
from skimage.segmentation import mark_boundaries
from PIL import Image, ImageDraw, ImageFont
from tifffile import imread

well = sys.argv[1]
smpl = sys.argv[2]
paramsFile = sys.argv[3]
with open(paramsFile, 'r') as inFile:
        params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.unpool import *

smplDir = params['unpoolTrackingQcParams']['showUnpoolDir']
dpWidth = params['unpoolTrackingQcParams']['dpWidth'] #Width of region to extract around each object
pctl4Normalization = 'q_' + str(params['unpoolTrackingQcParams']['dpNormPctl'])
nSamples = 50 
nPerCol = 10 
nCol = 5
sepWidth = params['unpoolTrackingQcParams']['sepWidth']
textHeight = params['unpoolTrackingQcParams']['textBoxHeight']
font = '/app/UbuntuMono-R.ttf'
fontsize = int(textHeight*params['unpoolTrackingQcParams']['fontSizeWrtBox'])
iteration = params['segmentationParams']['segIter']

mosaicDir_cp = os.path.join(params['analysisDir'], params['mosaicDir'])
maskDir_cp = os.path.join(params['analysisDir'], params['segmentationDir'], \
	'iter' + str(iteration) + params['segmentationParams']['maskDirSuffix'])
normalizers = pd.read_csv(os.path.join(mosaicDir_cp, \
	"wholePlate.rawIntensityQntls.txt"), sep = "\t")
normalizers = dict(zip(normalizers[['cycle', 'ch']].apply( \
		lambda x: 'C' + str(x[1]) + '-cyc' + str(x[0]), axis = 1).tolist(), \
		normalizers[pctl4Normalization].tolist()))

#-------------------------------------------
#Load tracking mates for this well and sample a set from this.
trackDir = os.path.join(params['analysisDir'], params['ashlarAnalysisDir'], \
	'45.track_cellpool_to_ashlar')
trackFileSuffix = '_tracking' + params['trackingParams']['stg2OutSuffix']
trackTbl = pd.DataFrame()
trackTbl = pd.read_csv(os.path.join(trackDir, well + '.cellpool_ashlar' + \
	trackFileSuffix), sep = '\t')

#-----------------------------------------------
#Draw a random sample
trackTbl = trackTbl.sample(n = nSamples)
cols = trackTbl.columns
trackTbl = trackTbl.apply(lambda x: [well + 'n' + str(y) for y in x], \
		axis = 1, result_type = 'expand')
trackTbl.columns = cols
for cyc in cols:
	trackTbl[cyc].to_csv(well + '-' + cyc + '.' + smpl + '.txt', \
		header = None, index = None)
#------------------------------------------------
mosaicSuffix = '.mosaic.trimmed.npy'
chsOfInterest = []
allIds = pd.DataFrame()
stitchingRef = 'HOECHST 33342' 

cycle = 1
chInfo = pd.read_csv(os.path.join(params['analysisDir'], \
		params['chToStainId']), sep = '\t')
mainCh = chInfo[chInfo.cycle == 'cycle' + str(cycle)].T
mainCh = mainCh[mainCh == stitchingRef].dropna().index[0][-1] + \
	'-cyc' + str(cycle)
chsOfInterest.append('C' + mainCh + mosaicSuffix)
for tl in ['ashlar', 'cellpool']:
	allIds[tl] = pd.read_csv( \
		well + '-' + tl + '.' + smpl + '.txt', header = None)[0].tolist()

thisMaskExt = os.path.splitext( \
        [x for x in os.listdir(maskDir_cp) if x.startswith(well) and \
                ('-cyc1' in x) and ('_cp_masks' in x)][0])[1]
ims_cp = unpoolObjects( \
	fileWithListOfObjIds = well + '-' + 'cellpool' + '.'  + smpl + '.txt', \
	dirWithImages = mosaicDir_cp, keepImIfPrefix = well, \
	keepImIfSuffix = mosaicSuffix, \
	objMaskFilePath = maskDir_cp, \
	objMaskId = well + '-cyc1_' + \
		params['segmentationParams']['compartment'] + '_cp_masks' + thisMaskExt, \
	returnAllMasksInDir = False, returnIndObjMask = True, showOnlyTgtObjs = False, \
	dirWithMasks = maskDir_cp, \
	keepMaskIfPrefix = well, dpWidth = dpWidth)

mosaicDir_ashlr = os.path.join(params['analysisDir'], params['ashlarAnalysisDir'], "20.mosaics")
im_ashlr = imread(os.path.join(mosaicDir_ashlr, well + '_' + \
	params['mosaickingParams']['ashlarParamSuffix'] + '.ome.tiff'))
with open(well + '-ashlar.npy', 'wb') as f:
	np.save(f, 65535*im_ashlr[0, :, :])

maskDir_ashlr = os.path.join(params['analysisDir'], params['ashlarAnalysisDir'], \
	params['segmentationDir'], \
	'iter' + str(iteration) + params['segmentationParams']['maskDirSuffix'])
ims_ashlr = unpoolObjects( \
	fileWithListOfObjIds = well + '-' + 'ashlar' + '.'  + smpl + '.txt', \
	dirWithImages = '.', keepImIfPrefix = well, \
	keepImIfSuffix = '-ashlar.npy', \
	objMaskFilePath = maskDir_ashlr, \
	objMaskId = well + '_' + params['mosaickingParams']['ashlarParamSuffix'] + '-cyc1_' + \
		params['segmentationParams']['compartment'] + '_cp_masks' + thisMaskExt, \
	returnAllMasksInDir = False, returnIndObjMask = True, showOnlyTgtObjs = False, \
	dirWithMasks = maskDir_ashlr, \
	keepMaskIfPrefix = well, dpWidth = dpWidth)

ims = {'images': {}, 'mainMasks' : ims_cp['mainMasks']}
ims['images']['cellpool'] = ims_cp['images']['C' + mainCh + '.mosaic.trimmed.npy']
ims['images']['ashlar'] = dict(zip(ims['images']['cellpool'].keys(),
		ims_ashlr['images']['ashlar.npy'].values()))
ims['masks'] = {}
cyc1_ids = allIds['cellpool'].tolist()
for cyc in allIds.columns:
	ims['masks'][cyc] = {}
	for obj in allIds['cellpool'].tolist():
		thisObj = allIds[cyc].loc[allIds.cellpool == obj].iloc[0]
		if (cyc == 'cellpool'):
			ims['masks'][cyc][obj] = ims_cp['mainMasks'][obj]
		else:
			ims['masks'][cyc][obj] = ims_ashlr['mainMasks'][thisObj]

normIms = OrderedDict({k1: OrderedDict({k2:np.where(v2 <= normalizers['C' + mainCh], v2, \
				normalizers['C' + mainCh])/normalizers['C' + mainCh] for \
	k2, v2 in v1.items()}) for \
	k1, v1 in ims['images'].items() if k1 in ['cellpool', 'ashlar']})

nStains = 2
outmage_chs = 3
outmage = np.zeros((nSamples * dpWidth, nStains * dpWidth, outmage_chs))
chsOfInterest = ['cellpool', 'ashlar']
for obj in range(nSamples):
	for stn in range(nStains):
		thisOriShape = normIms[chsOfInterest[stn]][cyc1_ids[obj]].shape
		thisImage = normIms[chsOfInterest[stn]][cyc1_ids[obj]].copy()
		thisCycObjId = allIds[chsOfInterest[stn]].loc[ \
			allIds.cellpool == cyc1_ids[obj]].iloc[0]
		thisImage = mark_boundaries(thisImage, \
			ims['masks'][chsOfInterest[stn]][cyc1_ids[obj]], \
				color = (1, 0, 0), mode = 'thick')
		outmage[((dpWidth * obj) + int((dpWidth - thisOriShape[0])/2)): \
			((dpWidth * obj) + int((dpWidth - thisOriShape[0])/2) + \
			thisOriShape[0]), ((dpWidth * stn) + \
			int((dpWidth - thisOriShape[1])/2)):((dpWidth * stn) + \
			int((dpWidth - thisOriShape[1])/2) + thisOriShape[1]), :] = \
				thisImage 

outmage[[x*dpWidth for x in range(nSamples)], :, :] = 1
outmage[:, [x*dpWidth for x in range(nStains)], :] = 1
outmage = (outmage*255).astype(np.uint8)

outmage2 = outmage[:(nPerCol*dpWidth), :, :]
for col in range(2, nCol + 1):
	outmage2 = np.concatenate((outmage2, np.full((nPerCol*dpWidth, sepWidth, 3), \
			fill_value = 255).astype(np.uint8), \
			outmage[(nPerCol*dpWidth*(col - 1)):(nPerCol*dpWidth*col), :, :]), \
			axis = 1)

outmage2 = np.concatenate((outmage2, np.full((textHeight, outmage2.shape[1], 3), \
			fill_value = 255).astype(np.uint8)), axis = 0)

listOfWells = pd.read_csv(os.path.join(params['analysisDir'], params['listOfWells']), \
		header = None)[0].tolist()
wellOrder = listOfWells.index(well)

outmage2 = Image.fromarray(outmage2)
write = ImageDraw.Draw(outmage2)
myFont = ImageFont.truetype(font, size = fontsize)
write.text((sepWidth, outmage2.size[1] - fontsize - 1), \
	str(nSamples) + ' random samples of tracking mates' + \
	(' out of ' + str(nLowConf) + ' ' if smpl == 'LowConf' else ' ') + 'in ' + \
	well, font=myFont, fill = (0,))

outmage2.save(os.path.join(f"{wellOrder:02d}" + '.' + smpl + '.png'))

