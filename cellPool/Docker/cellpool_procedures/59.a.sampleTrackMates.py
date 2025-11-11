import os
import numpy as np
import sys
import json
import pandas as pd
from collections import OrderedDict
from skimage.segmentation import mark_boundaries
from PIL import Image, ImageDraw, ImageFont

well = sys.argv[1]
smpl = sys.argv[2]
trackMateFileSuffix = sys.argv[3]
paramsFile = sys.argv[4]
with open(paramsFile, 'r') as inFile:
        params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.unpool import *
from cellpool.tracking import *

smplDir = params['unpoolTrackingQcParams']['showUnpoolDir']
dpWidth = params['unpoolTrackingQcParams']['dpWidth'] #Width of region to extract around each object
pctl4Normalization = 'q_' + str(params['unpoolTrackingQcParams']['dpNormPctl'])
nSamples = params['unpoolTrackingQcParams']['nSamplesPerWell']
nPerCol = params['unpoolTrackingQcParams']['layout'][0] 
nCol = params['unpoolTrackingQcParams']['layout'][1]
sepWidth = params['unpoolTrackingQcParams']['sepWidth']
textHeight = params['unpoolTrackingQcParams']['textBoxHeight']
font = '/app/UbuntuMono-R.ttf'
fontsize = int(textHeight*params['unpoolTrackingQcParams']['fontSizeWrtBox'])
iteration = params['segmentationParams']['segIter']

mosaicDir = os.path.join(params['analysisDir'], params['mosaicDir'])
maskDir = os.path.join(params['analysisDir'], params['segmentationDir'], \
	'iter' + str(iteration) + params['segmentationParams']['maskDirSuffix'])
normalizers = pd.read_csv(os.path.join(mosaicDir, \
	"wholePlate.rawIntensityQntls.txt"), sep = "\t")
normalizers = dict(zip(normalizers[['cycle', 'ch']].apply( \
		lambda x: 'C' + str(x[1]) + '-cyc' + str(x[0]), axis = 1).tolist(), \
		normalizers[pctl4Normalization].tolist()))

#-------------------------------------------
#Load tracking mates for this well and sample a set from this.
trackDir = os.path.join(params['analysisDir'], params['trackObjsDir'])
trackFileSuffix = '_tracking' + params['trackingParams'][trackMateFileSuffix]
trackTbl = pd.DataFrame()
for cyc in params['registrationParams']['registrationPairs']:
	if (trackTbl.shape[0] == 0):
		trackTbl = pd.read_csv(os.path.join(trackDir, well + '.' + \
			params['registrationParams']['registrationPairs'][cyc] + \
			'_' + cyc + trackFileSuffix), sep = '\t')
	else:
		trackTbl = trackTbl.merge(pd.read_csv(
			os.path.join(trackDir, well + '.' + \
			params['registrationParams']['registrationPairs'][cyc] + '_' + cyc + \
			trackFileSuffix), sep = '\t'))

if trackTbl.shape[0] <= 1:
	print(well + " : 0 tracking mates.")
	trackTbl.to_csv('null.allLabels.txt', header = None, index = None)
	Image.fromarray(np.zeros((2, 2, 3)).astype(np.uint8)).save('null.png')
	sys.exit()

#-------------------------------------------
#Load feature table for this well and subset tracking mates if smplSet is LowConf
nLowConf = trackTbl.shape[0]
nHighConf = trackTbl.shape[0]
ngbrsToConsider = params['unpoolTrackingQcParams']['ngbrsToConsiderHighConfSearch']
approxmtInfCost = params['trackingParams']['localRefParams']['approxmtInfCost']
if (trackMateFileSuffix == 'stg2OutSuffix'):
	localNgbrAlgnCostCutoff = params['unpoolTrackingQcParams']['localNgbrAlgnCostCutoff']
elif (trackMateFileSuffix == 'stg3OutSuffix'):
	localNgbrAlgnCostCutoff = params['trackingParams']['localRefParams']['sqrDistCutoff']

if ((smpl == 'LowConf') or (smpl == 'HighConf')):
	ftrDir = os.path.join(params['analysisDir'], params['segmentationDir'], 
		'iter' + str(iteration) + params['segmentationParams']['featureDirSuffix'])
	prop = params['unpoolTrackingQcParams']['prop']
	propDiffCutoff = params['unpoolTrackingQcParams']['propDiffCutoff'][ \
					1 if smpl == 'LowConf' else 0]
	comp = trackTbl.copy()
	for cyc in params['cycles']:
		ftr = pd.read_csv(os.path.join(ftrDir, well + '-cyc' + str(cyc) + \
			'.' + params['segmentationParams']['cellposeParams'][ \
				'useStain'].replace(" ", "_") + '.txt'), \
			sep = '\t')[['label', prop, 'centroid-0', 'centroid-1']]
		ftr.columns = ['cycle' + str(cyc), prop + '_cycle' + str(cyc), \
				'centroid-0' + '_cycle' + str(cyc), \
				'centroid-1' + '_cycle' + str(cyc)]
		comp = pd.merge(comp, ftr, on = 'cycle' + str(cyc))	

	for latterCyc in params['registrationParams']['registrationPairs']:
		prevCyc = params['registrationParams']['registrationPairs'][latterCyc]
		comp[prevCyc + '_vs_' + latterCyc] = np.log(comp[prop + '_' + latterCyc]/ \
				comp[prop + '_' + prevCyc])

	if (trackMateFileSuffix == 'stg2OutSuffix'):
		highConf = comp[[x for x in comp.columns if '_vs_cycle' in x]].apply( \
			lambda x: all([abs(v) < propDiffCutoff for v in x]), axis = 1)
		if (highConf.sum() > 10):
			trackTbl = trackTbl[highConf]
			trackTbl.reset_index(inplace = True, drop = True)

	comp = comp[comp.cycle1.isin(trackTbl.cycle1)]
	comp.reset_index(inplace = True, drop = True)
	for latterCyc in params['registrationParams']['registrationPairs']:
		prevCyc = params['registrationParams']['registrationPairs'][latterCyc]
		comp['ngbrAlgnCost_' + prevCyc + '_vs_' + latterCyc] = \
			localNgbrAlgnCost(comp[[latterCyc, \
						'centroid-0' + '_' + latterCyc, \
						'centroid-1' + '_' + latterCyc]], \
					comp[[prevCyc, \
						'centroid-0' + '_' + prevCyc, \
						'centroid-1' + '_' + prevCyc]], \
					latterCyc, prevCyc, trackTbl, \
					np.min([ngbrsToConsider, trackTbl.shape[0]]), \
					approxmtInfCost)

	highConf = comp[[x for x in comp.columns if ('_vs_cycle' in x) and \
				('ngbrAlgnCost_' in x)]].apply( \
		lambda x: all([abs(v) < localNgbrAlgnCostCutoff for v in x]), axis = 1)
	if (highConf.sum() > 1):
		comp = comp[highConf]
		trackTbl = trackTbl[trackTbl.cycle1.isin(comp.cycle1)]
		trackTbl.reset_index(inplace = True, drop = True)
	trackTbl.to_csv(well + '.highConfValidMates.txt', \
                sep = '\t', index = None)
	nHighConf = trackTbl.shape[0]

if (smpl == 'LowConf'):
	lowConf = comp[[x for x in comp.columns if ('_vs_cycle' in x) and \
				('ngbrAlgnCost_' in x)]].apply( \
		lambda x: any([abs(v) > localNgbrAlgnCostCutoff - 0.1 for v in x]), axis = 1)
	comp = comp[lowConf]
	if (comp.shape[0] > 0):
		trackTbl = trackTbl[trackTbl.cycle1.isin(comp.cycle1)]
		trackTbl.reset_index(inplace = True, drop = True)
		nLowConf = lowConf.sum()
		if trackTbl.shape[0] == 0:
			print(well + " : 0 low confidence tracking mates.")
			Image.fromarray(np.zeros((2, 2, 3)).astype(np.uint8)).save('null.png')
			sys.exit()
	else:
		print(well + " : 0 low confidence tracking mates.")
		Image.fromarray(np.zeros((2, 2, 3)).astype(np.uint8)).save('null.png')
		sys.exit()

if trackTbl.shape[0] == 0:
	print(well + " : 0 high confidence tracking mates.")
	Image.fromarray(np.zeros((2, 2, 3)).astype(np.uint8)).save('null.png')
	sys.exit()
#-----------------------------------------------
#Draw a random sample
trackTbl = trackTbl.sample(n = min([nSamples, nLowConf, nHighConf]))
cols = trackTbl.columns
trackTbl = trackTbl.apply(lambda x: [well + 'n' + str(y) for y in x], \
		axis = 1, result_type = 'expand')
trackTbl.columns = cols
for cyc in cols:
	trackTbl[cyc].to_csv(well + '-' + cyc.replace('cycle', 'cyc') + '.' + smpl + '.txt', \
		header = None, index = None)

trackTbl.to_csv(well + '.' + smpl + '.allLabels.txt', header = None, index = None)
#------------------------------------------------
mosaicSuffix = '.mosaic.trimmed.npy'
chsOfInterest = []
allIds = pd.DataFrame()
stitchingRef = params['registrationParams']['refChannel']
for cycle in params['cycles']:
	chInfo = pd.read_csv(os.path.join(params['analysisDir'], \
			params['chToStainId']), sep = '\t')
	mainCh = chInfo[chInfo.cycle == 'cycle' + str(cycle)].T
	mainCh = mainCh[mainCh == stitchingRef].dropna().index[0][-1] + \
		'-cyc' + str(cycle)
	chsOfInterest.append('C' + mainCh + mosaicSuffix)
	allIds['cycle' + str(cycle)] = pd.read_csv( \
		well + '-cyc' + str(cycle) + '.' + smpl + '.txt', header = None)[0].tolist()

thisMaskExt = os.path.splitext( \
        [x for x in os.listdir(maskDir) if x.startswith(well) and \
                ('-cyc1' in x) and ('_cp_masks' in x)][0])[1]
ims = unpoolObjects( \
	fileWithListOfObjIds = well + '-cyc1.' + smpl + '.txt', \
	dirWithImages = mosaicDir, keepImIfPrefix = well, \
	keepImIfSuffix = mosaicSuffix, \
	objMaskFilePath = maskDir, \
	objMaskId = well + '-cyc1_' + \
		params['segmentationParams']['compartment'] + '_cp_masks' + thisMaskExt, \
	returnAllMasksInDir = True, returnIndObjMask = True, showOnlyTgtObjs = False, \
	dirWithMasks = maskDir, \
	keepMaskIfPrefix = well, dpWidth = dpWidth)

cyc1_ids = allIds['cycle1'].tolist()
for cyc in allIds.columns:
	for obj in allIds['cycle1'].tolist():
		thisObj = allIds[cyc].loc[allIds.cycle1 == obj].iloc[0]
		thisMaskExt = os.path.splitext( \
			[x for x in os.listdir(maskDir) if x.startswith(well) and \
				('-cyc' + cyc.replace('cycle', '') in x) and \
				('_cp_masks' in x)][0])[1]
		thisMaskId = well + '-cyc' + cyc.replace('cycle', '') + '_' + \
			params['segmentationParams']['compartment'] + '_cp_masks' + thisMaskExt
		ims['masks'][thisMaskId][obj] = \
			np.where(ims['masks'][thisMaskId][obj] == \
					int(thisObj[4:]), int(thisObj[4:]), 0)

normIms = OrderedDict({k1: OrderedDict({k2:np.where(v2 <= normalizers[ \
				k1.replace(mosaicSuffix, '')], v2, \
				normalizers[k1.replace(mosaicSuffix, '')])/normalizers[ \
				k1.replace(mosaicSuffix, '')] for \
	k2, v2 in v1.items()}) for \
	k1, v1 in ims['images'].items() if k1 in chsOfInterest})

nStains = len(chsOfInterest)
outmage_chs = 3
outmage = np.zeros((nSamples * dpWidth, nStains * dpWidth, outmage_chs))

for obj in range(np.min([nSamples, trackTbl.shape[0]])):
	for stn in range(nStains):
		thisOriShape = normIms[chsOfInterest[stn]][cyc1_ids[obj]].shape
		thisImage = normIms[chsOfInterest[stn]][cyc1_ids[obj]].copy()
		thisCycObjId = allIds['cycle' + chsOfInterest[stn].split('-cyc')[1].replace( \
						mosaicSuffix, '')].loc[ \
			allIds.cycle1 == cyc1_ids[obj]].iloc[0]
		thisMaskExt = os.path.splitext( \
			[x for x in os.listdir(maskDir) if x.startswith(well) and \
				('-' + chsOfInterest[stn].split('-')[1].replace(mosaicSuffix, '') \
				in x) and ('_cp_masks' in x)][0])[1]
		thisImage = mark_boundaries(thisImage, \
			ims['masks'][well + '-' + chsOfInterest[stn].split('-')[1].replace( \
				mosaicSuffix, '') + '_' + \
				params['segmentationParams']['compartment'] + \
				'_cp_masks' + thisMaskExt][cyc1_ids[obj]], \
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

