import os
import numpy as np
import sys
import json
import pandas as pd
from collections import OrderedDict
from skimage.segmentation import mark_boundaries
from PIL import Image, ImageDraw, ImageFont
from tifffile import imread

wellOmeTiff = sys.argv[1]
paramsFile = sys.argv[2]
with open(paramsFile, 'r') as inFile:
        params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.unpool import *

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

chQtlDir = os.path.join(params['analysisDir'], params['mosaicDir'])
maskDir = os.path.join(params['analysisDir'], params['ashlarAnalysisDir'], \
	params['segmentationDir'], \
	'iter' + str(iteration) + params['segmentationParams']['maskDirSuffix'])
normalizers = pd.read_csv(os.path.join(chQtlDir, \
	"wholePlate.rawIntensityQntls.txt"), sep = "\t")
normalizers = dict(zip(normalizers[['cycle', 'ch']].apply( \
		lambda x: 'C' + str(x[1]) + '-cyc' + str(x[0]), axis = 1).tolist(), \
		normalizers[pctl4Normalization].tolist()))
wellOmeTiffName = os.path.basename(wellOmeTiff)
well = wellOmeTiffName[:3]
nucChs = [0, 2, 6, 10]
im = imread(wellOmeTiff)*65535
im = np.transpose(im, axes = [1, 2, 0])
mask = imageio.imread(os.path.join(maskDir, wellOmeTiffName.replace('.ome.tiff', '') + \
	'-cyc1_' + params['segmentationParams']['compartment'] + '_cp_masks.png'))

np.random.seed(90283409) #Randomly typed number
toKeep = np.random.choice([x for x in range(0, mask.max())], nSamples, replace = False)
mask[~np.isin(mask, toKeep)] = 0
dps = OrderedDict()

for ch in nucChs:
	dps['ch' + str(ch)] = expandBboxAndExtract(mask, im[:, :, ch], dpWidth)
	dps['ch' + str(ch)]['images'] = OrderedDict()
	for lbl in range(nSamples):
		dps['ch' + str(ch)]['images']['n' + \
			str(dps['ch' + str(ch)]['label'].iloc[lbl])] = \
				dps['ch' + str(ch)]['cellSnaps'][lbl]
	dps['ch' + str(ch)] = dps['ch' + str(ch)]['images'] 

indMasks = expandBboxAndExtract(mask, mask, dpWidth)

normalizers = dict(zip(['ch' + str(x) for x in nucChs], \
	[normalizers[x] for x in ['C1-cyc1', 'C2-cyc2', 'C1-cyc3', 'C1-cyc4']]))
normIms = OrderedDict({k1: OrderedDict({k2:np.where(v2 <= normalizers[k1], v2, \
				normalizers[k1])/normalizers[k1] for \
	k2, v2 in v1.items()}) for \
	k1, v1 in dps.items()})

nStains = len(nucChs)
outmage_chs = 3
outmage = np.zeros((nSamples * dpWidth, nStains * dpWidth, outmage_chs))

cyc1_ids = list(normIms['ch0'].keys())
for obj in range(nSamples):
	for stn in range(nStains):
		thisOriShape = normIms['ch' + str(nucChs[stn])][cyc1_ids[obj]].shape
		thisImage = normIms['ch' + str(nucChs[stn])][cyc1_ids[obj]].copy()
		thisImage = mark_boundaries(thisImage, \
			indMasks['cellSnaps'][obj], \
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
	str(nSamples) + ' random samples of mask transfer on ASHLAR outputs in ' + \
	well, font=myFont, fill = (0,))

outmage2.save(os.path.join(f"{wellOrder:02d}" + '.png'))

