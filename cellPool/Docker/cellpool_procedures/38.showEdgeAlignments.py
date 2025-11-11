import os
import numpy as np
import pandas as pd
import sys
from skimage.io import imsave
import json
from skimage.exposure import rescale_intensity
from PIL import Image, ImageDraw, ImageFont

wellToProcess = sys.argv[1]
paramsFile = sys.argv[2]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.utils import *
from cellpool.mosaic import *

stitchingRef = params['registrationParams']['refChannel']
tileBorderZoneHalfWidth = params['showStitchingParams']['tileBorderZoneHalfWidth']
fieldsPerWell = params['showStitchingParams']['fieldsPerWell']
sepWidth = params['showStitchingParams']['sepWidth'] #Relative to tileBorderZoneHalfWidth
fontsize = params['showStitchingParams']['fontSize']
font = '/app/UbuntuMono-R.ttf'

inDir = os.path.join(params['analysisDir'], params['mosaicDir'])

rawMosaicShapes = pd.read_csv(wellToProcess + '.rawMosaicShapes.txt', \
		sep = '\t')

chInfo = pd.read_csv(os.path.join(params['analysisDir'], params['chToStainId']), sep = '\t')
#-------------------------------------
for cycle in params['cycles']:
	nucCh = chInfo[chInfo.cycle == 'cycle' + str(cycle)].T
	nucCh = nucCh[nucCh == stitchingRef].dropna().index[0][-1] + \
		'-cyc' + str(cycle)

	fieldsData = pd.read_csv(os.path.join(params['analysisDir'], \
				params['wellMapFile']), sep = '\t')
	fieldsData = fieldsData.loc[fieldsData.cycle == cycle]
	stitchCoords = pd.read_csv(os.path.join(inDir, wellToProcess + \
                '-cyc' + str(cycle) + '.tileCentroidCoords.txt'), sep = '\t')

	#Add registration offsets between cycles
	currCycle = str(cycle)
	rootForReg = 'F' + f"{params['registrationParams']['rootField']:02d}"
	while 'cycle' + currCycle in params['registrationParams']['registrationPairs']:
		rgstrFile = os.path.join(inDir, wellToProcess + '_' + \
			rootForReg + '-cyc' + currCycle + '_wrt_cyc' + \
			params['registrationParams']['registrationPairs']['cycle' + \
							currCycle][-1] + '.coords.txt')
		rgstrCoords = pd.read_csv(rgstrFile, sep = '\t')
		stitchCoords.medX += round(rgstrCoords.X.iloc[0], 1)
		stitchCoords.medY += round(rgstrCoords.Y.iloc[0], 1)
		currCycle = params['registrationParams']['registrationPairs']['cycle' + \
										currCycle][-1]

	stitchCoords['imX1'] = int(rawMosaicShapes['cycle' + str(cycle)].iloc[0]/2) + \
					stitchCoords['medX'].astype(int) - \
		int(fieldsData.ImageSizeX.iloc[0]/2)
	stitchCoords['imX2'] = int(rawMosaicShapes['cycle' + str(cycle)].iloc[0]/2) + \
					stitchCoords['medX'].astype(int) + \
		int(fieldsData.ImageSizeX.iloc[0]/2)
	stitchCoords['imY1'] = int(rawMosaicShapes['cycle' + str(cycle)].iloc[0]/2) + \
					stitchCoords['medY'].astype(int) - \
		int(fieldsData.ImageSizeY.iloc[0]/2)
	stitchCoords['imY2'] = int(rawMosaicShapes['cycle' + str(cycle)].iloc[0]/2) + \
					stitchCoords['medY'].astype(int) + \
		int(fieldsData.ImageSizeY.iloc[0]/2)
	stitchCoords['stdR'] = stitchCoords[['stdX', 'stdY']].apply(lambda x: \
		(x[0]**2 + x[1]**2)**0.5, axis = 1)	

	fldsToShow = stitchCoords.nlargest(fieldsPerWell, 'stdR').fieldID.tolist()
	fldsToShow = [int(x) for x in fldsToShow]
	im = np.load(os.path.join(inDir, wellToProcess + '-C' + nucCh + \
			'.mosaic.npy'))
	sepColor = []
	for fld in fldsToShow:
		thisX1 = stitchCoords.imX1.iloc[fld - 1]
		thisX2 = stitchCoords.imX2.iloc[fld - 1]
		thisY1 = stitchCoords.imY1.iloc[fld - 1]
		thisY2 = stitchCoords.imY2.iloc[fld - 1]
		thisStitchIm = np.full((int(tileBorderZoneHalfWidth*(2 + \
			sepWidth))*5, thisY2 - thisY1 + \
			2*tileBorderZoneHalfWidth), fill_value = -1)

		thisStitchIm[:2*tileBorderZoneHalfWidth, :] = \
			im[(thisY1 - tileBorderZoneHalfWidth):(thisY1 + tileBorderZoneHalfWidth), \
			(thisX1 - tileBorderZoneHalfWidth):(thisX2 + tileBorderZoneHalfWidth)]
		thisStitchIm[int(tileBorderZoneHalfWidth*(2 + sepWidth)):( \
			int(tileBorderZoneHalfWidth*(2 + sepWidth)) + \
			2*tileBorderZoneHalfWidth), :] = \
			np.transpose(im[(thisY1 - tileBorderZoneHalfWidth):(thisY2 + \
				tileBorderZoneHalfWidth), \
			(thisX2 - tileBorderZoneHalfWidth):(thisX2 + tileBorderZoneHalfWidth)])
		thisStitchIm[int(tileBorderZoneHalfWidth*(2 + sepWidth)*2):( \
			int(tileBorderZoneHalfWidth*(2 + sepWidth)*2) + \
			2*tileBorderZoneHalfWidth), :] = \
			im[(thisY2 - tileBorderZoneHalfWidth):(thisY2 + tileBorderZoneHalfWidth), \
			(thisX1 - tileBorderZoneHalfWidth):(thisX2 + tileBorderZoneHalfWidth)]
		thisStitchIm[int(tileBorderZoneHalfWidth*(2 + sepWidth)*3):( \
			int(tileBorderZoneHalfWidth*(2 + sepWidth)*3) + \
			2*tileBorderZoneHalfWidth), :] = \
			np.transpose(im[(thisY1 - tileBorderZoneHalfWidth):(thisY2 + \
				tileBorderZoneHalfWidth), \
			(thisX1 - tileBorderZoneHalfWidth):(thisX1 + tileBorderZoneHalfWidth)])

		sepColor.append(np.quantile(thisStitchIm, 0.95))
		if (fld == fldsToShow[0]):
			thisCycleIm = thisStitchIm
		else:
			sepFields = np.full((int(tileBorderZoneHalfWidth*(2 + \
				sepWidth))*5, int(2*sepWidth*tileBorderZoneHalfWidth)), \
				fill_value = -1)
			thisCycleIm = np.concatenate((thisCycleIm, sepFields, thisStitchIm), axis=1)

	thisCycleIm = np.where(thisCycleIm == -1, np.max(sepColor), thisCycleIm)	

	#Annotate image with well and field labels
	thisCycleIm = Image.fromarray(rescale_intensity(np.clip(thisCycleIm, \
		a_min = 0, a_max = np.max(sepColor)), \
		out_range = (0, 255)).astype(np.uint8))
	write = ImageDraw.Draw(thisCycleIm)
	myFont = ImageFont.truetype(font, size = fontsize)
	write.text((tileBorderZoneHalfWidth, thisCycleIm.size[1] - tileBorderZoneHalfWidth), \
		wellToProcess + ' : showing ' + ', '.join([str(x) for x in fldsToShow]), \
		font=myFont, fill = (0,))

	wellOrder = pd.read_csv(os.path.join(params['analysisDir'], params['listOfWells']), \
			header = None)[0].tolist()
	thisWellOrder = [x + 1 for x in range(len(wellOrder)) if wellOrder[x] == wellToProcess][0]
	thisCycleIm.save(f"{thisWellOrder:02d}" + '.cycle' + str(cycle) + '.showOverlap.png')


