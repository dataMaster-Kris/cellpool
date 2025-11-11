import os
import numpy as np
import pandas as pd
import sys
import skimage.io
import json
from skimage.exposure import rescale_intensity

wellToProcess = sys.argv[1]
paramsFile = sys.argv[2]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.utils import *
from cellpool.mosaic import *

stitchingRef = params['registrationParams']['refChannel']
overlap = params['tileOverlap']
extraTrim = params['stitchingParams']['extraTrim']
tileBorderZoneHalfWidth = params['stitchingParams']['tileBorderZoneHalfWidth']

inDir = os.path.join(params['analysisDir'], params['mosaicDir'])

border = pd.read_csv(wellToProcess + '.trimBordersWrtRawMosaic.txt', \
		sep = '\t')
rawMosaicShapes = pd.read_csv(wellToProcess + '.rawMosaicShapes.txt', \
		sep = '\t')
#----------------------------------------------
borderHorizontals = border['horizontal'].tolist()
borderVerticals = border['vertical'].tolist()
center = [round((borderHorizontals[-1] - borderHorizontals[0] + 1)/2), \
		round((borderVerticals[-1] - borderVerticals[0] + 1)/2)]
radius = min(center) - extraTrim

#-------------------------------------
#Initialize an image with size the same as the raw mosaics, highlight overlap regions, ...
#..., trim and save 
for cycle in params['cycles']:
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

	borderIm = np.zeros(tuple(rawMosaicShapes['cycle' + str(cycle)]), dtype = np.int8)
	for fld in range(stitchCoords.shape[0]):
		thisX1 = stitchCoords.imX1.iloc[fld]
		thisX2 = stitchCoords.imX2.iloc[fld]
		thisY1 = stitchCoords.imY1.iloc[fld]
		thisY2 = stitchCoords.imY2.iloc[fld]
		borderIm[(thisY1 - tileBorderZoneHalfWidth):(thisY2 + tileBorderZoneHalfWidth), \
			(thisX1 - tileBorderZoneHalfWidth):(thisX1 + tileBorderZoneHalfWidth)] = 1
		borderIm[(thisY1 - tileBorderZoneHalfWidth):(thisY2 + tileBorderZoneHalfWidth), \
			(thisX2 - tileBorderZoneHalfWidth):(thisX2 + tileBorderZoneHalfWidth)] = 1
		borderIm[(thisY1 - tileBorderZoneHalfWidth):(thisY1 + tileBorderZoneHalfWidth), \
			(thisX1 - tileBorderZoneHalfWidth):(thisX2 + tileBorderZoneHalfWidth)] = 1
		borderIm[(thisY2 - tileBorderZoneHalfWidth):(thisY2 + tileBorderZoneHalfWidth), \
			(thisX1 - tileBorderZoneHalfWidth):(thisX2 + tileBorderZoneHalfWidth)] = 1

	borderIm = borderIm[border['horizontal'].iloc[0]:(border['horizontal'].iloc[-1] + 1), \
		border['vertical'].iloc[0]:(border['vertical'].iloc[-1] + 1)].copy()
	borderIm = borderIm[(center[0] - radius):(center[0] + radius), \
		(center[1] - radius):(center[1] + radius)]
	borderIm = borderIm.astype(float)
	outFile = os.path.join(wellToProcess + '-cyc' + str(cycle) + '.trimmed.tileOverlapZone.png') 
	skimage.io.imsave(outFile, borderIm.astype(np.uint8))


