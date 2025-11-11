import os
import numpy as np
import skimage.io
from scipy import ndimage
import pandas as pd
from matplotlib import pyplot as plt
import sys
import json

thisWell = sys.argv[1]
thisCycle = sys.argv[2]
paramsFile = sys.argv[3]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

nTilesAlongDiameter = params['stitchingParams']['nTilesPerAxisInPalette']
bitDepth = params['bitDepth']
vigCorrDir = os.path.join(params['analysisDir'], params['vigCorrectDir'])
mosaicDir = os.path.join(params['analysisDir'], params['mosaicDir'])
coordsFile = os.path.join(mosaicDir, thisWell + '-cyc' + thisCycle + '.tileCentroidCoords.txt')
targetFiles = [x for x in os.listdir(vigCorrDir) if x.startswith(thisWell) \
			and ('cyc' + thisCycle in x)]
uniqueChannels = np.unique([x.split('-C')[1].replace('.vigCorr.tiff', '') for x in targetFiles])

tileCoords = pd.read_csv(coordsFile, sep = "\t")
tileCoords.fieldID = tileCoords.fieldID.astype(int)
#Add registration offsets between cycles
currCycle = thisCycle
rootForReg = 'F' + f"{params['registrationParams']['rootField']:02d}"
while 'cycle' + currCycle in params['registrationParams']['registrationPairs']:
	rgstrFile = os.path.join(mosaicDir, thisWell + '_' + rootForReg + '-cyc' + currCycle + \
		'_wrt_cyc' + params['registrationParams']['registrationPairs']['cycle' + \
						currCycle][-1] + '.coords.txt')
	rgstrCoords = pd.read_csv(rgstrFile, sep = '\t')
	tileCoords.medX += round(rgstrCoords.X.iloc[0], 1)
	tileCoords.medY += round(rgstrCoords.Y.iloc[0], 1)
	currCycle = params['registrationParams']['registrationPairs']['cycle' + currCycle][-1]
#-------------------------------------
#Process channel by channel
for thisCh in uniqueChannels:
	outMosaic = os.path.join(thisWell + '-C' + thisCh + '.mosaic.npy')
	outMosaicPng = os.path.join(thisWell + '-C' + thisCh + '.mosaic.png') 
	#--------------------------------------
	#Load images
	#--------------------------------------
	images = {}
	thisFiles = [x for x in targetFiles if '-C' + thisCh in x]
	for i in thisFiles:
		images[i[4:7]] = skimage.io.imread(os.path.join(vigCorrDir, i))

	#--------------------------------------
	imageSize = images['F01'].shape[0]
	mosaicSize = images['F01'].shape[0] * nTilesAlongDiameter
	
	mosaic = np.zeros([mosaicSize, mosaicSize])
	centerStart = int((mosaicSize - imageSize)/2)
	centerEnd = int((mosaicSize + imageSize)/2)
	
	for thisField in tileCoords.fieldID:
		if (thisField == params['registrationParams']['rootField']):
			mosaic[centerStart:centerEnd, centerStart:centerEnd] = \
				images['F' + f'{thisField:02d}'].copy()
			continue
		if (np.isnan(tileCoords[tileCoords['fieldID'] == thisField].medY.iloc[0])):
			continue
		shiftX = tileCoords[tileCoords['fieldID'] == thisField].medX.iloc[0]
		shiftY = tileCoords[tileCoords['fieldID'] == thisField].medY.iloc[0]
		
		#Layer for 1st order interpolation
		intrpltLayer = np.zeros([imageSize + 2, imageSize + 2])
		intrpltLayer[1:(imageSize + 1), 1:(imageSize + 1)] = \
			images['F' + f'{thisField:02d}'].copy()
		intrpltLayer = ndimage.shift(intrpltLayer, \
				[round(shiftY - int(shiftY), 1), round(shiftX - int(shiftX), 1)], \
				order = 1)
		intrpltLayer = np.clip(intrpltLayer, 0, 2**bitDepth - 1)
		startY = centerStart + int(shiftY) - 1
		endY = centerStart + int(shiftY) + imageSize + 1
		startX = centerStart + int(shiftX) - 1
		endX = centerStart + int(shiftX) + imageSize + 1
		mosaic[startY:endY, startX:endX] = np.max(np.dstack([intrpltLayer, 
			mosaic[startY:endY, startX:endX]]), axis = 2)
	
	#--------------------------------
	with open(outMosaic, 'wb') as f:
		np.save(f, mosaic)

	#--------------------------------
	plt.figure(figsize=(100, 100))
	plt.imshow(mosaic, cmap = 'gray')
	plt.savefig(outMosaicPng, bbox_inches='tight')
	plt.close()





