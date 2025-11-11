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

inDir = os.path.join(params['analysisDir'], params['mosaicDir'])

#Load one image and all the tile coords.
#Get the borders
im = np.load(os.path.join(inDir, wellToProcess + '-C1-cyc1.mosaic.npy'))

border = {'vertical' : pd.DataFrame(), 'horizontal' : pd.DataFrame()}
for cyc in params['cycles']:
	fieldsData = pd.read_csv(os.path.join(params['analysisDir'], \
				params['wellMapFile']), sep = '\t')
	fieldsData = fieldsData.loc[fieldsData.cycle == cyc]
	tileGph = buildNeighborGraphOfTiles(fieldsData, \
			maxSeparation = params['mosaickingParams']['maxTileSepToTestOverlap'])
	fieldsData = addNghbrInfo(fieldsData, tileGph)
	tileCoords = pd.read_csv(os.path.join(inDir, wellToProcess + \
		'-cyc' + str(cyc) + '.tileCentroidCoords.txt'), sep = '\t')

	#Add registration offsets between cycles
	currCycle = str(cyc)
	rootForReg = 'F' + f"{params['registrationParams']['rootField']:02d}"
	while 'cycle' + currCycle in params['registrationParams']['registrationPairs']:
		rgstrFile = os.path.join(inDir, wellToProcess + '_' + \
			rootForReg + '-cyc' + currCycle + '_wrt_cyc' + \
			params['registrationParams']['registrationPairs']['cycle' + \
							currCycle][-1] + '.coords.txt')
		rgstrCoords = pd.read_csv(rgstrFile, sep = '\t')
		tileCoords.medX += round(rgstrCoords.X.iloc[0], 1)
		tileCoords.medY += round(rgstrCoords.Y.iloc[0], 1)
		currCycle = params['registrationParams']['registrationPairs']['cycle' + \
										currCycle][-1]

	borderVerticals = []
	borderHorizontals = []
	#Find borders with north side external
	for blockDistY in range(int(np.floor(np.min(fieldsData.BlockDistanceFromField1Y))), 0):
		edgeTiles = fieldsData.FieldID[(abs(fieldsData.BlockDistanceFromField1Y + \
			blockDistY * (1 - overlap)) < 0.001) & \
			np.isnan(fieldsData.ngbrToNorth)].tolist()
		if (len(edgeTiles) == 0):
			continue
		startY = []
		for thisField in edgeTiles:
			shiftY = tileCoords[tileCoords['fieldID'] == thisField].medY.iloc[0]
			centerStart = int((im.shape[0] - fieldsData.ImageSizeY.iloc[0])/2)
			startY.append(centerStart + int(shiftY) - 1)
		borderY = max(startY)
		borderHorizontals.append(borderY)

	#Find borders with south side external
	for blockDistY in range(int(np.ceil(np.max(fieldsData.BlockDistanceFromField1Y))), 0, -1):
		edgeTiles = fieldsData.FieldID[(abs(fieldsData.BlockDistanceFromField1Y + \
			blockDistY * (1 - overlap)) < 0.001) & \
			np.isnan(fieldsData.ngbrToSouth)].tolist()
		if (len(edgeTiles) == 0):
			continue
		endY = []
		for thisField in edgeTiles:
			shiftY = tileCoords[tileCoords['fieldID'] == thisField].medY.iloc[0]
			centerStart = int((im.shape[0] - fieldsData.ImageSizeY.iloc[0])/2)
			endY.append(centerStart + int(shiftY) + fieldsData.ImageSizeY.iloc[0] + 1)
		borderY = min(endY)
		borderHorizontals.append(borderY)

	#Find borders with west side external
	for blockDistX in range(int(np.ceil(np.max(fieldsData.BlockDistanceFromField1X))), 0, -1):
		edgeTiles = fieldsData.FieldID[(abs(fieldsData.BlockDistanceFromField1X + \
			blockDistX * (1 - overlap)) < 0.001) & \
			np.isnan(fieldsData.ngbrToWest)].tolist()
		if (len(edgeTiles) == 0):
			continue
		startX = []
		for thisField in edgeTiles:
			shiftX = tileCoords[tileCoords['fieldID'] == thisField].medX.iloc[0]
			centerStart = int((im.shape[0] - fieldsData.ImageSizeX.iloc[0])/2)
			startX.append(centerStart + int(shiftX) - 1)
		borderX = max(startX)
		borderVerticals.append(borderX)

	#Find borders with east side external
	for blockDistX in range(int(np.floor(np.min(fieldsData.BlockDistanceFromField1X))), 0):
		edgeTiles = fieldsData.FieldID[(abs(fieldsData.BlockDistanceFromField1X + \
			blockDistX * (1 - overlap)) < 0.001) & \
			np.isnan(fieldsData.ngbrToEast)].tolist()
		if (len(edgeTiles) == 0):
			continue
		endX = []
		for thisField in edgeTiles:
			shiftX = tileCoords[tileCoords['fieldID'] == thisField].medX.iloc[0]
			centerStart = int((im.shape[0] - fieldsData.ImageSizeX.iloc[0])/2)
			endX.append(centerStart + int(shiftX) + fieldsData.ImageSizeX.iloc[0] + 1)
		borderX = min(endX)
		borderVerticals.append(borderX)

	borderVerticals.sort()
	borderHorizontals.sort()
	border['vertical']['cycle' + str(cyc)] = borderVerticals
	border['horizontal']['cycle' + str(cyc)] = borderHorizontals

border['vertical'] = border['vertical'].apply( \
	lambda x: [min(x), max(x)][all([val < im.shape[0]/2 for val in x]) + 0], \
	axis = 1)
border['horizontal'] = border['horizontal'].apply( \
	lambda x: [min(x), max(x)][all([val < im.shape[0]/2 for val in x]) + 0], \
	axis = 1)

#Save for later use when marking cells in tile overlsp
pd.DataFrame.from_dict(border).to_csv(wellToProcess + \
	'.trimBordersWrtRawMosaic.txt', sep = '\t', index = None)
#----------------------------------------------
#Find pixels to make nan
whichToNan = np.zeros(im.shape)
borderHorizontals = border['horizontal'].tolist()
borderVerticals = border['vertical'].tolist()
whichToNan = whichToNan[borderHorizontals[0]:(borderHorizontals[-1] + 1), \
	borderVerticals[0]:(borderVerticals[-1] + 1)].copy()
borderVerticals = [x - min(borderVerticals) for x in borderVerticals]
borderHorizontals = [x - min(borderHorizontals) for x in borderHorizontals]

#Trim the left side
whichToNan[:borderHorizontals[1], :borderVerticals[3]] = np.nan
whichToNan[:borderHorizontals[2], :borderVerticals[2]] = np.nan
whichToNan[:borderHorizontals[3], :borderVerticals[1]] = np.nan
whichToNan[borderHorizontals[4]:borderHorizontals[5], :borderVerticals[1]] = np.nan
whichToNan[borderHorizontals[5]:borderHorizontals[6], :borderVerticals[2]] = np.nan
whichToNan[borderHorizontals[6]:borderHorizontals[7], :borderVerticals[3]] = np.nan

#Trim the right side
whichToNan[:borderHorizontals[1], borderVerticals[4]:] = np.nan
whichToNan[:borderHorizontals[2], borderVerticals[5]:] = np.nan
whichToNan[:borderHorizontals[3], borderVerticals[6]:] = np.nan
whichToNan[borderHorizontals[4]:borderHorizontals[5], borderVerticals[6]:] = np.nan
whichToNan[borderHorizontals[5]:borderHorizontals[6], borderVerticals[5]:] = np.nan
whichToNan[borderHorizontals[6]:borderHorizontals[7], borderVerticals[4]:] = np.nan

center = [round(whichToNan.shape[0]/2), round(whichToNan.shape[1]/2)]
radius = min(center) - extraTrim
whichToNan = whichToNan[(center[0] - radius):(center[0] + radius), \
		(center[1] - radius):(center[1] + radius)]
borderVerticals = [x - center[0] + radius for x in borderVerticals]
borderHorizontals = [x - center[1] + radius for x in borderHorizontals]
for i in range(whichToNan.shape[0]):
	for j in range(whichToNan.shape[1]):
		if ((i-center[0])**2 + (j-center[0])**2 >= radius**2):
			whichToNan[i, j] = np.nan

#-------------------------------------
#Trim and save
outDat = {}
qntls = np.linspace(0, 1, num = 101)
rawImShapes = {}
for cyc in params['cycles']:
	thisCycFiles = [x for x in os.listdir(inDir) if \
		x.endswith('.mosaic.npy') and x.startswith(wellToProcess) and \
		('cyc' + str(cyc) in x)]
	thisCycle = str(cyc)
	outDat[thisCycle] = {}
	for nextFile in thisCycFiles:
		im = np.load(os.path.join(inDir, nextFile))
		rawImShapes['cycle' + str(cyc)] = im.shape
		im = im[border['horizontal'].iloc[0]:(border['horizontal'].iloc[-1] + 1), \
			border['vertical'].iloc[0]:(border['vertical'].iloc[-1] + 1)]
		im = im[(center[0] - radius):(center[0] + radius), \
			(center[1] - radius):(center[1] + radius)]
		im[np.isnan(whichToNan)] = np.nan
		outFile = nextFile.replace('.mosaic.npy', '.mosaic.trimmed.npy')
		with open(outFile, 'wb') as f:
			np.save(f, im)
		thisQntls = np.nanquantile(im, q = qntls)
		thisCh = nextFile.split('-C')[1].split('-cyc')[0]
		outDat[thisCycle][thisCh] = \
			dict(zip(['q_' + str(round(x, 2)) \
						for x in qntls], thisQntls))
		outDat[thisCycle][thisCh]['ch'] = thisCh
		outDat[thisCycle][thisCh]['cycle'] = thisCycle

outDat = pd.concat([pd.DataFrame.from_dict(x).T for x in outDat.values()])
outDat.to_csv(wellToProcess + '.rawIntensityQuantiles.txt', \
	sep = '\t', index = False)

outBorderFile = wellToProcess + '.trimBordersWrtTrimmedMosaic.txt'
pd.DataFrame.from_dict({'vertical': borderVerticals, \
	'horizontal': borderHorizontals}).to_csv(outBorderFile, sep = '\t', index = None)

outRawShapesFile = wellToProcess + '.rawMosaicShapes.txt'
pd.DataFrame.from_dict(rawImShapes).to_csv(outRawShapesFile, \
	sep = '\t', index = None)

