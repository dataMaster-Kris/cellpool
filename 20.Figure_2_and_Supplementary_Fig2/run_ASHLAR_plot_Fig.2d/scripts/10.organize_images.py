import os
import sys
import json
import pandas as pd
import numpy as np
import itertools
from skimage.io import imsave

well = sys.argv[1]
paramsFile = sys.argv[2]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.mosaic import *

vigCorrDir = os.path.join(params['analysisDir'], params['vigCorrectDir'])
wells = pd.read_csv(os.path.join(params['analysisDir'], params['listOfWells']), \
	header = None).values.flatten()
wellFiles = os.listdir(vigCorrDir)

fieldsData = pd.read_csv(os.path.join(params['analysisDir'], params['wellMapFile']), sep = '\t')
tileGph = buildNeighborGraphOfTiles(fieldsData, \
		maxSeparation = params['mosaickingParams']['maxTileSepToTestOverlap'])
fieldsData = addNghbrInfo(fieldsData, tileGph)

fieldsData = fieldsData[['id', 'cycle', 'FieldID', 'ImageSizeX', 'ImageSizeY', \
	'BlockDistanceFromField1X', 'BlockDistanceFromField1Y']]
fieldsData['BlockDistanceFromField1X'] = np.sign(fieldsData['BlockDistanceFromField1X']) * \
	np.ceil(np.abs(fieldsData['BlockDistanceFromField1X']))
fieldsData['BlockDistanceFromField1Y'] = np.sign(fieldsData['BlockDistanceFromField1Y']) * \
	np.ceil(np.abs(fieldsData['BlockDistanceFromField1Y']))

fieldsData['BlockDistanceFromField1X'] = np.max(fieldsData['BlockDistanceFromField1X']) + \
	fieldsData['BlockDistanceFromField1X']
fieldsData['BlockDistanceFromField1Y'] = np.max(fieldsData['BlockDistanceFromField1Y']) - \
	fieldsData['BlockDistanceFromField1Y']

inChToElemId = pd.read_csv(os.path.join(params['analysisDir'], params['chToElemIdFile']), \
	sep = '\t')
outChToElemId = pd.read_csv(os.path.join(params['analysisDir'], params['ashlarAnalysisDir'], \
	params['reorderedChToElemId']), sep = '\t')
inOutPair = {}
for cyc in params['cycles']:
	inOutPair['cyc' + str(cyc)] = {}
	outElem = outChToElemId.loc[(outChToElemId.cycle == 'cycle' + \
		str(cyc))][['Channel.' + str(x) for x in range(1, 5)]].values[0]
	for ch in range(1, 5):
		inElem = inChToElemId.loc[inChToElemId.cycle == 'cycle' + str(cyc), \
			'Channel.' + str(ch)].iloc[0]
		outCh = [x for x in range(1, 5) if outElem[x - 1] == inElem]
		if len(outCh) == 1:
			outCh = outCh[0]
		else:
			continue
		inOutPair['cyc' + str(cyc)][str(ch)] = str(outCh)

thisWellFiles = [x for x in wellFiles if x.startswith(well)]
for cyc in params['cycles']:
	thisFldsData = fieldsData.loc[fieldsData.cycle == cyc]
	thisCycFiles = [x for x in thisWellFiles if 'cyc' + str(cyc) in x]
	imSizeX = int(np.unique(thisFldsData.ImageSizeX)[0])
	imSizeY = int(np.unique(thisFldsData.ImageSizeY)[0])
	fileDat = [x.split('_F')[1].split('-C') for x in thisCycFiles]
	channels = np.unique([x[1] for x in fileDat])
	for ch in channels:
		thisCycChFiles = [x for x in thisCycFiles if x.endswith(ch)]
		gridCoords = [x for x in itertools.product( \
			thisFldsData.BlockDistanceFromField1Y.unique(), \
			thisFldsData.BlockDistanceFromField1X.unique())]
		fileDat = [x.split('_F')[1].split('-C') for x in thisCycChFiles]
		for ix in range(len(thisCycChFiles)):
			inName = thisCycChFiles[ix]
			fldImg = fileDat[ix]
			row, col = fieldsData.loc[(fieldsData.FieldID == int(fldImg[0])) & \
				(fieldsData.cycle == int(cyc))][[ \
				'BlockDistanceFromField1Y', \
				'BlockDistanceFromField1X']].iloc[0].tolist()
			outName = well + '_R' + f"{int(row):02d}" + 'C' + \
				f"{int(col):02d}" + '-C' + inOutPair['cyc' + \
					str(cyc)][fldImg[1][0]] + fldImg[1][1:]
			gridCoords = [x for x in gridCoords if (int(row), int(col)) != x]
			os.symlink(os.path.join(vigCorrDir, inName), \
				os.path.join(params['analysisDir'], params['ashlarAnalysisDir'], \
					'10.vigCorrect', outName))

		blankIm = np.zeros((imSizeY, imSizeX))
		for ix in range(len(gridCoords)):
			row, col = gridCoords[ix]
			outName = well + '_R' + f"{int(row):02d}" + 'C' + \
				f"{int(col):02d}" + '-C' + inOutPair['cyc' + \
					str(cyc)][ch[0]] + ch[1:]
			imsave(outName, blankIm, \
				check_contrast = False)

