import os
import numpy as np
import pandas as pd
import skimage.io
from scipy import ndimage
import json
import sys

wellToProcess = sys.argv[1]
targetCycle = sys.argv[2]
paramsFile = sys.argv[3]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.utils import *
from cellpool.mosaic import *

maxShift = params['registrationParams']['maxShift']
rootForReg = params['registrationParams']['rootField']
refCycleForRegistration = params['registrationParams']['registrationPairs']['cycle' + targetCycle]
refTileCoordsFile = os.path.join(params['analysisDir'], params['mosaicDir'], wellToProcess + \
	'-cyc' + refCycleForRegistration[-1] + '.tileCentroidCoords.txt') 
tgtTileCoordsFile = os.path.join(params['analysisDir'], params['mosaicDir'], wellToProcess + \
	'-cyc' + targetCycle + '.tileCentroidCoords.txt')
vigCorrDir = os.path.join(params['analysisDir'], params['vigCorrectDir'])

outFileAllCoords = wellToProcess + '_F' + f'{rootForReg:02d}' + '-cyc' + \
		targetCycle + '_coords_wrt_cyc' + str(refCycleForRegistration)[-1] + '.' + \
		params['registrationParams']['refChannel'].replace(" ", "_") + '.intmdt.txt'
outFileSmryCoords = wellToProcess + '_F' + f'{rootForReg:02d}' + '-cyc' + \
		targetCycle + '_wrt_cyc' + str(refCycleForRegistration)[-1] + '.coords.txt'

#-----------------------------------------
#Load field metadata for cycles
#-----------------------------------------
fieldsData = pd.read_csv(os.path.join(params['analysisDir'], \
			params['wellMapFile']), sep = '\t')
tgtFieldsData = fieldsData.loc[fieldsData.cycle == int(targetCycle)]
refFieldsData = fieldsData.loc[fieldsData.cycle == int(str(refCycleForRegistration)[-1])]

#--------------------------------------
#Load images
#--------------------------------------
allFiles = os.listdir(vigCorrDir)
allFiles = [x for x in allFiles if \
	x.startswith(wellToProcess) and \
	x.endswith('.vigCorr.tiff')]

chInfo = pd.read_csv(os.path.join(params['analysisDir'], params['chToStainId']), sep = '\t')

refUniqCh = list(set([x[9:] for x in allFiles if ('cyc' + str(refCycleForRegistration)[-1] in x)]))
refUniqCh.sort()
refChInfo = chInfo.loc[chInfo.cycle == 'cycle' + str(refCycleForRegistration)[-1]]
refChInfo = pd.DataFrame.from_dict({x for x in zip([x.strip('.vigCorr.tiff') for x in refUniqCh], \
		refChInfo[['Channel.' + \
		str(x[0]) for x in refUniqCh]].iloc[0].tolist())})
refChInfo.columns = ['channel', 'channelName']
refChIdInRefCyc = refChInfo[refChInfo.channelName == \
			params['registrationParams']['refChannel']].channel.iloc[0]

tgtUniqCh = list(set([x[9:] for x in allFiles if ('cyc' + targetCycle in x)]))
tgtUniqCh.sort()
tgtChInfo = chInfo.loc[chInfo.cycle == 'cycle' + targetCycle]
tgtChInfo = pd.DataFrame.from_dict({x for x in zip([x.strip('.vigCorr.tiff') for x in tgtUniqCh], \
		tgtChInfo[['Channel.' + \
		str(x[0]) for x in tgtUniqCh]].iloc[0].tolist())})
tgtChInfo.columns = ['channel', 'channelName']
refChIdInTgtCyc = tgtChInfo[tgtChInfo.channelName == \
			params['registrationParams']['refChannel']].channel.iloc[0]

targetFiles = [x for x in allFiles if ('-C' + refChIdInTgtCyc + '.vigCorr.tiff' in x)]
tgtImages = {}
for i in targetFiles:
	tgtImages[i.split('-')[0][4:7]] = \
		skimage.io.imread(os.path.join(vigCorrDir, i))

refFiles = [x for x in allFiles if ('-C' + refChIdInRefCyc + '.vigCorr.tiff' in x)]
refImages = {}
for i in refFiles:
	refImages[i.split('-')[0][4:7]] = \
		skimage.io.imread(os.path.join(vigCorrDir, i))
#-------------------------------------
#-------------------------------------
refTileCoords = pd.read_csv(refTileCoordsFile, sep = '\t') 
tgtTileCoords = pd.read_csv(tgtTileCoordsFile, sep = '\t')
layeredTileGph = buildTileGraphForRegistration(tgtTileCoords, refTileCoords)
rgstrCoords = registerWrtRefCycle(tgtImages, refImages, layeredTileGph)

rgstrWrtRefRoot = rgstrTgtCycleRootWrtRefRoot(rgstrCoords, refTileCoords, \
			tgtTileCoords, 'F' + f'{rootForReg:02d}', layeredTileGph)

rgstrWrtRefRoot.to_csv(outFileAllCoords, sep = '\t', index = None)

pd.DataFrame(rgstrWrtRefRoot[['X', 'Y']].apply( \
		lambda x: np.median([y for y in x if abs(y) < maxShift]))).T.to_csv( \
		outFileSmryCoords, sep = '\t', index = None)




