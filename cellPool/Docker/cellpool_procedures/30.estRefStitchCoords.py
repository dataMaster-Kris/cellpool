import os
import numpy as np
import pandas as pd
import skimage.io
from scipy import ndimage
import json
import sys
import networkx as nx

thisCycle = int(sys.argv[1])
thisWell = sys.argv[2]
paramsFile = sys.argv[3]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.utils import *
from cellpool.mosaic import *

maxSeparationBetweenTiles = params['mosaickingParams']['maxTileSepToTestOverlap']
vigCorrDir = os.path.join(params['analysisDir'], params['vigCorrectDir'])

#-----------------------------------------
fieldsData = pd.read_csv(os.path.join(params['analysisDir'], \
			params['wellMapFile']), sep = '\t')
fieldsData = fieldsData.loc[fieldsData.cycle == thisCycle]
fieldsData.FieldID = fieldsData.FieldID.astype('int32')
allFiles = os.listdir(vigCorrDir)
allFiles = [x for x in allFiles if \
		x.startswith(thisWell) and \
		('cyc' + str(thisCycle) in x) and \
		x.endswith('.vigCorr.tiff')]

uniqCh = set([x[9:] for x in allFiles])

for thisCh in uniqCh:
	thisFldsData = fieldsData.copy()
	tileGph = buildNeighborGraphOfTiles(thisFldsData, \
			maxSeparation = maxSeparationBetweenTiles)

	#--------------------------------------
	#Load images
	#--------------------------------------
	thisFiles = [x for x in allFiles if x.endswith(thisCh)]
	images = {}
	for i in thisFiles:
		images[i[4:7]] = skimage.io.imread( \
			os.path.join(vigCorrDir, i))

	#-------------------------------------
	coordEsts = estimateTileCoordsWrtNghbr(images, thisFldsData, tileGph)
	coordEsts.to_csv(thisWell + '-C' + thisCh.strip('.vigCorr.tiff') + \
			'.coordEsts_est1.intmdt.txt', sep = "\t", index = False, na_rep='nan')

	coordEsts = filterUnrealisticEstimates(coordEsts)
	coordEsts.to_csv(thisWell + '-C' + thisCh.strip('.vigCorr.tiff') + \
			'.coordEsts_est1_filtered.intmdt.txt', sep = "\t", \
			index = False, na_rep='nan')

	coordEsts = cropAndEstimateCoordsWrtNghbr(images, coordEsts)
	coordEsts.to_csv(thisWell + '-C' + thisCh.strip('.vigCorr.tiff') + \
			'.coordEsts_est2.intmdt.txt', sep = "\t", index = False, na_rep='nan')

	coordEsts = filterUnrealisticEstimates(coordEsts)
	coordEsts.to_csv(thisWell + '-C' + thisCh.strip('.vigCorr.tiff') + \
			'.coordEsts_est2_filtered.intmdt.txt', sep = "\t", \
			index = False, na_rep='nan')

	tileGph = pruneTileGph(tileGph, coordEsts)
	nx.write_graphml_lxml(tileGph, thisWell + '-C' + thisCh.strip('.vigCorr.tiff') + \
			'.pruned_tileGph.intmdt.graphml')




