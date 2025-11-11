import os
import numpy as np
import pandas as pd
import json
import sys
import sklearn.linear_model
import networkx as nx

thisCycle = sys.argv[1]
thisWell = sys.argv[2]
paramsFile = sys.argv[3]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.utils import *
from cellpool.mosaic import *

maxDepth = params['mosaickingParams']['maxDepth']
maxPaths = params['mosaickingParams']['maxPaths']
rootForStitching = params['mosaickingParams']['rootForStitching']
rndSeed = params['mosaickingParams']['rndSeed'] #Used to select tiles for fitting with linear model
minPathTgt = params['mosaickingParams']['minPathTgt']
chPriorityOrder = params['priorityChForStitching']['cycle' + thisCycle]

vigCorrDir = os.path.join(params['analysisDir'], params['vigCorrectDir'])
allFiles = os.listdir(vigCorrDir)
allFiles = [x for x in allFiles if \
		x.startswith(thisWell) and \
		('cyc' + thisCycle in x) and \
		x.endswith('.vigCorr.tiff')]
uniqCh = list(set([x[9:] for x in allFiles]))
uniqCh.sort()

outFile = os.path.join(thisWell + '-cyc' + thisCycle + '.tileCentroidCoords.txt')

chInfo = pd.read_csv(os.path.join(params['analysisDir'], params['chToStainId']), sep = '\t')
chInfo = chInfo.loc[chInfo.cycle == 'cycle' + thisCycle]
chInfo = pd.DataFrame.from_dict({x for x in zip([x.strip('.vigCorr.tiff') for x in uniqCh], \
		chInfo[['Channel.' + \
		str(x[0]) for x in uniqCh]].iloc[0].tolist())})
chInfo.columns = ['channel', 'channelName']

priorityOrderChId = [chInfo.channel[chInfo.channelName == x].iloc[0] for \
	x in chPriorityOrder]
priorityOrderChId += list(set(chInfo.channel).difference(priorityOrderChId))
priorityOrderChId = [x for x in priorityOrderChId]
#-----------------------------------------
#Merge tile graphs and coordinate estimates from previous iterations
#-----------------------------------------
tileGph = nx.read_graphml(thisWell + '-C' + \
		priorityOrderChId[0] + '.pruned_tileGph.intmdt.graphml', node_type = int)
if thisWell + '-C' + priorityOrderChId[0] not in params['mosaickingParams']['dfltNgbrLocs']:
    coordsWrtNghbr = pd.read_csv(thisWell+ '-C' + \
                priorityOrderChId[0] + '.coordEsts_est2_filtered.intmdt.txt', sep = '\t')
else:
    coordsWrtNghbr = pd.read_csv(os.path.join(params['analysisDir'], params['mosaicDir'], \
                                    params['mosaickingParams']['mosaicIntermediateFilesDir'], \
                                    params['mosaickingParams']['dfltNgbrLocs'][thisWell + \
                                    '-C' + priorityOrderChId[0]] + '-C' + priorityOrderChId[0] + \
                                    '.coordEsts_est2_filtered.intmdt.txt'), sep = '\t') 
centroidLocs = estimateTileCoordsWrtImRoot(coordsWrtNghbr, tileGph, \
		root = rootForStitching, maxDepth = maxDepth, maxPaths = maxPaths)
prioritySmryCntrdLocs = smrzCoordEsts(centroidLocs, coordsWrtNghbr, rootForStitching)
prioritySmryCntrdLocs.reset_index(drop = True, inplace = True)
prioritySmryCntrdLocs['coordSource'] = chInfo[chInfo.channel == \
		priorityOrderChId[0]]['channelName'].iloc[0]

dctns = ['North', 'East', 'West', 'South']
for nextCh in priorityOrderChId[1:]:
	if not prioritySmryCntrdLocs.medX.isna().any():
		break

	nextTileGph = nx.read_graphml(thisWell + '-C' + \
		nextCh + '.pruned_tileGph.intmdt.graphml', node_type = int) 
	nextCoordsWrtNghbr = pd.read_csv(thisWell + '-C' + \
		nextCh + '.coordEsts_est2_filtered.intmdt.txt', sep = '\t')

	newCoordsWrtNghbr = coordsWrtNghbr.copy()	
	for dctn in dctns:
		for axis in ['X', 'Y']:
			newCoordsWrtNghbr['locEst' + axis + 'WrtNghbrTo' + dctn] = \
				np.where(newCoordsWrtNghbr['locEst' + axis + 'WrtNghbrTo' + \
						dctn].isna(), \
					nextCoordsWrtNghbr['locEst' + axis + 'WrtNghbrTo' + \
						dctn], \
					newCoordsWrtNghbr['locEst' + axis + 'WrtNghbrTo' + dctn])

	fltrdCoordsWrtNghbr = filterUnrealisticEstimates(newCoordsWrtNghbr)
	if (not fltrdCoordsWrtNghbr.locEstXWrtNghbrToWest.isna().all()) and \
		(not fltrdCoordsWrtNghbr.locEstXWrtNghbrToNorth.isna().all()):
		tileGph = nx.compose(tileGph, nextTileGph)
		coordsWrtNghbr = fltrdCoordsWrtNghbr
	
	tileGph = pruneTileGph(tileGph, coordsWrtNghbr)

	centroidLocs = estimateTileCoordsWrtImRoot(coordsWrtNghbr, tileGph, \
			root = rootForStitching, maxDepth = maxDepth, maxPaths = maxPaths)
	smryCntrdLocs = smrzCoordEsts(centroidLocs, coordsWrtNghbr, rootForStitching)
	smryCntrdLocs.reset_index(drop = True, inplace = True)

	smryCntrdLocs['coordSource'] = chInfo[chInfo.channel == nextCh]['channelName'].iloc[0]

	naInPriority = set(prioritySmryCntrdLocs.fieldID[prioritySmryCntrdLocs.medX.isna()])
	fewPathsInPriority = set(prioritySmryCntrdLocs.fieldID[ \
		prioritySmryCntrdLocs.nPaths < minPathTgt])
	morePathsInCurrent = set(smryCntrdLocs.fieldID[ \
		smryCntrdLocs.nPaths >= minPathTgt]).intersection(fewPathsInPriority)
	toReplace = naInPriority.union(morePathsInCurrent)
	
	prioritySmryCntrdLocs.loc[prioritySmryCntrdLocs.fieldID.isin(toReplace)] = \
		smryCntrdLocs.loc[smryCntrdLocs.fieldID.isin(toReplace)].copy()

S = [tileGph.subgraph(c).copy() for c in nx.connected_components(tileGph)]

#----------------------------------------------
#Linear model to link the connected subgraphs to the main graph with field 1
#----------------------------------------------
np.random.seed(rndSeed)
dctnPairs = {'North': 'South', 'East': 'West', 'West': 'East', 'South': 'North'}

while len(S) != 1:
	nodesConnectedToF01 = list([S[x].nodes for x in \
					range(len(S)) if 1 in S[x]][0])
	ngbrsOfWellConnectedNodes = np.unique(coordsWrtNghbr.loc[ \
					coordsWrtNghbr.FieldID.isin(nodesConnectedToF01), \
			['ngbrTo' + x for x in dctns]].to_numpy().flatten())
	ngbrsOfWellConnectedNodes = [int(x) for x in ngbrsOfWellConnectedNodes if \
			(~np.isnan(x)) & (x not in nodesConnectedToF01)]

	tilesToPlaceByLnrModel = [list(set(S[x]).intersection(ngbrsOfWellConnectedNodes)) \
				for x in range(len(S)) if 1 not in S[x]]
	tilesToPlaceByLnrModel = [np.random.choice(x, 1)[0] for x in tilesToPlaceByLnrModel if \
				len(x) > 0]

	datForRegression = coordsWrtNghbr.copy()
	datForRegression = datForRegression.loc[~datForRegression[['locEstXWrtNghbrTo' + x \
			for x in ['North', 'East', 'West', 'South']]].apply( \
			lambda x: np.isnan(x).all(), axis = 1)]

	if (datForRegression.shape[0] > 1):
		for dctn in dctns:
			datForRegression = coordsWrtNghbr.copy()
			datForRegression = datForRegression.loc[ \
				~datForRegression['locEstXWrtNghbrTo' + dctn].isna()]
			lrModel = sklearn.linear_model.LinearRegression()
			lrModel.fit(datForRegression[['PositionX', 'PositionY']], \
				datForRegression[['locEstXWrtNghbrTo' + dctn, \
					'locEstYWrtNghbrTo' + dctn]])

			for field in tilesToPlaceByLnrModel:
				if (~coordsWrtNghbr.loc[coordsWrtNghbr.FieldID == field, \
					'ngbrTo' + dctn].isna().iloc[0]):
					thisPred = [round(x, 1) for x in \
						lrModel.predict(coordsWrtNghbr.loc[ \
						coordsWrtNghbr.FieldID == field, \
						['PositionX', 'PositionY']])[0]]

					coordsWrtNghbr.loc[coordsWrtNghbr.FieldID == field, \
						['locEstXWrtNghbrTo' + dctn, \
						'locEstYWrtNghbrTo' + dctn]] = thisPred
					
					fieldNgbr = coordsWrtNghbr.loc[coordsWrtNghbr.FieldID == \
						field, 'ngbrTo' + dctn].iloc[0]
					coordsWrtNghbr.loc[coordsWrtNghbr.FieldID == fieldNgbr, \
						['locEstXWrtNghbrTo' + dctnPairs[dctn], \
						'locEstYWrtNghbrTo' + dctnPairs[dctn]]] = \
						[-x for x in thisPred]

					tileGph.add_edge(field, int(coordsWrtNghbr.loc[ \
							coordsWrtNghbr.FieldID == field, \
							'ngbrTo' + dctn].iloc[0]))

	S = [tileGph.subgraph(c).copy() for c in nx.connected_components(tileGph)]

#-------------------------------------------
#Graph traversal and estimating tile coordinates wrt f01
#-------------------------------------------
centroidLocs = estimateTileCoordsWrtImRoot(coordsWrtNghbr, tileGph, \
		root = rootForStitching, maxDepth = maxDepth, maxPaths = maxPaths)
smryCntrdLocs = smrzCoordEsts(centroidLocs, coordsWrtNghbr, rootForStitching)
smryCntrdLocs.reset_index(drop = True, inplace = True)

if smryCntrdLocs.medX.isnull().any():
	tileGph = nx.minimum_spanning_tree(tileGph)
	centroidLocs = estimateMstTileCoordsWrtImRoot(coordsWrtNghbr, tileGph, \
		root = rootForStitching)
	smryCntrdLocs = smrzCoordEsts(centroidLocs, coordsWrtNghbr, rootForStitching)
	smryCntrdLocs.reset_index(drop = True, inplace = True)

smryCntrdLocs['coordSource'] = 'LinearModel'

naInPriority = set(prioritySmryCntrdLocs.fieldID[prioritySmryCntrdLocs.medX.isna()])
prioritySmryCntrdLocs.loc[prioritySmryCntrdLocs.fieldID.isin(naInPriority)] = \
	smryCntrdLocs.loc[smryCntrdLocs.fieldID.isin(naInPriority)].copy()

prioritySmryCntrdLocs.to_csv(outFile, sep = "\t", index = False, na_rep='nan')



