import os
import pandas as pd
import numpy as np
import json
import sys
from scipy.spatial import KDTree
from scipy.optimize import linear_sum_assignment

childCyc = "ashlar"
well = sys.argv[1]
paramsFile = sys.argv[2]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

#sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
#from cellpool.tracking import *

iteration = params['segmentationParams']['segIter']
ngbrsToConsider = params['trackingParams']['ngbrsToConsider'] #Tracking labels should be within 5 nearest neigbors in parent cycle
approxmtInfCost = params['trackingParams']['approxmtInfCost'] #used to set infinite costs before running the Kuhn-Munkres algorithm
maxEuclideanDist = params['trackingParams']['maxEuclideanDist'] #Expect label pairs within this distance between cycles
inSuffix = params['trackingParams']['stg1OutSuffix']
outSuffix = params['trackingParams']['stg2OutSuffix']
stainId = params['segmentationParams']['cellposeParams']['useStain'].replace(' ', '_')
ftrDir_cp = os.path.join(params['analysisDir'], params['segmentationDir'], \
	'iter' + str(iteration) + params['segmentationParams']['featureDirSuffix'])
ftrDir_ashlr = os.path.join(params['analysisDir'], params['ashlarAnalysisDir'], \
	params['segmentationDir'], 'iter' + str(iteration) + \
	params['segmentationParams']['featureDirSuffix'])
fileName_cp = well + '-cyc1.' + stainId + '.txt'
fileName_ashlr = well + '_' + params['mosaickingParams']['ashlarParamSuffix'] + '.ch0.csv'

childTbl = pd.read_csv(os.path.join(ftrDir_ashlr, fileName_ashlr), sep = ",")
prntTbl = pd.read_csv(os.path.join(ftrDir_cp, fileName_cp), sep = "\t")

trackingDir = os.path.join(params['analysisDir'], params['ashlarAnalysisDir'], \
		'45.track_cellpool_to_ashlar')
parent = "cellpool" 

stg1_tracking = pd.read_csv(os.path.join(trackingDir, well + '.' + \
	parent + '_' + childCyc + '_tracking' + inSuffix), sep = "\t")
stg1_ptntl_errors = {}
stg1_ptntl_errors[parent] = pd.read_csv(os.path.join(trackingDir, well + '.' + \
	parent + '_' + childCyc + '_potential_errors' + inSuffix), sep = '\t')
stg1_ptntl_errors[childCyc] = stg1_tracking.loc[ \
	stg1_tracking[parent].isin(stg1_ptntl_errors[parent][parent].tolist()) \
				][childCyc].tolist()

#--------------------------------------------------
def trackEuclidean(childTbl, prntTbl, childCycle, parentCycle, \
	trackMatesStage1, errorsStage1, ngbrsToConsider = 5, \
	approxmtInfCost = 1000000, maxEuclideanDist = 50):
	trackTbl = pd.DataFrame()
	tree = KDTree(prntTbl[['centroid-0', 'centroid-1']])
	dd, ii = tree.query(childTbl[['centroid-0', 'centroid-1']], \
		k = ngbrsToConsider)

	#1000 approximates infinity in the cost matrix because...
	#... inf values are not allowed by scipy
	cost = np.full((prntTbl.shape[0], childTbl.shape[0]), approxmtInfCost)
	allowedPairs = {}
	validNewParentLabels = set(prntTbl.label.tolist()).difference( \
			trackMatesStage1[parentCycle]).union( \
			errorsStage1[parentCycle][parentCycle].tolist())
	for ix_c in range(childTbl.shape[0]):
		childFtr = childTbl.orientation.iloc[ix_c]
		childX = childTbl['centroid-1'].iloc[ix_c]
		childY = childTbl['centroid-0'].iloc[ix_c]
		if (childTbl.label.iloc[ix_c] in \
			trackMatesStage1[childCycle].tolist()) & \
			(childTbl.label.iloc[ix_c] not in errorsStage1[childCycle]):
			allowedPairs[childTbl.label.iloc[ix_c]] = \
				trackMatesStage1[trackMatesStage1[childCycle] == \
					childTbl.label.iloc[ix_c]][parentCycle].tolist()
		else:
			allowedPairs[childTbl.label.iloc[ix_c]] = \
				set(prntTbl.label.iloc[[ii[ix_c][x] for x \
					in range(ngbrsToConsider) if \
					dd[ix_c][x] <= maxEuclideanDist \
					]].tolist()).intersection(validNewParentLabels)
		for ix_p in ii[ix_c]:
			if (childTbl.label.iloc[ix_c] not in errorsStage1[childCycle]):
				if (childTbl.label.iloc[ix_c] not in \
						trackMatesStage1[childCycle].tolist()):
					thisCost = approxmtInfCost
				elif (trackMatesStage1.loc[trackMatesStage1[childCycle] == \
					childTbl.label.iloc[ix_c]][parentCycle].iloc[0] == \
						prntTbl.label.iloc[ix_p]):
					thisCost = 0
				else:
					thisCost = approxmtInfCost
			else:
				if (prntTbl.label.iloc[ix_p] in \
					errorsStage1[parentCycle][parentCycle].tolist()):
					expChildX = errorsStage1[parentCycle]['expDispX'].loc[ \
						errorsStage1[parentCycle][parentCycle] == \
							prntTbl.label.iloc[ix_p]].iloc[0] + \
						prntTbl['centroid-1'].iloc[ix_p]
					expChildY = errorsStage1[parentCycle]['expDispY'].loc[ \
						errorsStage1[parentCycle][parentCycle] == \
							prntTbl.label.iloc[ix_p]].iloc[0] + \
						prntTbl['centroid-0'].iloc[ix_p]
					thisCost = ((childX - expChildX)**2 + \
						(childY - expChildY)**2)**0.5
				else:
					thisCost = approxmtInfCost
			cost[ix_p, ix_c] = abs(thisCost)
	prnt_ind, child_ind = linear_sum_assignment(cost)

	prntLabels = prntTbl.label.iloc[prnt_ind].copy()
	prntLabels.reset_index(drop=True, inplace = True)
	childLabels = childTbl.label.iloc[child_ind].copy()
	childLabels.reset_index(drop=True, inplace = True)

	trackTbl[parentCycle] = prntLabels
	trackTbl[childCycle] = childLabels
	rowsWithAllowedPairs = []
	for pairIx in range(trackTbl.shape[0]):
		if trackTbl[parentCycle].iloc[pairIx] in \
			allowedPairs[trackTbl[childCycle].iloc[pairIx]]:
			rowsWithAllowedPairs.append(pairIx)
	trackTbl = trackTbl.iloc[rowsWithAllowedPairs].copy()
	trackTbl.reset_index(drop=True, inplace = True)
	return(trackTbl)

#---------------------------------------------------

trackTbl = trackEuclidean(childTbl, prntTbl, childCyc, parent, \
	stg1_tracking, stg1_ptntl_errors, \
	ngbrsToConsider, approxmtInfCost, maxEuclideanDist) 
trackTbl.to_csv(well + '.' + parent + '_' + childCyc + '_tracking' + outSuffix, \
	sep = '\t', index = False)

