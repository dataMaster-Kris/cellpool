import os
import pandas as pd
import numpy as np
import json
import sys
from scipy.spatial import KDTree
from scipy.optimize import linear_sum_assignment

well = sys.argv[1]
paramsFile = sys.argv[2]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.tracking import *

iteration = params['segmentationParams']['segIter']
outSuffix = params['trackingParams']['stg1OutSuffix']
angleDiffCutoff = params['trackingParams']['angleDiffCutoff'] #in radians
ngbrsToConsider = params['trackingParams']['ngbrsToConsider'] #Tracking labels should be within 5 nearest neigbors in parent cycle
approxmtInfCost = params['trackingParams']['approxmtInfCost'] #used to set infinite costs before running the Kuhn-Munkres algorithm
maxEuclideanDist = params['trackingParams']['maxEuclideanDist'] #Expect label pairs within this distance between cycles
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

trackTbl, nghbrProfiles = trackOrientation(childTbl, prntTbl, ['cellpool', 'ashlar'], \
	angleDiffCutoff, ngbrsToConsider, approxmtInfCost, maxEuclideanDist)

trackTbl.to_csv(well + '.' + 'cellpool_ashlar_tracking' + outSuffix, \
	sep = '\t', index = False)

pd.DataFrame.from_dict(nghbrProfiles).T.to_csv( \
	well + '.' + 'cellpool_ashlar_nearest_nghbrs_stats' + outSuffix, \
	sep = '\t', index = False)




