import os
import pandas as pd
import numpy as np
import json
import sys
from scipy.spatial import KDTree
from scipy.optimize import linear_sum_assignment

childCyc = sys.argv[1]
well = sys.argv[2]
paramsFile = sys.argv[3]
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
ftrDir = os.path.join(params['analysisDir'], params['segmentationDir'], \
	'iter' + str(iteration) + params['segmentationParams']['featureDirSuffix'])
files = [x for x in os.listdir(ftrDir) if x.startswith(well)]

parent = params['registrationParams']['registrationPairs']['cycle' + childCyc]

stainId = params['segmentationParams']['cellposeParams']['useStain'].replace(' ', '_')
childTbl = pd.read_csv(os.path.join(ftrDir, \
		[x for x in files if x.startswith(well) and ('-cyc' + childCyc in x) and \
			(stainId in x)][0]), \
	sep = "\t")
prntTbl = pd.read_csv(os.path.join(ftrDir, \
		[x for x in files if x.startswith(well) and ('-cyc' + parent[-1] in x) and \
			(stainId in x)][0]), \
	sep = "\t")

trackTbl, nghbrProfiles = trackOrientation(childTbl, prntTbl, [parent, 'cycle' + childCyc], \
	angleDiffCutoff, ngbrsToConsider, approxmtInfCost, maxEuclideanDist)

trackTbl.to_csv(well + '.' + parent + '_cycle' + childCyc + '_tracking' + outSuffix, \
	sep = '\t', index = False)

pd.DataFrame.from_dict(nghbrProfiles).T.to_csv( \
	well + '.' + parent + '_cycle' + childCyc + '_nearest_nghbrs_stats' + outSuffix, \
	sep = '\t', index = False)




