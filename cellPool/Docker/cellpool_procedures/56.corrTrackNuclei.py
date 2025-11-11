import os
import pandas as pd
import numpy as np
import json
import sys

childCyc = sys.argv[1]
well = sys.argv[2]
paramsFile = sys.argv[3]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.tracking import *

iteration = params['segmentationParams']['segIter']
ngbrsToConsider = params['trackingParams']['ngbrsToConsider'] #Tracking labels should be within 5 nearest neigbors in parent cycle
approxmtInfCost = params['trackingParams']['approxmtInfCost'] #used to set infinite costs before running the Kuhn-Munkres algorithm
maxEuclideanDist = params['trackingParams']['maxEuclideanDist'] #Expect label pairs within this distance between cycles
inSuffix = params['trackingParams']['stg1OutSuffix']
outSuffix = params['trackingParams']['stg2OutSuffix']
ftrDir = os.path.join(params['analysisDir'], params['segmentationDir'], \
	'iter' + str(iteration) + params['segmentationParams']['featureDirSuffix'])
files = [x for x in os.listdir(ftrDir) if x.startswith(well)]

trackingDir = os.path.join(params['analysisDir'], params['trackObjsDir'])
parent = params['registrationParams']['registrationPairs']['cycle' + childCyc]

stg1_tracking = pd.read_csv(os.path.join(trackingDir, well + '.' + \
	parent + '_cycle' + childCyc + '_tracking' + inSuffix), sep = "\t")
stg1_ptntl_errors = {}
stg1_ptntl_errors[parent] = pd.read_csv(os.path.join(trackingDir, well + '.' + \
	parent + '_cycle' + childCyc + '_potential_errors' + inSuffix), sep = '\t')
stg1_ptntl_errors['cycle' + childCyc] = stg1_tracking.loc[ \
	stg1_tracking[parent].isin(stg1_ptntl_errors[parent][parent].tolist()) \
				]['cycle' + childCyc].tolist()

stainId = params['segmentationParams']['cellposeParams']['useStain'].replace(' ', '_')
childTbl = pd.read_csv(os.path.join(ftrDir, \
		[x for x in files if x.startswith(well) and ('-cyc' + childCyc in x) and \
			(stainId in x)][0]), \
	sep = "\t")
prntTbl = pd.read_csv(os.path.join(ftrDir, \
		[x for x in files if x.startswith(well) and ('-cyc' + parent[-1] in x) and \
			(stainId in x)][0]), \
	sep = "\t")

trackTbl = trackEuclidean(childTbl, prntTbl, childCyc, parent, \
	stg1_tracking, stg1_ptntl_errors, \
	ngbrsToConsider, approxmtInfCost, maxEuclideanDist) 
trackTbl.to_csv(well + '.' + parent + '_cycle' + childCyc + '_tracking' + outSuffix, \
	sep = '\t', index = False)

