import os
import pandas as pd
import numpy as np
import json
import sys
from scipy.spatial import KDTree
from scipy.optimize import linear_sum_assignment

childCyc = 'cycle' + str(sys.argv[1])
well = sys.argv[2]
paramsFile = sys.argv[3]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.tracking import *

iteration = params['segmentationParams']['segIter']
outSuffix = params['trackingParams']['stg3OutSuffix']
sqrDistCutoff = params['trackingParams']['localRefParams']['sqrDistCutoff']
ngbrsToConsider = params['trackingParams']['localRefParams']['ngbrsToConsider'] #Tracking labels should be within N nearest neigbors in parent cycle
approxmtInfCost = params['trackingParams']['localRefParams']['approxmtInfCost']
trackingQcDir = params['trackingParams']['localRefParams']['trackingQcDir']

ftrDir = os.path.join(params['analysisDir'], params['segmentationDir'], \
	'iter' + str(iteration) + params['segmentationParams']['featureDirSuffix'])
files = [x for x in os.listdir(ftrDir) if x.startswith(well + '-cyc') and \
	x.endswith(params['segmentationParams']['cellposeParams'][ \
		'useStain'].replace(' ', '_') + '.txt')]

validMates = pd.read_csv(os.path.join(params['analysisDir'], params['trackObjsDir'], \
	well + '.highConfValidMates.txt'), sep = '\t')
parent = params['registrationParams']['registrationPairs'][childCyc]
trackTbl = validMates[[parent, childCyc]].copy()
trackTbl['state'] = 'N' #U: used as ref already #N: not used yet

childTbl = pd.read_csv(os.path.join(ftrDir, \
		[x for x in files if x.startswith(well + '-cyc' + childCyc[-1])][0]), \
	sep = "\t")
prntTbl = pd.read_csv(os.path.join(ftrDir, \
		[x for x in files if x.startswith(well + '-cyc' + parent[-1])][0]), \
	sep = "\t")

trackTbl = trackLocalRef(childTbl, prntTbl, childCyc, parent, trackTbl, \
	ngbrsToConsider, approxmtInfCost, sqrDistCutoff)

trackTbl[[parent, childCyc]].astype(np.int64).to_csv( \
		well + '.' + parent + '_' + childCyc + '_tracking' + \
		outSuffix, \
		sep = '\t', index = False)

