import os
import pandas as pd
import numpy as np
import json
import sys

well = sys.argv[1]
paramsFile = sys.argv[2]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.tracking import *

iteration = params['segmentationParams']['segIter']
nNgbrs = params['trackingParams']['nNgbrsToSavePerImage']
ftrDir = os.path.join(params['analysisDir'], params['ashlarAnalysisDir'], \
	params['segmentationDir'], \
        'iter' + str(iteration) + params['segmentationParams']['featureDirSuffix'])
fileName = well + '_' + params['mosaickingParams']['ashlarParamSuffix'] + '.ch0.csv'

tbl = pd.read_csv(os.path.join(ftrDir, fileName), sep = ",")

findNghbrsInSameImage(tbl, nNgbrs).to_csv( \
	well + '_' + params['mosaickingParams']['ashlarParamSuffix'] + \
	'_nearest_nghbrs_stats.txt', sep = '\t', index = False)


