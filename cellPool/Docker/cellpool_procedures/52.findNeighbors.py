import os
import pandas as pd
import numpy as np
import json
import sys

cycle = sys.argv[1]
well = sys.argv[2]
paramsFile = sys.argv[3]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.tracking import *

iteration = params['segmentationParams']['segIter']
nNgbrs = params['trackingParams']['nNgbrsToSavePerImage']
ftrDir = os.path.join(params['analysisDir'], params['segmentationDir'], \
        'iter' + str(iteration) + params['segmentationParams']['featureDirSuffix'])
files = [x for x in os.listdir(ftrDir) if x.startswith(well)]
stainId = params['segmentationParams']['cellposeParams']['useStain'].replace(' ', '_')

tbl = pd.read_csv(os.path.join(ftrDir, \
		[x for x in files if x.startswith(well) and ('-cyc' + cycle in x) and \
			(stainId in x)][0]), sep = "\t")

findNghbrsInSameImage(tbl, nNgbrs).to_csv( \
	well + '-cyc' + cycle + '_nearest_nghbrs_stats.txt', \
	sep = '\t', index = False)














