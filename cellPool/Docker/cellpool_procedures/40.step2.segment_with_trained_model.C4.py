import numpy as np
from cellpose import core, utils, io, models, metrics
import sys
import os
import pandas as pd
import json

cycle = sys.argv[1]
well = sys.argv[2]
paramsFile = sys.argv[3]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

mosaicDir = os.path.join(params['analysisDir'], params['mosaicDir'])
modelFile = os.path.join(params['analysisDir'], params['segmentationDir'], \
		'iter' + str(params['segmentationParams']['segIter']) + \
		params['segmentationParams']['trainedModelDirSuffix'], \
		params['segmentationParams']['cellposeParams']['model'])

model = models.CellposeModel(gpu = params['segmentationParams']['useGPU'], \
		pretrained_model = modelFile)

diameter = model.diam_labels

chInfo = pd.read_csv(os.path.join(params['analysisDir'], params['chToStainId']), sep = '\t')
nucCh = chInfo[chInfo.cycle == 'cycle' + str(cycle)].T
nucCh = nucCh[nucCh == \
	params['segmentationParams']['cellposeParams']['useStain']].dropna().index[0][-1] + \
	'-cyc' + str(cycle)

im = np.load(os.path.join(mosaicDir, well + '-C' + nucCh + '.mosaic.trimmed.npy'))
im = np.nan_to_num(im)

#See https://cellpose.readthedocs.io/_/downloads/en/latest/pdf/ ...
#... for parameter interpretations
masks, flows, styles = model.eval(im, \
	channels = params['segmentationParams']['cellposeParams']['channels'], \
	diameter = diameter, \
	flow_threshold = params['segmentationParams']['cellposeParams']['flow_threshold'], \
	cellprob_threshold = params['segmentationParams']['cellposeParams']['cellprob_threshold'])

io.save_masks(im, \
	masks, flows, well + '-cyc' + cycle + '_' + params['segmentationParams']['compartment'], \
	png=True, \
	tif=False, \
	save_flows=False, \
	save_outlines=False, \
	savedir = '.', \
	save_txt = False)


