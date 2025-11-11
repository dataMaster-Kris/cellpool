import numpy as np
from cellpose import core, utils, io, models, metrics
import sys
import os
import pandas as pd
import json
from tifffile import imread

wellOmeTiff = sys.argv[1]
paramsFile = sys.argv[2]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

modelFile = os.path.join(params['analysisDir'], params['segmentationDir'], \
		'iter' + str(params['segmentationParams']['segIter']) + \
		params['segmentationParams']['trainedModelDirSuffix'], \
		params['segmentationParams']['cellposeParams']['model'])

model = models.CellposeModel(gpu = params['segmentationParams']['useGPU'], \
		pretrained_model = modelFile)

diameter = model.diam_labels

im = imread(wellOmeTiff)
im = np.nan_to_num(im[0, :, :])

#See https://cellpose.readthedocs.io/_/downloads/en/latest/pdf/ ...
#... for parameter interpretations
masks, flows, styles = model.eval(im, \
	channels = params['segmentationParams']['cellposeParams']['channels'], \
	diameter = diameter, \
	flow_threshold = params['segmentationParams']['cellposeParams']['flow_threshold'], \
	cellprob_threshold = params['segmentationParams']['cellposeParams']['cellprob_threshold'])

io.save_masks(im, \
	masks, flows, wellOmeTiff.replace('.ome.tiff', '') + \
		'-cyc1_' + params['segmentationParams']['compartment'], \
	png=True, \
	tif=False, \
	save_flows=False, \
	save_outlines=False, \
	savedir = '.', \
	save_txt = False)


