import os
import numpy as np
import pandas as pd
import sys
import json
from skimage.io import imsave
from sklearn.feature_extraction import image

cycle = sys.argv[1]
well = sys.argv[2]
paramsFile = sys.argv[3]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

thisChannels = params['segmentationParams']['trainingChannels']['cycle' + cycle]

chInfo = pd.read_csv(os.path.join(params['analysisDir'], params['chToStainId']), sep = '\t')
chInfo = chInfo.loc[chInfo.cycle == 'cycle' + cycle]
thisChannelsOrder = [chInfo.T[chInfo.T == x].dropna().index[0][-1] for x in thisChannels]

mosaicDir = os.path.join(params['analysisDir'], params['mosaicDir'])

imFileNames = [well + '-C' + x + '-cyc' + cycle + '.mosaic.trimmed.npy' for x in \
		thisChannelsOrder]

regSize = params['segmentationParams']['patchSizeForTraining']
nRegs = params['segmentationParams']['nPatchesPerWellForTraining']

for thisName in imFileNames:
	thisIm = np.load(os.path.join(mosaicDir, thisName))
	if (thisName == imFileNames[0]):
		im = thisIm
	else:
		im = np.dstack((im, thisIm))

im = np.nan_to_num(im)
patches = image.extract_patches_2d(im, (regSize, regSize), max_patches = nRegs)

for ix in range(len(patches)):
	imsave(well + '-cyc' + cycle + '.patch' + str(ix) + '.tiff', patches[ix])

