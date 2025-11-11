import skimage.io
import skimage.exposure
import os
import json
import numpy as np
import pandas as pd
import sys

cycle = int(sys.argv[1])
field = int(sys.argv[2])
paramsFile = sys.argv[3]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.utils import *

cycleInformation = getCycleDirMatch(params['rawDataDir'], \
	params['cycleIdDelimiters'], params['cycleIdFormat'], \
	params['keepDirIfPrefix'], params['keepDirIfSuffix'], params['keepDirIfContains'], \
	params['throwDirIfPrefix'], params['throwDirIfSuffix'], params['throwDirIfContains'])

rawDirPath = os.path.join(params['rawDataDir'], \
		cycleInformation['dirs'][cycleInformation.cycle == cycle].iloc[0], \
		params['imgPathRelativeToCycDir'])
imFiles = os.listdir(rawDirPath)
imFiles = [x for x in imFiles if x.endswith(params['imgFormat'])]
imFilesInfo = getImageFileInfo(imFiles, params['fileNameFormat'])
imFiles = imFilesInfo[imFilesInfo.Field == field]

with open(os.path.join(params['analysisDir'], params['listOfWells']), 'r') as fopen:
	 wells = fopen.readlines()
wells = [x.rstrip('\n') for x in wells]

for thisWell in wells:
	thisChannels = np.unique(imFiles.Channel)
	thisPlanes = np.unique(imFiles.Z)
	for channel in thisChannels:
		toMaxProject = imFiles.filename[(imFiles.Well == thisWell) & \
				(imFiles.Channel == channel)].tolist()
		im = np.expand_dims(skimage.io.imread(fname = os.path.join(rawDirPath, \
				toMaxProject[0])), axis = 2)
		for anotherZPlane in toMaxProject[1:]:
			imThis = skimage.io.imread(fname = \
					os.path.join(rawDirPath, \
					anotherZPlane))
			im = np.dstack((im, imThis))
		im = np.max(im, axis = 2)
		skimage.io.imsave(fname = os.path.join(thisWell + '_F' + f"{field:02d}" + \
				'-C' + str(channel)) + '-cyc' + str(cycle) + '.tiff', \
				arr = im, check_contrast = False)


