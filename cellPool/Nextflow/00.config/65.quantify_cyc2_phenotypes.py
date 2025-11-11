import os
import numpy as np
import pandas as pd
import sys
from PIL import Image
Image.MAX_IMAGE_PIXELS = 1000000000 #Turn off DecompressionBombWarning for large mosaics
import imageio
from skimage.segmentation import expand_labels
from skimage.measure import regionprops_table
import json

cyc = '2'
well = sys.argv[1]
paramsFile = sys.argv[2]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.utils import *
from cellpool.mosaic import *

wellBorderWidth = params['featureExtractionParams']['wellBorderWidth']
mosaicDir = os.path.join(params['analysisDir'], params['mosaicDir'])
maskDir = os.path.join(params['analysisDir'], params['segmentationDir'], \
	'iter' + str(params['segmentationParams']['segIter']) + \
	params['segmentationParams']['maskDirSuffix'])
properties = tuple(params['featureExtractionParams']['regionProperties'])

#----------------------------------------------------
thisMaskExt = os.path.splitext( \
        [x for x in os.listdir(maskDir) if x.startswith(well) and \
                ('-cyc' + str(cyc) in x) and \
                ('_cp_masks' in x)][0])[1]
mask = imageio.imread(os.path.join(maskDir, well + '-cyc' + cyc + '_' + \
                params['segmentationParams']['compartment'] + '_cp_masks' + thisMaskExt))
mask2 = expand_labels(mask, distance = 3) - mask
chInfo = pd.read_csv(os.path.join(params['analysisDir'], params['chToStainId']), sep = '\t')

for nextCh in range(1, 5):
        chIm = np.load(os.path.join(mosaicDir, well + '-C' + str(nextCh) + \
                                                '-cyc' + str(cyc) + '.mosaic.trimmed.npy'))
        props = regionprops_table(mask2 if (nextCh == 1) else mask, intensity_image = chIm, \
                                properties = properties, extra_properties = (distribSummarize,))
        keysOfProps = [x for x in props if x != 'distribSummarize']
        outTbl = pd.concat([pd.DataFrame({key: value for key, value \
                                    in props.items() \
                                    if key in keysOfProps}), \
                            pd.DataFrame(props['distribSummarize'].tolist())], axis = 1)
        outTbl.to_csv(os.path.join(well + '-cyc' + str(cyc) + \
                '.' + 'ch' + str(nextCh) + '_phenotypic.txt'), \
                sep = '\t', index = False)






