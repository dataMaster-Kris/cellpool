import pandas as pd
import os
import sys
import json
import numpy as np
import skimage.filters

paramsFile = sys.argv[1]
with open(paramsFile, 'r') as inFile:
        params = json.load(inFile)

dat = pd.read_csv('integrated_intensity_table.txt.gz', sep = '\t', compression = 'gzip')
chToElem = pd.read_csv('chToElem.txt', sep = '\t')

bkbnCh = chToElem.channel.loc[chToElem.element == \
	params['debarcodingParams']['backboneProt']].iloc[0]

if (params['debarcodingParams']['backboneProtGatesNoVir'] != ""):
	tcMap = pd.read_csv(os.path.join(params['analysisDir'], \
		params['tcMapFile']), sep = '\t')
	dat['well'] = dat.objId.str[:3]
	noVirDat = dat[dat.well.isin(tcMap.loc[tcMap.poolId == 0].well)]
	
	bkbn = params['debarcodingParams']['backboneProtGatesNoVir']
	bkbnT1 = np.nanquantile(np.log2(noVirDat[bkbnCh].values.tolist()), bkbn[0])
	bkbnT2 = np.nanquantile(np.log2(noVirDat[bkbnCh].values.tolist()), bkbn[1])
	bkbn = pd.DataFrame([2 ** bkbnT1, 2 ** bkbnT2])
else:
	thresholdFunc = getattr(skimage.filters, params['debarcodingParams']['thresholdFunc'])
	bkbnT2 = 2 ** thresholdFunc(np.log2(dat[bkbnCh].values.tolist()))
	bkbnT1 = bkbnT2/2
	bkbn = pd.DataFrame([bkbnT1, bkbnT2])

bkbn.to_csv('backbone_thresholds.txt', header = False, index = False)


