import pandas as pd
import sys
import os
from itertools import combinations
import json

paramsFile = sys.argv[1]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

nEpitopesPerBarcode = params['debarcodingParams']['nEpiPerCombo']
bkbn = params['debarcodingParams']['backboneProt']
bcdOutfile = 'barcodes.refmt.txt'

barcodes = pd.read_csv(os.path.join(params['analysisDir'], \
	params['barcodes']), sep = '\t')
barcodes = barcodes[[x for x in barcodes.columns if x != 'id']]
barcodes.columns = ['epitope' + str(x) for x in range(nEpitopesPerBarcode)]
elemCh = pd.read_csv('chToElem.txt', sep = '\t')

elemToCh = dict(zip(elemCh.element, elemCh.channel))
chToElem = dict(zip(elemCh.channel, elemCh.element))

#Assign unique ids to barcodes
barcodes['id'] = 1
barcodes['signalChannels'] = 0
barcodes['noiseChannels'] = 0
barcodes['compatibleStates'] = 0
barcodes['signalChannels'] = barcodes['signalChannels'].astype(dtype = 'object')
barcodes['noiseChannels'] = barcodes['noiseChannels'].astype(dtype = 'object')
barcodes['compatibleStates'] = barcodes['compatibleStates'].astype(dtype = 'object')
for x in range(barcodes.shape[0]):
	barcodes.at[x, 'signalChannels'] = [elemToCh[elem] for elem in barcodes.iloc[x] if \
			elem in elemToCh.keys()]
	barcodes.at[x, 'noiseChannels'] = list(chToElem.keys() - \
			barcodes.at[x, 'signalChannels'] - \
			{elemToCh[bkbn]})
	thisEpitopes = barcodes[['epitope' + str(n) for n in \
		range(nEpitopesPerBarcode)]].iloc[x]
	thisCode = elemCh.element[elemCh.element != bkbn].isin(thisEpitopes).tolist()
	barcodes.loc[x, 'id'] = int(sum([int(thisCode[ix]) * 10**(ix) for ix \
			in range(len(thisCode))]))
	#Find classification states that are compatible with this barcode
	thisCompatibleStates = []
	for nToChoose in range(nEpitopesPerBarcode + 1):
		thisCompatibleStates += [y for y in combinations(thisEpitopes, nToChoose)]
	thisCompatibleCodes = [elemCh.element[elemCh.element != bkbn].isin(y).tolist() \
				for y in thisCompatibleStates]
	thisCompatibleCodes = [int(sum([int(y[ix]) * 10**(ix) for ix \
			in range(len(y))])) for y in thisCompatibleCodes]
	thisCompatibleCodes = [y for y in thisCompatibleCodes if \
		y != barcodes.loc[x, 'id']]
	barcodes.at[x, 'compatibleStates'] = thisCompatibleCodes

barcodes.to_csv(bcdOutfile, sep = '\t', index = None)

