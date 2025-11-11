import pandas as pd
from sklearn.cluster import KMeans
import sys
import random
import os
from itertools import combinations
import json

thisCut = int(sys.argv[1])
paramsFile = sys.argv[2]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.barcode import *

nEpitopesPerBarcode = params['debarcodingParams']['nEpiPerCombo']
bkbn = params['debarcodingParams']['backboneProt']
nCuts = params['debarcodingParams']['nCuts']
inDir = os.path.join(params['analysisDir'], params['barcodeAnalysisDir'])
minPoints4Bckgd = params['debarcodingParams']['minPts4Bckgd']
perCutIntnstyOutFileSuffix = '_intDat.txt'
bckgdProbsOutFileSuffix = '_background_probs_per_barcode.txt'
epiSignalInitOutFileSuffix = '_init_mu.txt'

intnsty = pd.read_csv(os.path.join(inDir, 'integrated_intensity_table.txt.gz'), sep = '\t')
calls = pd.read_csv(os.path.join(inDir, 'calls_with_pop_level_epitope_gating.txt.gz'), sep = '\t')
barcodes = pd.read_csv(os.path.join(inDir, 'barcodes.refmt.txt'), sep = '\t')
elemCh = pd.read_csv(os.path.join(inDir, 'chToElem.txt'), sep = '\t')

elemToCh = dict(zip(elemCh.element, elemCh.channel))
chToElem = dict(zip(elemCh.channel, elemCh.element))
barcodes.signalChannels = barcodes.signalChannels.apply(eval)
barcodes.noiseChannels = barcodes.noiseChannels.apply(eval)
barcodes.compatibleStates = barcodes.compatibleStates.apply(eval)

intnsty = intnsty[['objId'] + elemCh.channel.tolist()]
intnsty[elemCh.channel.tolist()] = np.log2(intnsty[elemCh.channel.tolist()])
uninfctdIntnsty = intnsty[calls[elemToCh[bkbn]] == 0]
uninfctdCalls = calls[calls[elemToCh[bkbn]] == 0]
infctdIntnsty = intnsty[calls[elemToCh[bkbn]] == 1]
infctdCalls = calls[calls[elemToCh[bkbn]] == 1]

bkbnQcut = pd.qcut(infctdIntnsty[elemToCh[bkbn]], q = nCuts, labels = range(1, nCuts + 1))

intDat = infctdIntnsty[bkbnQcut == thisCut].copy()
callDat = infctdCalls[bkbnQcut == thisCut].copy()
intDat['classificationStatus'] = callDat[elemCh.channel[elemCh.element != bkbn]].apply( \
	lambda x: classificationStatus(x, nEpitopesPerBarcode), axis = 1)
intDat.drop(intDat[intDat.classificationStatus == -999].index, inplace = True)
intDat.dropna(inplace = True)
intDat.reset_index(drop = True, inplace = True)
intDat.to_csv('cut_' + str(thisCut) + \
	perCutIntnstyOutFileSuffix, sep = '\t', index = None)
#---------------------------------
#Get background distributions independent of the barcode
#---------------------------------
bckgdDists = {}
bckgdP = {}
for bcdIx in range(barcodes.shape[0]):
	thisNoiseChannels = barcodes.noiseChannels.iloc[bcdIx]
	thisDat4BckgdCalc = intDat[intDat.classificationStatus == \
		barcodes.id.iloc[bcdIx]]
	if (thisDat4BckgdCalc.shape[0] < minPoints4Bckgd):
		print(bcdIx)
		thisDat4BckgdCalc = pd.concat([uninfctdIntnsty, \
			intDat[intDat.classificationStatus == barcodes.id.iloc[bcdIx]]], \
			ignore_index = True)
	bckgdDists[bcdIx] = getEmpiricalProbabilityFunction( \
				thisDat4BckgdCalc[thisNoiseChannels])
	bckgdP[bcdIx] = bckgdPDF(intDat[thisNoiseChannels], bckgdDists[bcdIx])

bckgdP = pd.DataFrame.from_dict(bckgdP)
bckgdP.to_csv('cut_' + str(thisCut) + bckgdProbsOutFileSuffix, \
	sep = '\t', index = None)
#---------------------------------
#Initialization
#---------------------------------
epiSignalMeans = {}
for epi in elemToCh:
	if epi == bkbn:
		continue
	kmeans = KMeans(n_clusters = 2, random_state = 0, \
		n_init = 10).fit(np.array(intDat[elemToCh[epi]]).reshape(-1, 1))
	epiSignalMeans[elemToCh[epi]] = kmeans.cluster_centers_.max()

pd.DataFrame.from_dict({key:[epiSignalMeans[key]] for key \
	in epiSignalMeans}).to_csv('cut_' + str(thisCut) + epiSignalInitOutFileSuffix, \
	sep = '\t', index = None)



