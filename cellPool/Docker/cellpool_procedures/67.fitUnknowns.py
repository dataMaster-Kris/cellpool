import pandas as pd
import sys
import random
import os
import time
import json

start_time = time.time()

cut = int(sys.argv[1])
init = int(sys.argv[2])
paramsFile = sys.argv[3]
with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

sys.path.append(params['reqdCellpoolContainerPath_LET_ME_BE'])
from cellpool.barcode import *

print("Processing cut " + str(cut) + " init " + str(init))

nEpitopesPerBarcode = params['debarcodingParams']['nEpiPerCombo']
bkbn = params['debarcodingParams']['backboneProt']
nCuts = params['debarcodingParams']['nCuts']
inDir = os.path.join(params['analysisDir'], params['barcodeAnalysisDir'])
minPoints4Bckgd = params['debarcodingParams']['minPts4Bckgd']
maxIter = params['debarcodingParams']['maxIter']
resume = params['debarcodingParams']['resumeEM']
resumeFrom = params['debarcodingParams']['resumeEMFrom']
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

intDat = pd.read_csv(os.path.join(inDir, 'cut_' + str(cut) + \
	perCutIntnstyOutFileSuffix), sep = '\t')
bckgdP = pd.read_csv(os.path.join(inDir, 'cut_' + str(cut) + \
	bckgdProbsOutFileSuffix), sep = '\t')
bckgdP.columns = [int(x) for x in bckgdP.columns]

if (resume):
	#---------------------------------
	#Load the parameter states from the last iteration
	#---------------------------------
	pi, delta, mu, Sigma, totalLogL, ix = loadLastIter(inDir, cut, init, \
			nEpitopesPerBarcode, intDat.shape[0], resumeFrom)
	p, e1, e2, logL = {}, {}, {}, {}
	if ((ix == 0) or (ix == -1)):
		resume = False	

if (not resume):
	#---------------------------------
	#Initialization
	#---------------------------------
	epiSignalMeans = pd.read_csv(os.path.join(inDir, \
		'cut_' + str(cut) + epiSignalInitOutFileSuffix), sep = '\t')
	delta, mu, Sigma, pi = {}, {}, {}, {}
	for bcdIx in range(barcodes.shape[0]):
		delta[bcdIx] = np.array(random.sample(range(-10, 10), nEpitopesPerBarcode))
		mu[bcdIx] = np.array([epiSignalMeans[ch].iloc[0] for ch in \
				barcodes.signalChannels.iloc[bcdIx]])
		Sigma[bcdIx] = np.eye(nEpitopesPerBarcode)
		pi[bcdIx] = 1/barcodes.shape[0]

	pd.DataFrame.from_dict({key: [pi[key]] for key in pi}).to_csv( \
		'cut_' + str(cut) + '_init_' + str(init) + '_pi.txt', \
		sep = '\t', index = None)
	pd.DataFrame.from_dict(delta).to_csv( \
		'cut_' + str(cut) + '_init_' + str(init) + '_delta.txt', \
		sep = '\t', index = None)
	pd.DataFrame.from_dict({key: np.ravel(Sigma[key]) for key in Sigma}).to_csv( \
		'cut_' + str(cut) + '_init_' + str(init) + '_Sigma.txt', \
		sep = '\t', index = None)

	#---------------------------------
	#Run expectation-maximization
	#---------------------------------
	p, e1, e2, logL, totalLogL, ix = {}, {}, {}, {}, [], 0

while (ix <= maxIter):
	print([ix, totalLogL[-1]] if len(totalLogL) > 0 else ix)
	#------------------------
	#E-step
	#------------------------
	for bcdIx in range(barcodes.shape[0]):
		print(bcdIx)
		bcd = barcodes.id.iloc[bcdIx]
		p[bcdIx], e1[bcdIx], e2[bcdIx] = fit_mvsnmix_E(intDat, \
				delta[bcdIx], mu[bcdIx], Sigma[bcdIx], bcd, \
				barcodes, bckgdP[bcdIx])

	r = pd.DataFrame(np.column_stack(tuple(pi[x]*p[x] for x in \
		range(barcodes.shape[0]))))
	r = r.divide(r.sum(axis = 1), axis = 0) 

	for bcdIx in range(barcodes.shape[0]):
		print(bcdIx)
		bcd = barcodes.id.iloc[bcdIx]
		#-----------------------
		#M-step
		#-----------------------
		thisDat = intDat[barcodes.signalChannels.iloc[bcdIx]]
		thisDat.columns = [x for x in range(thisDat.shape[1])]
		pi[bcdIx], delta[bcdIx], mu[bcdIx], Sigma[bcdIx] = \
			fit_mvsn_mix_M(thisDat, mu[bcdIx], r[bcdIx], e1[bcdIx], e2[bcdIx])
		#--------------------
		#Compute log likelihood
		logL[bcdIx] = fit_mvsn_mix_logL1(intDat, mu[bcdIx], delta[bcdIx], \
			Sigma[bcdIx], pi[bcdIx], r[bcdIx], bcd, barcodes, bckgdP[bcdIx])

	r['objId'] = intDat.objId
	pd.DataFrame.from_dict({key: [pi[key]] for key in pi}).to_csv( \
		'cut_' + str(cut) + '_init_' + str(init) + '_iter_' + str(ix) + '_pi.txt', \
		sep = '\t', index = None)
	pd.DataFrame.from_dict(mu).to_csv( \
		'cut_' + str(cut) + '_init_' + str(init) + '_iter_' + str(ix) + '_mu.txt', \
		sep = '\t', index = None)
	pd.DataFrame.from_dict(delta).to_csv( \
		'cut_' + str(cut) + '_init_' + str(init) + '_iter_' + str(ix) + '_delta.txt', \
		sep = '\t', index = None)
	pd.DataFrame.from_dict({key: np.ravel(Sigma[key]) for key in Sigma}).to_csv( \
		'cut_' + str(cut) + '_init_' + str(init) + '_iter_' + str(ix) + '_Sigma.txt', \
		sep = '\t', index = None)
	pd.DataFrame.from_dict({key: [logL[key]] for key in logL}).to_csv( \
		'cut_' + str(cut) + '_init_' + str(init) + '_iter_' + str(ix) + '_logL.txt', \
		sep = '\t', index = None)
	r.to_csv('cut_' + str(cut) + '_init_' + str(init) + '_iter_' + str(ix) + '_r.txt', \
		sep = '\t', index = None)

	totalLogL.append(sum([logL[x] for x in range(barcodes.shape[0])]))
	ix += 1

pd.DataFrame(totalLogL).to_csv('cut_' + str(cut) + \
	'_init_' + str(init) + '_totalLogL.txt', sep = '\t', index = None) 

print('Time taken to finish execution: ' + str(time.time() - start_time) + ' seconds')



