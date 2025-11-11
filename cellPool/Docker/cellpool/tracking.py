import pandas as pd
import numpy as np
from scipy.spatial import KDTree
from scipy.optimize import linear_sum_assignment

def trackOrientation(childTbl, prntTbl, imPairIds, angleDiffCutoff = 0.2, \
	ngbrsToConsider = 5, approxmtInfCost = 1000000, maxEuclideanDist = 50):

	trackTbl = pd.DataFrame()
	tree = KDTree(prntTbl[['centroid-0', 'centroid-1']])
	dd, ii = tree.query(childTbl[['centroid-0', 'centroid-1']], \
		k = ngbrsToConsider)

	#approxmtInfCost approximates infinity in the cost matrix because...
	#... inf values are not allowed by scipy
	cost = np.full((prntTbl.shape[0], childTbl.shape[0]), approxmtInfCost)
	allowedPairs = {}
	nghbrProfiles = {}
	ngbrTblColIds = ['chldLabel']
	ngbrTblColIds.extend(['prntNgbrLabel' + str(x) for x in range(ngbrsToConsider)])
	ngbrTblColIds.extend(['prntNgbrDist' + str(x) for x in range(ngbrsToConsider)])
	ngbrTblColIds.extend(['prntNgbrAnglDiff' + str(x) for x in range(ngbrsToConsider)])
	for ix_c in range(childTbl.shape[0]):
		childFtr = childTbl.orientation.iloc[ix_c]
		thisNgbrs = [childTbl.label.iloc[ix_c]]
		thisNgbrs.extend(prntTbl.label.iloc[ii[ix_c]])
		thisNgbrs.extend(dd[ix_c])
		allowedPairs[childTbl.label.iloc[ix_c]] = \
			prntTbl.label.iloc[[ii[ix_c][x] for x \
				in range(ngbrsToConsider) if \
				dd[ix_c][x] <= maxEuclideanDist]].tolist()
		for ix_p in ii[ix_c]:
			thisCost = prntTbl.orientation.iloc[ix_p] - childFtr
			if thisCost <= angleDiffCutoff:
				cost[ix_p, ix_c] = abs(thisCost)
			thisNgbrs.append(thisCost)
		nghbrProfiles[ix_c] = dict(zip(ngbrTblColIds, thisNgbrs))
	prnt_ind, child_ind = linear_sum_assignment(cost)

	prntLabels = prntTbl.label.iloc[prnt_ind].copy()
	prntLabels.reset_index(drop=True, inplace = True)
	childLabels = childTbl.label.iloc[child_ind].copy()
	childLabels.reset_index(drop=True, inplace = True)

	trackTbl[imPairIds[0]] = prntLabels
	trackTbl[imPairIds[1]] = childLabels
	rowsWithAllowedPairs = []
	for pairIx in range(trackTbl.shape[0]):
		if trackTbl[imPairIds[0]].iloc[pairIx] in \
			allowedPairs[trackTbl[imPairIds[1]].iloc[pairIx]]:
			rowsWithAllowedPairs.append(pairIx)

	trackTbl = trackTbl.iloc[rowsWithAllowedPairs].copy()
	trackTbl.reset_index(drop=True, inplace = True)
	return([trackTbl, nghbrProfiles])


def findNghbrsInSameImage(tbl, nNgbrs = 31):
	tree = KDTree(tbl[['centroid-0', 'centroid-1']])
	dd, ii = tree.query(tbl[['centroid-0', 'centroid-1']], \
		k = nNgbrs)

	nghbrProfiles = {}
	ngbrTblColIds = ['label']
	ngbrTblColIds.extend(['ngbrLabel' + str(x) for x in range(nNgbrs)])
	ngbrTblColIds.extend(['ngbrDist' + str(x) for x in range(nNgbrs)])
	ngbrTblColIds.extend(['ngbrAnglDiff' + str(x) for x in range(nNgbrs)])
	for ix_c in range(tbl.shape[0]):
		ftr = tbl.orientation.iloc[ix_c]
		thisNgbrs = [tbl.label.iloc[ix_c]]
		thisNgbrs.extend(tbl.label.iloc[ii[ix_c]])
		thisNgbrs.extend(dd[ix_c])
		for ix_p in ii[ix_c]:
			thisCost = tbl.orientation.iloc[ix_p] - ftr
			thisNgbrs.append(thisCost)
		nghbrProfiles[ix_c] = dict(zip(ngbrTblColIds, thisNgbrs))
	
	return(pd.DataFrame.from_dict(nghbrProfiles).T)


def trackEuclidean(childTbl, prntTbl, childCycle, parentCycle, \
	trackMatesStage1, errorsStage1, ngbrsToConsider = 5, \
	approxmtInfCost = 1000000, maxEuclideanDist = 50):
	trackTbl = pd.DataFrame()
	tree = KDTree(prntTbl[['centroid-0', 'centroid-1']])
	dd, ii = tree.query(childTbl[['centroid-0', 'centroid-1']], \
		k = ngbrsToConsider)

	#1000 approximates infinity in the cost matrix because...
	#... inf values are not allowed by scipy
	cost = np.full((prntTbl.shape[0], childTbl.shape[0]), approxmtInfCost)
	allowedPairs = {}
	validNewParentLabels = set(prntTbl.label.tolist()).difference( \
			trackMatesStage1[parentCycle]).union( \
			errorsStage1[parentCycle][parentCycle].tolist())
	for ix_c in range(childTbl.shape[0]):
		childFtr = childTbl.orientation.iloc[ix_c]
		childX = childTbl['centroid-1'].iloc[ix_c]
		childY = childTbl['centroid-0'].iloc[ix_c]
		if (childTbl.label.iloc[ix_c] in \
			trackMatesStage1['cycle' + childCycle].tolist()) & \
			(childTbl.label.iloc[ix_c] not in errorsStage1['cycle' + childCycle]):
			allowedPairs[childTbl.label.iloc[ix_c]] = \
				trackMatesStage1[trackMatesStage1['cycle' + childCycle] == \
					childTbl.label.iloc[ix_c]][parentCycle].tolist()
		else:
			allowedPairs[childTbl.label.iloc[ix_c]] = \
				set(prntTbl.label.iloc[[ii[ix_c][x] for x \
					in range(ngbrsToConsider) if \
					dd[ix_c][x] <= maxEuclideanDist \
					]].tolist()).intersection(validNewParentLabels)
		for ix_p in ii[ix_c]:
			if (childTbl.label.iloc[ix_c] not in errorsStage1['cycle' + childCycle]):
				if (childTbl.label.iloc[ix_c] not in \
						trackMatesStage1['cycle' + childCycle].tolist()):
					thisCost = approxmtInfCost
				elif (trackMatesStage1.loc[trackMatesStage1['cycle' + childCycle] == \
					childTbl.label.iloc[ix_c]][parentCycle].iloc[0] == \
						prntTbl.label.iloc[ix_p]):
					thisCost = 0
				else:
					thisCost = approxmtInfCost
			else:
				if (prntTbl.label.iloc[ix_p] in \
					errorsStage1[parentCycle][parentCycle].tolist()):
					expChildX = errorsStage1[parentCycle]['expDispX'].loc[ \
						errorsStage1[parentCycle][parentCycle] == \
							prntTbl.label.iloc[ix_p]].iloc[0] + \
						prntTbl['centroid-1'].iloc[ix_p]
					expChildY = errorsStage1[parentCycle]['expDispY'].loc[ \
						errorsStage1[parentCycle][parentCycle] == \
							prntTbl.label.iloc[ix_p]].iloc[0] + \
						prntTbl['centroid-0'].iloc[ix_p]
					thisCost = ((childX - expChildX)**2 + \
						(childY - expChildY)**2)**0.5
				else:
					thisCost = approxmtInfCost
			cost[ix_p, ix_c] = abs(thisCost)
	prnt_ind, child_ind = linear_sum_assignment(cost)

	prntLabels = prntTbl.label.iloc[prnt_ind].copy()
	prntLabels.reset_index(drop=True, inplace = True)
	childLabels = childTbl.label.iloc[child_ind].copy()
	childLabels.reset_index(drop=True, inplace = True)

	trackTbl[parentCycle] = prntLabels
	trackTbl['cycle' + childCycle] = childLabels
	rowsWithAllowedPairs = []
	for pairIx in range(trackTbl.shape[0]):
		if trackTbl[parentCycle].iloc[pairIx] in \
			allowedPairs[trackTbl['cycle' + childCycle].iloc[pairIx]]:
			rowsWithAllowedPairs.append(pairIx)
	trackTbl = trackTbl.iloc[rowsWithAllowedPairs].copy()
	trackTbl.reset_index(drop=True, inplace = True)
	return(trackTbl)

def trackLocalRef(childTbl, prntTbl, childCycle, parentCycle, trackTbl, \
	ngbrsToConsider = 5, approxmtInfCost = 1000000, sqrDistCutoff = 9):
	childTree = KDTree(childTbl[['centroid-0', 'centroid-1']])
	prntTree = KDTree(prntTbl[['centroid-0', 'centroid-1']])

	while ~np.all(trackTbl.state == 'U') & (trackTbl.shape[0] != childTbl.shape[0]):
		childRef, prntRef = trackTbl[trackTbl.state == \
					'N'].sample()[[childCycle, parentCycle]].iloc[0]

		pRefFtr = prntTbl[prntTbl.label == prntRef]
		cRefFtr = childTbl[childTbl.label == childRef]

		cdd, cii = childTree.query(cRefFtr[['centroid-0', 'centroid-1']], \
			k = ngbrsToConsider)
		pdd, pii = prntTree.query(pRefFtr[['centroid-0', 'centroid-1']], \
			k = ngbrsToConsider)

		pNgbrs = prntTbl.iloc[pii[0]].copy()
		cNgbrs = childTbl.iloc[cii[0]].copy()
		pNgbrs = pNgbrs[~pNgbrs.label.isin(trackTbl[parentCycle])]
		cNgbrs = cNgbrs[~cNgbrs.label.isin(trackTbl[childCycle])]

		pNgbrs[['centroid-0', 'centroid-1']] -= \
			pRefFtr[['centroid-0', 'centroid-1']].iloc[0].tolist()
		cNgbrs[['centroid-0', 'centroid-1']] -= \
			cRefFtr[['centroid-0', 'centroid-1']].iloc[0].tolist()

		cost = np.full((pNgbrs.shape[0], cNgbrs.shape[0]), approxmtInfCost)

		for ix_c in range(cNgbrs.shape[0]):
			childX = cNgbrs['centroid-1'].iloc[ix_c]
			childY = cNgbrs['centroid-0'].iloc[ix_c]
			for ix_p in range(pNgbrs.shape[0]):
				prntX = pNgbrs['centroid-1'].iloc[ix_p]
				prntY = pNgbrs['centroid-0'].iloc[ix_p]
				thisCost = (prntX - childX) ** 2 + (prntY - childY) ** 2
				cost[ix_p, ix_c] = thisCost

		prnt_ind, child_ind = linear_sum_assignment(cost)

		#Remove all assignments that do not meet the sqrDistCutoff criteria
		asgnmtCosts = [cost[x[0], x[1]] for x in zip(prnt_ind, child_ind)]
		toRemove = [x for x in range(len(asgnmtCosts)) if asgnmtCosts[x] > sqrDistCutoff]
		prnt_ind = list(prnt_ind)
		child_ind = list(child_ind)
		for i in sorted(toRemove, reverse=True):
			del prnt_ind[i]
			del child_ind[i]

		newTracked = pd.DataFrame()
		newTracked[parentCycle] = pNgbrs.label.iloc[prnt_ind].tolist()
		newTracked[childCycle] = cNgbrs.label.iloc[child_ind].tolist()
		newTracked['state'] = 'N'

		trackTbl.loc[trackTbl[parentCycle] == prntRef, 'state'] = 'U'
		trackTbl = pd.concat([trackTbl, newTracked], ignore_index = True)

	return(trackTbl)

def localNgbrAlgnCost(childTbl, prntTbl, childCycle, parentCycle, trackTbl, \
	ngbrsToConsider = 5, approxmtInfCost = 1000000):

	childTbl.columns = ['label', 'centroid-0', 'centroid-1']
	prntTbl.columns = ['label', 'centroid-0', 'centroid-1']
	childTbl = childTbl[childTbl.label.isin(trackTbl[childCycle])]
	prntTbl = prntTbl[prntTbl.label.isin(trackTbl[parentCycle])]
	childTbl.reset_index(inplace = True, drop = True)
	prntTbl.reset_index(inplace = True, drop = True)

	childTree = KDTree(childTbl[['centroid-0', 'centroid-1']])
	prntTree = KDTree(prntTbl[['centroid-0', 'centroid-1']])
	costs = []
	for ix in range(trackTbl.shape[0]):
		childRef, prntRef = trackTbl[[childCycle, parentCycle]].iloc[ix]

		pRefFtr = prntTbl[prntTbl.label == prntRef]
		cRefFtr = childTbl[childTbl.label == childRef]

		cdd, cii = childTree.query(cRefFtr[['centroid-0', 'centroid-1']], \
			k = ngbrsToConsider)
		cNgbrs = childTbl.iloc[cii[0]].copy()
		cNgbrs[['centroid-0', 'centroid-1']] -= \
			cRefFtr[['centroid-0', 'centroid-1']].iloc[0].tolist()

		thisCost = 0
		for ix_c in range(cNgbrs.shape[0]):
			childLabel = cNgbrs['label'].iloc[ix_c]
			childX = cNgbrs['centroid-1'].iloc[ix_c]
			childY = cNgbrs['centroid-0'].iloc[ix_c]
			prntLabel = trackTbl[trackTbl[childCycle] == childLabel \
						][parentCycle].iloc[0]
			prntX = prntTbl[prntTbl.label == prntLabel]['centroid-1'].iloc[0] - \
				pRefFtr['centroid-1'].iloc[0]
			prntY = prntTbl[prntTbl.label == prntLabel]['centroid-0'].iloc[0] - \
				pRefFtr['centroid-0'].iloc[0]
			thisCost += (prntX - childX) ** 2 + (prntY - childY) ** 2
		
		costs.append(thisCost/(ngbrsToConsider - 1))

	return(costs)

