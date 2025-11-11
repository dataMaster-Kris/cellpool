import numpy as np
from scipy.stats import multivariate_normal as mvn, norm
from scipy.stats._multivariate import _squeeze_output
import pandas as pd
import os

#Direct parametrization without reference to a generative process
def mvsn_pdf1(x, shape, mu = None, cov = None):
	#Make numpy arrays storing the inputs
	dmnsn = len(shape)
	shape = np.asarray(shape)
	mu = np.zeros(dmnsn) if mu is None else np.asarray(mu)
	x_centered = x - mu
	cov = np.eye(dmnsn) if cov is None else np.asarray(cov)

	#Calculate probability density at the input x
	mean_mvn = np.zeros(dmnsn)
	x_centered = mvn._process_quantiles(x_centered, dmnsn)
	pdf = mvn(mean_mvn, cov, allow_singular = True).logpdf(x_centered)
	cdf = norm(0, 1).logcdf(np.dot(x_centered, shape))
	return(_squeeze_output(2*np.exp(pdf + cdf)))

#Indirect parametrization with reference to generative process in Azzalini and Capitanio, 1999
def mvsn_pdf2(x, delta, mu = None, Sigma = None, return_log = False):
	#Make numpy arrays storing the inputs
	dmnsn = len(delta)
	delta = np.asarray(delta)
	mu = np.zeros(dmnsn) if mu is None else np.asarray(mu)
	x_centered = x - mu
	Sigma = np.eye(dmnsn) if Sigma is None else np.asarray(Sigma)

	#Calculate probability density at the input x
	mean_mvn = np.zeros(dmnsn)
	x_centered = mvn._process_quantiles(x_centered, dmnsn)
	cov = Sigma + (delta[np.newaxis].T @ delta[np.newaxis])
	cov_inv = np.linalg.inv(cov)
	shape = (1/np.sqrt(1 - delta @ cov_inv @ delta)) * delta @ cov_inv
	pdf = mvn(mean_mvn, cov, allow_singular = True).logpdf(x_centered)
	cdf = norm(0, 1).logcdf(np.dot(x_centered, shape))
	return(_squeeze_output((np.log(2) + pdf + cdf) if return_log else 2*np.exp(pdf + cdf)))

#Direct parametrization
#Sampling based on stochastic process described by Azzalini and Capitanio, 1999
def mvsn_smpl1(n, shape, mu = None, cov = None):
	#Make numpy arrays storing the inputs
	dmnsn = len(shape)
	shape = np.asarray(shape)[np.newaxis].T
	mu = np.zeros(dmnsn) if mu is None else np.asarray(mu)
	cov = np.eye(dmnsn) if cov is None else np.asarray(cov)

	#Set up complementary multivariate normal covariance matrix
	delta = (1/np.sqrt(1 + (shape.T @ cov @ shape))) * (cov @ shape)
	mean_mvn = np.zeros(dmnsn + 1)
	cov_mvn = np.block([	[np.ones(1),	delta.T	], \
				[delta,		cov	]	])

	#Draw n samples
	x = mvn(mean_mvn, cov_mvn, allow_singular = True).rvs(n)
	x0, x1 = x[:, 0], x[:, 1:]
	flpSigns = x0 <= 0
	x1[flpSigns] = -1 * x1[flpSigns]
	return(x1 + mu)

#Indirect parametrization
#Sampling based on stochastic process described by Azzalini and Capitanio, 1999
def mvsn_smpl2(n, delta, mu = None, Sigma = None):
	#Make numpy arrays storing the inputs
	dmnsn = len(delta)
	delta = np.asarray(delta)
	mu = np.zeros(dmnsn) if mu is None else np.asarray(mu)
	Sigma = np.eye(dmnsn) if Sigma is None else np.asarray(Sigma)
	Omega = Sigma + (delta[np.newaxis].T @ delta[np.newaxis])

	#Set up complementary multivariate normal covariance matrix
	mean_mvn = np.zeros(dmnsn + 1)
	cov_mvn = np.block([	[np.ones(1),			delta	], \
				[delta[np.newaxis].T,		Omega	]	])

	#Draw n samples
	x = mvn(mean_mvn, cov_mvn, allow_singular = True).rvs(n)
	x0, x1 = x[:, 0], x[:, 1:]
	flpSigns = x0 <= 0
	x1[flpSigns] = -1 * x1[flpSigns]
	return(x1 + mu)

def fit_mvsnmix_E(y, delta, xi, Sigma, bcdGroup, barcodes, bkgd_p):
	#------------------------
	#E-step
	#------------------------
	#MVSN probability density at each observation
	p = pd.Series([np.nan]*y.shape[0], dtype = 'float64')
	p[y.classificationStatus.isin(barcodes.id.tolist())] = 0
	p[y.classificationStatus == bcdGroup] = 1
	bcdCompatibleStates = barcodes.compatibleStates[barcodes.id == bcdGroup].iloc[0]
	p[y.classificationStatus.isin(bcdCompatibleStates)] = \
		mvsn_pdf2(y.loc[y.classificationStatus.isin(bcdCompatibleStates), \
				barcodes.signalChannels[barcodes.id == bcdGroup].iloc[0]], \
				delta, xi, Sigma) * \
		bkgd_p[y.classificationStatus.isin(bcdCompatibleStates)]
	p[np.isnan(p)] = 0 #All remaining nans are due to incompatible barcodes
	if ((p == 0).all()): #Assume that every barcode is present even if very low counts
		p += (1/p.shape[0]) * 10**(-10)

	Omega = Sigma + (delta[np.newaxis].T @ delta[np.newaxis])
	Omega_inv = np.linalg.inv(Omega)
	delta_matmul_Omega_inv = delta @ Omega_inv
	x_cntrd_scld = (y[barcodes.signalChannels[barcodes.id == bcdGroup].iloc[0]] - \
				xi).dot(delta_matmul_Omega_inv)
	sig = np.sqrt(1 - delta_matmul_Omega_inv @ delta)
	e1Nr = sig*norm(0, 1).pdf(x_cntrd_scld/sig)
	e1Dr = norm(0, 1).cdf(x_cntrd_scld/sig) 
	e1 = x_cntrd_scld + np.divide(e1Nr, e1Dr)
	e1[~np.isfinite(e1)] = e1[np.isfinite(e1)].max()
	e2 = (x_cntrd_scld ** 2) + (sig ** 2) + \
		(x_cntrd_scld * sig * norm(0, 1).pdf(x_cntrd_scld * sig)/norm(0, \
			1).cdf(x_cntrd_scld * sig))
	e2[~np.isfinite(e2)] = e2[np.isfinite(e2)].max()
	return([p, e1, e2])

def fit_mvsn_mix_M(y, xi, r, e1, e2):
	#-----------------------
	#M-step
	#-----------------------
	pi = np.sum(r)/r.shape[0]
	y_cntrd = y - xi
	delta = np.asarray(y_cntrd.apply(lambda x: e1 * r * x).sum()/np.sum(r*e2))	
	xi = np.array((y - pd.DataFrame(np.array(e1)[np.newaxis].T @ \
			delta[np.newaxis])).mul(r, axis = 0).sum()/np.sum(r))
	y_cntrd = np.array(y_cntrd)
	Sigma_term1 = (np.multiply(y_cntrd.reshape(-1, y_cntrd.shape[1], 1), \
			y_cntrd.reshape(-1, 1, y_cntrd.shape[1])) * \
			np.array(r).reshape(-1, 1, 1)).sum(axis = 0)
	Sigma_term2 = (np.multiply(delta.reshape(-1, delta.shape[0], 1), \
			y_cntrd.reshape(-1, 1, y_cntrd.shape[1])) * \
			np.array(r * e1).reshape(-1, 1, 1)).sum(axis = 0) 
	Sigma_term3 = Sigma_term2.T
	Sigma_term4 = np.sum(r * e2) * delta[np.newaxis].T @ delta[np.newaxis]
	Sigma = (1/np.sum(r)) * (Sigma_term1 - Sigma_term2 - Sigma_term3 + Sigma_term4)
	return([pi, delta, xi, Sigma])

#In terms of direct parametrization
def fit_mvsn_mix_logL1(y, xi, delta, Sigma, pi, r, bcdGroup, barcodes, bkgd_p):
	pi = 10**(-30) if (pi == 0) else pi
	Omega = Sigma + (delta[np.newaxis].T @ delta[np.newaxis])
	Omega_inv = np.linalg.inv(Omega)
	alpha = (1/np.sqrt(1 - delta @ Omega_inv @ delta)) * delta @ Omega_inv
	logP = mvsn_pdf2(y[barcodes.signalChannels[barcodes.id == bcdGroup].iloc[0]], \
				delta, xi, Sigma, return_log = True) + \
		np.log(bkgd_p)
	logL1 = np.sum(logP * r)
	logL2 = np.log(pi) * np.sum(r)
	return(logL1 + logL2)

#Indirect parametrization with reference to generative process in Azzalini and Capitanio, 1999
def fit_mvsn_mix_logL2(y, xi, delta, Sigma, e1, e2, pi, r):
	logL1 = - 0.5 * np.log(np.linalg.det(Sigma)) * np.sum(r)
	#U is the MVN distributed component in the generative process, see Pyne et al, 2009
	U = (y - xi) - pd.DataFrame(e1.apply(lambda x: delta*x).tolist())
	Sigma_inv = np.linalg.inv(Sigma)
	logL2 = - 0.5 * np.sum(U.apply(lambda x: np.asarray(x) @ Sigma_inv @ np.asarray(x), \
				axis = 1) * r)
	logL3 = - 0.5 * np.sum(e2 * r)
	logL4 = np.log(pi) * np.sum(r)
	return(logL1 + logL2 + logL3 + logL4)

def classificationStatus(calls, nEpiPerBarcode):
	if np.nansum(calls) > nEpiPerBarcode:
		return(-999) #Represents multiple infections
	else:
		return(np.nansum([calls[x] * 10**(x) for x in range(len(calls))]))

def getEmpiricalProbabilityFunction(inDf, nbins = 100):
	outDfs = {}
	for chnl in inDf.columns:
		hist = np.histogram(inDf[chnl], bins = nbins)
		histDf = pd.DataFrame()
		histDf['binStart'] = np.concatenate(([2 * hist[1][0] - hist[1][1]], hist[1]))
		histDf['binEnd'] = np.concatenate((hist[1], [2 * hist[1][-1] - hist[1][-2]]))
		histDf['binCenter'] = histDf[['binStart', 'binEnd']].mean(axis = 1)
		#Adding small number because finite sample is used to estimate ...
		#... empirical probability function. 
		#... Adding 1/nbins amounts to one observation added to the overall sample.
		#... Adding 1 may be adding too many observations overall.
		histDf['massFraction'] = np.concatenate(([0], hist[0], [0])) + 1 
		histDf['massFraction'] /= histDf.massFraction.sum()
		outDfs[chnl] = histDf
	return(outDfs)	

def bckgdPDF(vals, bkgdDist):
	for chnl in vals.columns:
		if (chnl == vals.columns[0]):
			pbck = vals[chnl].apply(lambda x: bkgdDist[chnl]['massFraction'][ \
				np.abs(bkgdDist[chnl]['binCenter'] - x).idxmin()])
		else:
			pbck = pbck * vals[chnl].apply(lambda x: bkgdDist[chnl]['massFraction'][ \
				np.abs(bkgdDist[chnl]['binCenter'] - x).idxmin()])
	return(pbck)

def loadLastIter(emDir, cut, init, nEpisPerBarcode, nObs, resumeFrom = -1):
	prefix = 'cut_' + str(cut) + '_init_' + str(init) + '_iter'
	files = os.listdir(emDir)
	files = [x for x in files if x.startswith(prefix)]
	if len(files) == 0:
		return([{}, {}, {}, {}, {}, -1])
	startedIters = [x.split('iter_')[1] for x in files]
	startedIters = np.unique([x.split('_')[0] for x in startedIters])
	lastIter = np.max([int(x) for x in startedIters]) if resumeFrom == -1 else resumeFrom

	lastIterDone = False #Assume that last iteration was not properly saved
	while ((not lastIterDone) and (lastIter >= 0)):
		pi = pd.read_csv(os.path.join(emDir, 'cut_' + str(cut) + \
			'_init_' + str(init) + '_iter_' + str(lastIter) + '_pi.txt'), sep = '\t')
		if (pi.shape[0] == 0):
			lastIter -= 1
			continue
		mu = pd.read_csv(os.path.join(emDir, 'cut_' + str(cut) + \
			'_init_' + str(init) + '_iter_' + str(lastIter) + '_mu.txt'), sep = '\t')
		delta = pd.read_csv(os.path.join(emDir, 'cut_' + str(cut) + \
			'_init_' + str(init) + '_iter_' + str(lastIter) + '_delta.txt'), sep = '\t')
		Sigma = pd.read_csv(os.path.join(emDir, 'cut_' + str(cut) + \
			'_init_' + str(init) + '_iter_' + str(lastIter) + '_Sigma.txt'), sep = '\t')
		logL = pd.read_csv(os.path.join(emDir, 'cut_' + str(cut) + \
			'_init_' + str(init) + '_iter_' + str(lastIter) + '_logL.txt'), sep = '\t')
		r = pd.read_csv(os.path.join(emDir, 'cut_' + str(cut) + \
			'_init_' + str(init) + '_iter_' + str(lastIter) + '_r.txt'), sep = '\t')
		if ((mu.shape[0] != nEpisPerBarcode) | (delta.shape[0] != nEpisPerBarcode) | \
			(Sigma.shape[0] != nEpisPerBarcode**2) | (logL.shape[0] != 1) | \
			(r.shape[0] != nObs)):
			lastIter -= 1
			continue
	
		pi = dict(zip([int(x) for x in pi.columns], \
			[pi[x].iloc[0] for x in pi.columns]))
		mu = dict(zip([int(x) for x in mu.columns], \
			[np.array(mu[x]) for x in mu.columns]))
		delta = dict(zip([int(x) for x in delta.columns], \
			[np.array(delta[x]) for x in delta.columns]))
		Sigma = dict(zip([int(x) for x in Sigma.columns], \
			[np.reshape(np.array(Sigma[x]), (nEpisPerBarcode, -1)) for x \
				in Sigma.columns]))
		lastIterDone = True

	if lastIter == -1:
		return([{}, {}, {}, {}, {}, -1])
	totalLogLFile = os.path.join(emDir, 'cut_' + str(cut) + \
		'_init_' + str(init) + '_totalLogL.txt')
	if os.path.isfile(totalLogLFile):
		totalLogL = pd.read_csv(totalLogLFile, sep = '\t')
	totalLogL = pd.DataFrame(np.array([]), columns = ['0'])
	totalLogL = totalLogL['0'].tolist() if (totalLogL.shape[0] == lastIter + 1) else (
			totalLogL['0'].tolist()[:(lastIter + 1)] if \
			(totalLogL.shape[0] >= lastIter + 1) else \
			(totalLogL['0'].tolist() + [np.nan for x in \
				range(lastIter + 1 - totalLogL.shape[0])]))
	return([pi, delta, mu, Sigma, totalLogL, lastIter])







