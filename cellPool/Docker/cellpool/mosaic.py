import pandas as pd
import networkx as nx
import sys
import numpy as np
from skimage.registration import phase_cross_correlation
from scipy import ndimage
from scipy.fft import fft2
import itertools

#---------------------------------------------------------------
#---------------------------------------------------------------
def addNghbrInfo(fieldsData, tileGph):
	"""
	fieldsData: metadata about imaging tiles in each well 
	tileGph: networkx object with graph of tile connections.
	"""
	dirOpposites = {'West': 'East', 'East': 'West', 'North': 'South', 'South': 'North'}  
	for nghbrs in tileGph.edges:
		im1Id = int(fieldsData.FieldID[fieldsData.FieldID == nghbrs[0]].iloc[0])
		im2Id = int(fieldsData.FieldID[fieldsData.FieldID == nghbrs[1]].iloc[0])
		if (fieldsData[fieldsData.FieldID == nghbrs[1]].PositionX.iloc[0] == \
			fieldsData[fieldsData.FieldID == nghbrs[0]].PositionX.iloc[0]):
			if (fieldsData[fieldsData.FieldID == nghbrs[1]].PositionY.iloc[0] < \
				fieldsData[fieldsData.FieldID == nghbrs[0]].PositionY.iloc[0]):
				im2LocWrtIm1 = 'South'
			else:
				im2LocWrtIm1 = 'North'
		else: #if PositionX is not equal then, PositionY must be. 
			if (fieldsData[fieldsData.FieldID == nghbrs[1]].PositionX.iloc[0] < \
				fieldsData[fieldsData.FieldID == nghbrs[0]].PositionX.iloc[0]):
				im2LocWrtIm1 = 'West'
			else:
				im2LocWrtIm1 = 'East'
		fieldsData.loc[fieldsData.FieldID == nghbrs[1], \
			'ngbrTo' + dirOpposites[im2LocWrtIm1]] = im1Id
		fieldsData.loc[fieldsData.FieldID == nghbrs[0], \
			'ngbrTo' + im2LocWrtIm1] = im2Id
	
	return(fieldsData)

#---------------------------------------------------------------
#---------------------------------------------------------------
def buildLayeredTileGraph(tgtCycleGph, refTileCoords): 
	"""
	Currently, assuming that target and reference cycles have identical tiling. If not, 
	could use the metadata for the two cycles and allowed maxSeparation to build the 
	layered graph.
	"""

	mapping = dict(zip(tgtCycleGph.nodes, ['tgt_' + str(x) for x in tgtCycleGph.nodes]))
	layeredTileGph = nx.relabel_nodes(tgtCycleGph, mapping, copy=True)
	refTiles = [int(x) for x in refTileCoords[ \
		~refTileCoords['medX'].isnull()].fieldID.str.replace('f', '')]
	layeredTileGph.add_edge('ref_1', 'tgt_1')
	edgesToAdd = [('ref_1', 'ref_' + str(x)) for x in refTiles if x != 1]
	for x in refTiles:
		edgesToAdd.append(('ref_' + str(x), 'tgt_' + str(x)))

	layeredTileGph.add_edges_from(edgesToAdd)

	return(layeredTileGph)


#---------------------------------------------------------------
#Receives a metadata table as read from Harmony exported XML.
#Assumes a 96-well plate.
#Outputs a graph linking neighbors.
#---------------------------------------------------------------
def buildNeighborGraphOfTiles(fieldsData, maxSeparation = 1.2): 
	"""
	maxSeparation: maximum separation (in units of tile size) to still consider neighbors.
	"""

	#Block distance of 1 means adjacent block
	fieldsData['BlockDistanceFromField1X'] = fieldsData['PositionX'] / \
		(fieldsData['ImageResolutionX'] * fieldsData['ImageSizeX'])
	fieldsData['BlockDistanceFromField1Y'] = fieldsData['PositionY'] / \
		(fieldsData['ImageResolutionY'] * fieldsData['ImageSizeY'])

	tileGraph = nx.Graph()
	tileGraph.add_nodes_from([x for x in fieldsData.FieldID])
	for i in range(int(fieldsData.FieldID.min()), int(fieldsData.FieldID.max())):
		for j in range(i+1, int(fieldsData.FieldID.max()) + 1):
			thisFields = fieldsData[fieldsData.FieldID.isin([i, j])]
			thisBlockDistsFromField1X = thisFields.BlockDistanceFromField1X.tolist()
			thisBlockDistsX = abs(thisBlockDistsFromField1X[0] - \
				thisBlockDistsFromField1X[1])
			thisBlockDistsFromField1Y = thisFields.BlockDistanceFromField1Y.tolist()
			thisBlockDistsY = abs(thisBlockDistsFromField1Y[0] - \
				thisBlockDistsFromField1Y[1])
			if (thisBlockDistsX + thisBlockDistsY < maxSeparation):
				tileGraph.add_edge(i, j)

	return(tileGraph)


#---------------------------------------------------------------
#---------------------------------------------------------------
def buildTileGraphForRegistration(tgtTileCoords, refTileCoords): 
	"""
	Currently, assuming that target and reference cycles have identical tiling. If not, 
	could use the metadata for the two cycles and allowed maxSeparation to build the 
	layered graph.
	"""
	G = nx.Graph()
	refTiles = ['ref_F' + f'{int(x):02d}' for x in \
			refTileCoords[~refTileCoords['medX'].isnull()].fieldID]
	tgtTiles = ['tgt_F' + f'{int(x):02d}' for x in \
			tgtTileCoords[~tgtTileCoords['medX'].isnull()].fieldID]
	edgesToAdd = [('ref_F01', x) for x in refTiles if x != 'ref_F01'] + \
		[('tgt_F01', x) for x in tgtTiles if x != 'tgt_F01'] + \
		[(x, 'tgt_' + x[4:]) for x in refTiles]
	G.add_edges_from(edgesToAdd)
	return(G)


#---------------------------------------------------------------
#---------------------------------------------------------------
def cropAndEstimateCoordsWrtNghbr(images, coordEsts, cropFraction = 0.25):
	"""
	images: directory with fields (aka tiles) as keys and image as a numpy array
	coordEsts: output from filtering coordinate estimates 
	"""
	directions = {'West': 'East', 'East': 'West', 'North': 'South', 'South': 'North'}
	imageSize = coordEsts.ImageSizeX.unique()[0] #Assuming square tiles
	cropSize = round(imageSize * cropFraction)
	adjustForInfTilingByPhsCorr = {'West' : {'X': imageSize, 'Y': 0}, \
		'East': {'X': -imageSize, 'Y': 0}, 'North': {'X': 0, 'Y': imageSize}, \
		'South': {'X': 0, 'Y': -imageSize}} 
	for row in range(coordEsts.shape[0]):
		#If the field in this row has already been placed properly, skip to next iteration.
		if (not np.isnan(coordEsts.iloc[row][['locEstYWrtNghbrTo' + x for \
			x in directions]].tolist()).all()):
			continue
	
		im1Id = int(coordEsts.iloc[row].FieldID)
		for nghbrDir in directions:
			thisNghbr = coordEsts.iloc[row]['ngbrTo' + nghbrDir]
			if (np.isnan(thisNghbr)):
				continue
			
			im2Id = int(coordEsts[coordEsts.FieldID == thisNghbr].FieldID.iloc[0])
			if (nghbrDir == 'North'):
				im1 = images['F' + f"{im1Id:02d}"][:cropSize, :]
				im2 = images['F' + f"{im2Id:02d}"][-cropSize:, :]
			elif (nghbrDir == 'South'):
				im1 = images['F' + f"{im1Id:02d}"][-cropSize:, :]
				im2 = images['F' + f"{im2Id:02d}"][:cropSize, :]
			elif (nghbrDir == 'West'):
				im1 = images['F' + f"{im1Id:02d}"][:, :cropSize]
				im2 = images['F' + f"{im2Id:02d}"][:, -cropSize:]
			elif (nghbrDir == 'East'):
				im1 = images['F' + f"{im1Id:02d}"][:, -cropSize:]
				im2 = images['F' + f"{im2Id:02d}"][:, :cropSize]
			
			thisOffset = register(im1, im2)[0]
		
			#Register function returns how to translate im2 for best alignment with im1.
			coordEsts.loc[coordEsts.FieldID == im2Id, \
				'locEstYWrtNghbrTo' + directions[nghbrDir]] = thisOffset[0] + \
				adjustForInfTilingByPhsCorr[directions[nghbrDir]]['Y']
			coordEsts.loc[coordEsts.FieldID == im2Id, \
				'locEstXWrtNghbrTo' + directions[nghbrDir]] = thisOffset[1] + \
				adjustForInfTilingByPhsCorr[directions[nghbrDir]]['X']
			coordEsts.loc[coordEsts.FieldID == im1Id, \
				'locEstYWrtNghbrTo' + nghbrDir] = -thisOffset[0] + \
				adjustForInfTilingByPhsCorr[nghbrDir]['Y']
			coordEsts.loc[coordEsts.FieldID == im1Id, \
				'locEstXWrtNghbrTo' + nghbrDir] = -thisOffset[1] + \
				adjustForInfTilingByPhsCorr[nghbrDir]['X']
	
	return(coordEsts)

#---------------------------------------------------------------
#---------------------------------------------------------------
def estimateTileCoordsWrtImRoot(coordEsts, tileGph, root = 1, maxDepth = 15, maxPaths = 10):
	centroidLocs = {}
	directions = ['North', 'South', 'East', 'West']
	for toPlace in coordEsts.FieldID:
		thisPaths = []
		currDepth = 1 if (maxDepth > 15) else maxDepth
		if nx.has_path(tileGph, root, toPlace):
			while (len(thisPaths) < maxPaths) & (currDepth <= maxDepth):
				allSimplePaths = list(nx.all_simple_paths(tileGph, root, \
					toPlace, cutoff = currDepth))
				allSimplePaths = sorted(allSimplePaths, key = len)
				thisPaths.extend(allSimplePaths[:maxPaths])
				thisPaths = [list(x) for x in set(tuple(x) for x in thisPaths)]
				currDepth += 1

		thisPaths = sorted(thisPaths, key = len)[:maxPaths]
		centroidLocs[toPlace] = []
		for path in thisPaths:
			thisX, thisY = [0, 0]
			previous_node = path[0]
			for thisNode in path[1:]:
				thisCoordEsts = coordEsts[coordEsts.FieldID == thisNode]
				thisNodeNghbrs = dict(zip(directions, \
					[thisCoordEsts['ngbrTo' + x].iloc[0] for x in directions]))
				prevNodeDirWrtThisNode = [x for x in thisNodeNghbrs.keys() if \
					thisNodeNghbrs[x] == previous_node][0]
				thisX += thisCoordEsts['locEstXWrtNghbrTo' + \
					prevNodeDirWrtThisNode].iloc[0]
				thisY += thisCoordEsts['locEstYWrtNghbrTo' + \
					prevNodeDirWrtThisNode].iloc[0]
				previous_node = thisNode
			centroidLocs[toPlace].append([path, thisX, thisY])
	return(centroidLocs)

#---------------------------------------------------------------
#---------------------------------------------------------------
def estimateMstTileCoordsWrtImRoot(coordEsts, tileGph, root = 1):
	centroidLocs = {}
	directions = ['North', 'South', 'East', 'West']
	for toPlace in coordEsts.FieldID:
		path = nx.shortest_path(tileGph, root, toPlace)
		centroidLocs[toPlace] = []
		thisX, thisY = [0, 0]
		previous_node = path[0]
		for thisNode in path[1:]:
			thisCoordEsts = coordEsts[coordEsts.FieldID == thisNode]
			thisNodeNghbrs = dict(zip(directions, \
				[thisCoordEsts['ngbrTo' + x].iloc[0] for x in directions]))
			prevNodeDirWrtThisNode = [x for x in thisNodeNghbrs.keys() if \
				thisNodeNghbrs[x] == previous_node][0]
			thisX += thisCoordEsts['locEstXWrtNghbrTo' + \
				prevNodeDirWrtThisNode].iloc[0]
			thisY += thisCoordEsts['locEstYWrtNghbrTo' + \
				prevNodeDirWrtThisNode].iloc[0]
			previous_node = thisNode
		centroidLocs[toPlace].append([path, thisX, thisY])

	return(centroidLocs)

#---------------------------------------------------------------
#---------------------------------------------------------------
def estimateTileCoordsWrtNghbr(images, fieldsData, tileGph):
	"""
	images: directory with fields (aka tiles) as keys and image as a numpy array
	fieldsData: metadata about the image tiles in a well 
	tileGph: networkx object with graph of tile connections
	"""
	dirOpposites = {'West': 'East', 'East': 'West', 'North': 'South', 'South': 'North'}  
	for nghbrs in tileGph.edges:
		im1Id = fieldsData[fieldsData.FieldID == nghbrs[0]].FieldID.iloc[0]
		im2Id = fieldsData[fieldsData.FieldID == nghbrs[1]].FieldID.iloc[0]
		if (fieldsData[fieldsData.FieldID == nghbrs[1]].PositionX.iloc[0] == \
			fieldsData[fieldsData.FieldID == nghbrs[0]].PositionX.iloc[0]):
			if (fieldsData[fieldsData.FieldID == nghbrs[1]].PositionY.iloc[0] < \
				fieldsData[fieldsData.FieldID == nghbrs[0]].PositionY.iloc[0]):
				im2LocWrtIm1 = 'South'
			else:
				im2LocWrtIm1 = 'North'
		else: #if PositionX is not equal then, PositionY must be. 
			if (fieldsData[fieldsData.FieldID == nghbrs[1]].PositionX.iloc[0] < \
				fieldsData[fieldsData.FieldID == nghbrs[0]].PositionX.iloc[0]):
				im2LocWrtIm1 = 'West'
			else:
				im2LocWrtIm1 = 'East'
		thisOffset = register(images['F' + f"{im1Id:02d}"], \
					images['F' + f"{im2Id:02d}"])[0]
		#Register function returns how to translate im2 so it is best aligned with im1.
		fieldsData.loc[fieldsData.FieldID == nghbrs[1], \
			'locEstYWrtNghbrTo' + dirOpposites[im2LocWrtIm1]] = thisOffset[0]
		fieldsData.loc[fieldsData.FieldID == nghbrs[1], \
			'locEstXWrtNghbrTo' + dirOpposites[im2LocWrtIm1]] = thisOffset[1]
		fieldsData.loc[fieldsData.FieldID == nghbrs[1], \
			'ngbrTo' + dirOpposites[im2LocWrtIm1]] = im1Id
		fieldsData.loc[fieldsData.FieldID == nghbrs[0], \
			'locEstYWrtNghbrTo' + im2LocWrtIm1] = -thisOffset[0]
		fieldsData.loc[fieldsData.FieldID == nghbrs[0], \
			'locEstXWrtNghbrTo' + im2LocWrtIm1] = -thisOffset[1]
		fieldsData.loc[fieldsData.FieldID == nghbrs[0], \
			'ngbrTo' + im2LocWrtIm1] = im2Id
	
	imageSizeX = fieldsData.ImageSizeX.unique()[0]
	imageSizeY = fieldsData.ImageSizeY.unique()[0]

	#Register function returns the minimum offset assuming infinitely repeated tiling of im1.
	fieldsData.locEstXWrtNghbrToWest += imageSizeX
	fieldsData.locEstXWrtNghbrToEast -= imageSizeX
	fieldsData.locEstYWrtNghbrToNorth += imageSizeY
	fieldsData.locEstYWrtNghbrToSouth -= imageSizeY
	return(fieldsData)

#---------------------------------------------------------------
#---------------------------------------------------------------
def estimateTileCoordsWrtRefCycleRoot(coordEsts, refTileCoords, layeredTileGph, \
					root = 'ref_1', maxDepth = 15, maxPaths = 10):
	centroidLocs = {}
	directions = ['ToNorth', 'ToSouth', 'ToEast', 'ToWest', 'InRefCycle']
	for toPlace in coordEsts.FieldID:
		thisPaths = []
		currDepth = 0
		while (len(thisPaths) < maxPaths) & (currDepth < maxDepth):
			thisPaths.extend(list(nx.all_simple_paths(layeredTileGph, root, \
				'tgt_' + str(toPlace), cutoff = currDepth))[:maxPaths])
			thisPaths = [list(x) for x in set(tuple(x) for x in thisPaths)]
			currDepth += 1
		thisPaths = sorted(thisPaths, key = len)[:maxPaths]
		centroidLocs[toPlace] = []
		for path in thisPaths:
			thisX, thisY = [0, 0]
			previous_node = path[0]
			thisFieldNum = int(previous_node[4:])
			thisRefCoords = refTileCoords[refTileCoords.fieldID == \
				'f' + f"{thisFieldNum:02}"]
			thisX += thisRefCoords.medX.iloc[0]
			thisY += thisRefCoords.medY.iloc[0]
			for thisNode in path[1:]:
				if (thisNode.startswith('ref_')):
					thisFieldNum = int(thisNode[4:])
					thisRefCoords = refTileCoords[refTileCoords.fieldID == \
						'f' + f"{thisFieldNum:02}"]
					thisX += thisRefCoords.medX.iloc[0]
					thisY += thisRefCoords.medY.iloc[0]
					previous_node = thisNode
					continue
				thisCoordEsts = coordEsts[coordEsts.FieldID == int(thisNode[4:])]
				thisNodeNghbrs = dict(zip(directions, \
					[thisCoordEsts['ngbr' + x].iloc[0] for x in directions]))
				prevNodeDirWrtThisNode = [x for x in thisNodeNghbrs.keys() if \
					thisNodeNghbrs[x] == previous_node][0]
				thisX += thisCoordEsts['locEstXWrtNghbr' + \
					prevNodeDirWrtThisNode].iloc[0]
				thisY += thisCoordEsts['locEstYWrtNghbr' + \
					prevNodeDirWrtThisNode].iloc[0]
				previous_node = int(thisNode[4:])
			centroidLocs[toPlace].append([path, thisX, thisY])
	return(centroidLocs)

#---------------------------------------------------------------
#---------------------------------------------------------------
def filterUnrealisticEstimates(coordEsts, maxDevFromMedTranslation = 3, filterWrtRefCycle = False, \
				maxTranslation = 20, expectedDispOfNgbrs = 'NA'):
	"""
		coordEsts: dataframe with estimated locations of tiles based on phase correlation.
	"""
	#Remove neighbors with excessive estimated translation.
	#Assumes that there are minor mechanical errors from a calibrated camera movement.
	if (filterWrtRefCycle):
		coordEsts.loc[(abs(coordEsts.locEstXWrtNghbrInRefCycle - \
			np.nanmedian(coordEsts.locEstXWrtNghbrInRefCycle)) > \
			maxDevFromMedTranslation) | (abs(coordEsts.locEstYWrtNghbrInRefCycle - \
			np.nanmedian(coordEsts.locEstYWrtNghbrInRefCycle)) > \
			maxDevFromMedTranslation), ['locEstXWrtNghbrInRefCycle', \
			'locEstYWrtNghbrInRefCycle']] = np.nan
		coordEsts.loc[(abs(coordEsts.locEstXWrtNghbrInRefCycle) > maxTranslation) | \
			(abs(coordEsts.locEstYWrtNghbrInRefCycle) > maxTranslation), \
			['locEstXWrtNghbrInRefCycle', 'locEstYWrtNghbrInRefCycle']] = np.nan
		return(coordEsts)

	if expectedDispOfNgbrs == 'NA':
		coordEsts.loc[(abs(coordEsts.locEstXWrtNghbrToWest - \
			np.nanmedian(coordEsts.locEstXWrtNghbrToWest)) > \
			maxDevFromMedTranslation) | \
			(abs(coordEsts.locEstYWrtNghbrToWest) > maxDevFromMedTranslation), \
			['locEstXWrtNghbrToWest', 'locEstYWrtNghbrToWest']] = np.nan
		coordEsts.loc[(abs(coordEsts.locEstXWrtNghbrToEast - \
			np.nanmedian(coordEsts.locEstXWrtNghbrToEast)) > \
			maxDevFromMedTranslation) | \
			(abs(coordEsts.locEstYWrtNghbrToEast) > maxDevFromMedTranslation), \
			['locEstXWrtNghbrToEast', 'locEstYWrtNghbrToEast']] = np.nan
		coordEsts.loc[(abs(coordEsts.locEstYWrtNghbrToNorth - \
			np.nanmedian(coordEsts.locEstYWrtNghbrToNorth)) > \
			maxDevFromMedTranslation) | \
			(abs(coordEsts.locEstXWrtNghbrToNorth) > maxDevFromMedTranslation), \
			['locEstXWrtNghbrToNorth', 'locEstYWrtNghbrToNorth']] = np.nan
		coordEsts.loc[(abs(coordEsts.locEstYWrtNghbrToSouth - \
			np.nanmedian(coordEsts.locEstYWrtNghbrToSouth)) > \
			maxDevFromMedTranslation) | \
			(abs(coordEsts.locEstXWrtNghbrToSouth) > maxDevFromMedTranslation), \
			['locEstXWrtNghbrToSouth', 'locEstYWrtNghbrToSouth']] = np.nan
	else:
		coordEsts.loc[(abs(abs(coordEsts.locEstXWrtNghbrToWest) - expectedDispOfNgbrs) > \
			maxDevFromMedTranslation) | \
			(abs(coordEsts.locEstYWrtNghbrToWest) > maxDevFromMedTranslation), \
			['locEstXWrtNghbrToWest', 'locEstYWrtNghbrToWest']] = np.nan
		coordEsts.loc[(abs(abs(coordEsts.locEstXWrtNghbrToEast) - expectedDispOfNgbrs) > \
			maxDevFromMedTranslation) | \
			(abs(coordEsts.locEstYWrtNghbrToEast) > maxDevFromMedTranslation), \
			['locEstXWrtNghbrToEast', 'locEstYWrtNghbrToEast']] = np.nan
		coordEsts.loc[(abs(abs(coordEsts.locEstYWrtNghbrToNorth) - expectedDispOfNgbrs) > \
			maxDevFromMedTranslation) | \
			(abs(coordEsts.locEstXWrtNghbrToNorth) > maxDevFromMedTranslation), \
			['locEstXWrtNghbrToNorth', 'locEstYWrtNghbrToNorth']] = np.nan
		coordEsts.loc[(abs(abs(coordEsts.locEstYWrtNghbrToSouth) - expectedDispOfNgbrs) > \
			maxDevFromMedTranslation) | \
			(abs(coordEsts.locEstXWrtNghbrToSouth) > maxDevFromMedTranslation), \
			['locEstXWrtNghbrToSouth', 'locEstYWrtNghbrToSouth']] = np.nan

	return(coordEsts)

#---------------------------------------------------------------
#---------------------------------------------------------------
def pruneTileGph(tileGph, coordEsts):
	"""
		coordEsts: dataframe with estimated locations of tiles based on phase correlation.
	"""
	
	edges_to_remove = []
	for row in range(coordEsts.shape[0]):
		thisField = coordEsts.FieldID.iloc[row]
		badNghbr = [coordEsts['ngbrTo' + x].iloc[row] for x in \
			['North', 'South', 'East', 'West'] if \
			np.isnan(coordEsts['locEstXWrtNghbrTo' + x].iloc[row])]
		badNghbr = [x for x in badNghbr if not np.isnan(x)]
		edges_to_remove.append([(thisField, int(x)) for x in badNghbr])
	
	#Flatten edges_to_remove which is a list of lists.
	edges_to_remove = [edge for sublist in edges_to_remove for edge in sublist]

	tileGph.remove_edges_from(edges_to_remove)
	return(tileGph)

#---------------------------------------------------------------
#---------------------------------------------------------------
def register(im1, im2, sigma = 3, upsample_factor = 10, normalization = None):
	im1_LoG = ndimage.gaussian_laplace(im1, sigma)
	im2_LoG = ndimage.gaussian_laplace(im2, sigma)
	
	im1_fft = fft2(im1_LoG)
	im2_fft = fft2(im2_LoG)
	
	#Returns displacement that will best align im2 wrt im1.
	#1st entry in shift is the y offset with positive values meaning im2 should be shifted down.
	#2nd entry in shift is the x offset with positive values meaning im2 should be shifted right.
	shift, error, diffphase = phase_cross_correlation(im1_fft, im2_fft, \
					upsample_factor = upsample_factor, space = 'fourier', \
					normalization = normalization)
	
	return([shift, error, diffphase])

#---------------------------------------------------------------
#---------------------------------------------------------------
def registerWrtRefCycle(tgtImages, refImages, layeredTileGph):
	#Assumes that each tile in tgtImages has a single neighbor in refImages
	regCoords = {}
	for thisField in tgtImages:
		thisFieldNode = 'tgt_' + thisField
		if (thisFieldNode not in layeredTileGph.nodes):
			continue
		nghbr = [x for x in layeredTileGph.neighbors(thisFieldNode) \
			if x.startswith('ref_')]
		if (len(nghbr) == 1):
			nghbr = nghbr[0]
		else:
			continue
		nghbrFieldId = nghbr[4:]
		thisOffset = register(refImages[nghbrFieldId], \
					tgtImages[thisField])[0]
		regCoords[thisField] = {'tgtField': thisField, \
				'refField': nghbrFieldId, \
				'locEstYWrtNghbrInRefCycle': thisOffset[0], \
				'locEstXWrtNghbrInRefCycle': thisOffset[1]}

	regCoords = pd.DataFrame.from_dict(regCoords).T
	regCoords.reset_index(drop = True, inplace = True)
	return(regCoords)

#---------------------------------------------------------------
#---------------------------------------------------------------
def registerWrtRefCycleWithSimultaneousStitchAndReg(tgtImages, refImages, \
		coordEsts, layeredTileGph):
	#Assumes that each tile in tgtImages has a single neighbor in refImages
	for thisField in tgtImages:
		thisFieldNode = 'tgt_' + str(int(thisField[1:]))
		nghbr = [x for x in layeredTileGph.neighbors(thisFieldNode) \
			if x.startswith('ref_')][0]
		nghbrFieldId = int(nghbr[4:])
		thisOffset = register(refImages['f' + f"{nghbrFieldId:02}"], \
					tgtImages[thisField])[0]
		coordEsts.loc[coordEsts.FieldID == int(thisField[1:]), \
			'ngbrInRefCycle'] = nghbr
		coordEsts.loc[coordEsts.FieldID == int(thisField[1:]), \
			'locEstYWrtNghbrInRefCycle'] = thisOffset[0]
		coordEsts.loc[coordEsts.FieldID == int(thisField[1:]), \
			'locEstXWrtNghbrInRefCycle'] = thisOffset[1]
	return(coordEsts)

#---------------------------------------------------------------
#---------------------------------------------------------------
def rgstrTgtCycleRootWrtRefRoot(regCoords, refStchCoords, tgtStchCoords, \
		root, tileGph):
	allPaths = list(nx.all_simple_paths(tileGph, 'ref_' + root, 'tgt_' + root))
	locEstsX = []
	locEstsY = []
	pathCenterFld = []
	for path in allPaths:
		if len(path) == 2:
			locEstsY.append(regCoords[regCoords.tgtField == \
				root]['locEstYWrtNghbrInRefCycle'].iloc[0])
			locEstsX.append(regCoords[regCoords.tgtField == \
				root]['locEstXWrtNghbrInRefCycle'].iloc[0])
			pathCenterFld.append(root)
			continue

		fldBtwnRoots = path[1][4:]
		thisX = refStchCoords[refStchCoords.fieldID == int(fldBtwnRoots[1:])]['medX'].iloc[0] + \
			regCoords[regCoords.tgtField == fldBtwnRoots][ \
				'locEstXWrtNghbrInRefCycle'].iloc[0] - \
			tgtStchCoords[tgtStchCoords.fieldID == int(fldBtwnRoots[1:])]['medX'].iloc[0]
		thisY = refStchCoords[refStchCoords.fieldID == int(fldBtwnRoots[1:])]['medY'].iloc[0] + \
			regCoords[regCoords.tgtField == fldBtwnRoots][ \
				'locEstYWrtNghbrInRefCycle'].iloc[0] - \
			tgtStchCoords[tgtStchCoords.fieldID == int(fldBtwnRoots[1:])]['medY'].iloc[0]
		locEstsX.append(thisX)
		locEstsY.append(thisY)
		pathCenterFld.append(fldBtwnRoots)

	locEsts = pd.DataFrame(pathCenterFld, columns = ['pathCenterField'])
	locEsts['Y'] = locEstsY
	locEsts['X'] = locEstsX
	return(locEsts)

#---------------------------------------------------------------
#---------------------------------------------------------------
def smrzCoordEsts(centroidLocs, coordEsts, root): 
	smryCntrdLocs = {}
	for thisField in coordEsts.FieldID:
		thisFieldImId = int(coordEsts.FieldID[coordEsts['FieldID'] == thisField].iloc[0])
		if (thisField == root):
			smryCntrdLocs[thisField] = {'fieldID' : thisFieldImId, \
				'medX' : 0, 'medY' : 0, 'stdX' : 0, 'stdY' : 0, 'nPaths' : 0}
			continue
		if (len(centroidLocs[thisField]) == 0):
			smryCntrdLocs[thisField] = {'fieldID' : thisFieldImId, \
					'medX' : np.nan, 'medY' : np.nan, \
					'stdX' : np.nan, 'stdY' : np.nan, 'nPaths' : 0}
			continue
		thisMedianX = np.nanmedian([round(x[1], 1) for x in centroidLocs[thisField]])
		thisMedianY = np.nanmedian([round(x[2], 1) for x in centroidLocs[thisField]])
		if (len(centroidLocs[thisField]) == 1):
			thisStdX, thisStdY = [np.nan, np.nan]
		else:
			thisStdX = np.nanstd([round(x[1], 1) for x in centroidLocs[thisField]])
			thisStdY = np.nanstd([round(x[2], 1) for x in centroidLocs[thisField]])
		smryCntrdLocs[thisField] = {'fieldID' : thisFieldImId, \
					'medX' : round(thisMedianX, 1), \
					'medY' : round(thisMedianY, 1), \
					'stdX' : round(thisStdX, 1), \
					'stdY' : round(thisStdY, 1), \
					'nPaths' : len(centroidLocs[thisField])}
	
	smryCntrdLocs = pd.DataFrame.from_dict(smryCntrdLocs).T
	return(smryCntrdLocs)






