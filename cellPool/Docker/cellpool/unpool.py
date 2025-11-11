from skimage.io import imread, imsave
import os
import numpy as np
from collections import OrderedDict
from skimage.measure import regionprops_table
import pandas as pd
import sys
from skimage.segmentation import find_boundaries, mark_boundaries
from skimage.exposure import rescale_intensity
from PIL import Image, ImageDraw, ImageFont
Image.MAX_IMAGE_PIXELS = 1000000000 #Turn off DecompressionBombWarning for large mosaics
import imageio
from cellpool.utils import *

def expandBboxAndExtract(inMask, image, dpWidth):
	props = regionprops_table(inMask, intensity_image = image, \
		properties = ('bbox', 'label'))
	props = pd.DataFrame.from_dict(props)
	props['xmin'] = ((props['bbox-1'] + props['bbox-3'])/2 - \
		dpWidth/2).apply(np.floor).astype(int)
	props.loc[props.xmin < 0, 'xmin'] = 0
	props['xmax'] = props['xmin'] + dpWidth
	props.loc[props.xmax > image.shape[1], 'xmax'] = image.shape[1]
	props['ymin'] = ((props['bbox-0'] + props['bbox-2'])/2 - \
		dpWidth/2).apply(np.floor).astype(int)
	props.loc[props.ymin < 0, 'ymin'] = 0
	props['ymax'] = props['ymin'] + dpWidth
	props.loc[props.ymax > image.shape[0], 'ymax'] = image.shape[0]

	#Sometimes the above steps result in sides < dpWidth. Following guards against this.  
	def imgPullAndFrame(xmin, xmax, ymin, ymax, image = image, dpWidth = dpWidth):
		imFrame = np.zeros((dpWidth, dpWidth), dtype = int)
		imFrame[:(ymax - ymin), :(xmax - xmin)] = image[ymin:ymax, xmin:xmax].copy()
		np.nan_to_num(imFrame, copy = False)
		return(imFrame)

	cellSnaps = props[['xmin', 'xmax', 'ymin', 'ymax']].apply( \
		lambda x: imgPullAndFrame(x[0], x[1], x[2], x[3]), \
		axis = 1)
	return({"label" : props.label, "cellSnaps": cellSnaps})


def unpoolObjectsInSameWell(listOfObjIds, dirWithImages, objMaskFilePath, well, objMaskId, \
	keepImIfSuffix = None, keepImIfPrefix = None, keepImIfContains = None, \
	throwImIfSuffix = None, throwImIfPrefix = None, throwImIfContains = None, \
	keepIm = None, keepMaskIfSuffix = None, keepMaskIfPrefix = None, \
	keepMaskIfContains = None, throwMaskIfPrefix = None, throwMaskIfSuffix = None, \
	throwMaskIfContains = None, keepMask = None, dpWidth = 300, \
	labelSize = None, dirWithMasks = None, \
	showOnlyTgtObjs = True, returnAllMasksInDir = False, returnIndObjMask = False): 

	imFiles = os.listdir(dirWithImages) if keepIm is None else keepIm
	imFiles = [x for x in imFiles if x.startswith(well)]
	imFiles = filterList(imFiles, keepImIfPrefix, keepImIfSuffix, keepImIfContains, \
		throwImIfPrefix, throwImIfSuffix, throwImIfContains)

	if (returnAllMasksInDir & (dirWithMasks is None)):
		raise ValueError('If returnAllMasksInDir is True, must input dirWithMasks.')

	if (returnAllMasksInDir):
		maskFiles = os.listdir(dirWithMasks) if keepMask is None else keepMask
		maskFiles = filterList(maskFiles, keepMaskIfPrefix, keepMaskIfSuffix, \
			keepMaskIfContains, throwMaskIfPrefix, throwMaskIfSuffix, \
			throwMaskIfContains)
	objMask = imageio.imread(os.path.join(objMaskFilePath, objMaskId))
	currObjs = [int(x.split('n')[1]) for x in listOfObjIds]
	objMask[~np.isin(objMask, currObjs)] = 0
	
	if returnIndObjMask:
		#objBdries = find_boundaries(objMask, mode = 'inner')
		#objBdries = np.where(objBdries, objMask, 0)
		if showOnlyTgtObjs:
			props = regionprops_table(objMask, intensity_image = objMask, \
				properties = ('image_intensity', 'label'))
			indMasks = OrderedDict(zip([well + 'n' + str(x) for x in props['label']], \
				props['image_intensity']))
			#for idx in props['label']:
			#	bdries[well + 'n' + str(idx)][bdries[well + 'n' \
			#		+ str(idx)] != 0] = 1 
		else:
			indMasks = expandBboxAndExtract(objMask, objMask, dpWidth)
			for idx in range(len(indMasks['label'])):
				indMasks['cellSnaps'][idx][indMasks['cellSnaps'][idx] != \
								indMasks['label'][idx]] = 0 
				indMasks['cellSnaps'][idx][indMasks['cellSnaps'][idx] != 0] = 1
					
			indMasks = OrderedDict(zip([well + 'n' + str(x) for \
				x in indMasks['label']], indMasks['cellSnaps']))

	ims = OrderedDict()
	for nextFile in imFiles:
		oriIm = np.load(os.path.join(dirWithImages, nextFile))
		if showOnlyTgtObjs:
			props = regionprops_table(objMask, intensity_image = oriIm, \
				properties = ('image_intensity', 'label'))
			ims[nextFile[4:]] = OrderedDict(zip([well + 'n' + str(x) for \
				x in props['label']], props['image_intensity']))
		else:
			ims[nextFile[4:]] = expandBboxAndExtract(objMask, oriIm, dpWidth)
			ims[nextFile[4:]] = OrderedDict(zip([well + 'n' + str(x) for \
				x in ims[nextFile[4:]]['label']], ims[nextFile[4:]]['cellSnaps']))

	if (returnAllMasksInDir):
		masks = OrderedDict()
		for nextFile in maskFiles:
			oriMask = imageio.imread(os.path.join(dirWithMasks, nextFile))
			masks[nextFile] = expandBboxAndExtract(objMask, oriMask, dpWidth)
			masks[nextFile] = OrderedDict(zip([well + 'n' + str(x) for \
				x in masks[nextFile]['label']], masks[nextFile]['cellSnaps']))
	
	outDict = {"images": ims} 
	if returnAllMasksInDir: 
		outDict["masks"] = masks
	if returnIndObjMask:
		outDict["mainMasks"] = indMasks

	return(outDict)
	
#Function that reads in obj labels from a file, splits them by well and passes on to unpool by well 
def unpoolObjects(fileWithListOfObjIds, dirWithImages, objMaskFilePath, objMaskId = None, \
	keepImIfSuffix = None, keepImIfPrefix = None, keepImIfContains = None, \
	throwImIfSuffix = None, throwImIfPrefix = None, throwImIfContains = None, \
	keepIm = None, keepMaskIfSuffix = None, keepMaskIfPrefix = None, \
	keepMaskIfContains = None, throwMaskIfPrefix = None, throwMaskIfSuffix = None, \
	throwMaskIfContains = None, keepMask = None, dpWidth = 300, \
	labelSize = None, dirWithMasks = None, \
	showOnlyTgtObjs = True, returnAllMasksInDir = False, returnIndObjMask = False, \
	keepImIfWellPrefix = False, objMaskIdWellPrefix = False, objMaskIdContains = None, \
	objMaskIdFindExt = False): 
	
	allObjs = pd.read_csv(fileWithListOfObjIds, sep = '\t', header = None)[[0]]
	wells = np.unique([x.split('n')[0] for x in allObjs[0].tolist()])
	wellByWellIms = OrderedDict()
	for well in wells:
		thisWellObjs = [x for x in allObjs[0].tolist() if x.startswith(well)]
		keepImIfPrefix = well if keepImIfWellPrefix else keepImIfPrefix
		objMaskIdPrefix = well + objMaskIdContains if objMaskIdWellPrefix else None
		objMaskIdSuffix = os.path.splitext( \
			[x for x in os.listdir(objMaskFilePath) if \
				x.startswith(objMaskIdPrefix)][0])[1] if \
					objMaskIdFindExt else None
		thisObjMaskId = objMaskIdPrefix + objMaskIdSuffix if \
				objMaskId is None else objMaskId
		wellByWellIms[well] = unpoolObjectsInSameWell(thisWellObjs, \
			dirWithImages, objMaskFilePath, well, thisObjMaskId, \
			keepImIfSuffix, keepImIfPrefix, keepImIfContains, \
			throwImIfSuffix, throwImIfPrefix, throwImIfContains, \
			keepIm, keepMaskIfSuffix, keepMaskIfPrefix, \
			keepMaskIfContains, throwMaskIfPrefix, throwMaskIfSuffix, \
			throwMaskIfContains, keepMask, dpWidth, \
			labelSize, dirWithMasks, \
			showOnlyTgtObjs, returnAllMasksInDir, returnIndObjMask) 
	consolidated = OrderedDict()
	consolidated["images"] = OrderedDict()
	imKeys = [x for x in wellByWellIms[wells[0]]['images']]
	for key in imKeys:
		consolidated["images"][key] = OrderedDict()
		for nextObjId in allObjs[0].tolist():
			consolidated["images"][key][nextObjId] = \
				wellByWellIms[nextObjId[:3]]['images'][key][nextObjId]

	if (returnAllMasksInDir):
		consolidated["masks"] = OrderedDict()
		maskKeys = [x for x in wellByWellIms[wells[0]]['masks']]
		for key in maskKeys:
			consolidated["masks"][key] = OrderedDict()
			for nextObjId in allObjs[0].tolist():
				consolidated["masks"][key][nextObjId] = \
					wellByWellIms[nextObjId[:3]]['masks'][key][nextObjId]

	if returnIndObjMask:	
		consolidated["mainMasks"] = OrderedDict()
		for nextObjId in allObjs[0].tolist():
			consolidated["mainMasks"][nextObjId] = \
				wellByWellIms[nextObjId[:3]]['mainMasks'][nextObjId]
	return(consolidated)


#Function receives dict of images and boundaries and plots the image matrix
def showCellImageMatrix(inmages, maskImgs = None, dpWidth = 300, trimFromStainId = None, \
	markBoundaries = False, showObjIds = False, fontsize = 20, vadjust = 10, \
	hadjust = 10, \
	font = '/app/UbuntuMono-R.ttf'):
	if (markBoundaries & (maskImgs is None)):
		raise ValueError('If markBoundaries is True, maskImgs cannot be None.')

	stains = [x for x in inmages]
	objs = [x for x in inmages[stains[0]]]
	nStains = len(stains)
	nObjs = len(inmages[stains[0]])
		
	if trimFromStainId is None:
		stainSuffix = dict(zip(stains, stains))
	else:
		suffix = ['_' + x.replace(trimFromStainId, '') for x in stains]
		stainSuffix = dict(zip(stains, suffix))
	
	outmage_chs = 3 if markBoundaries else 1
	outmage = np.zeros((nObjs * dpWidth, nStains * dpWidth, outmage_chs))
	for obj in range(nObjs):
		for stn in range(nStains):
			thisOriShape = inmages[stains[stn]][objs[obj]].shape
			thisImage = inmages[stains[stn]][objs[obj]].copy()
			if markBoundaries:
				thisImage = mark_boundaries(thisImage, \
					maskImgs[objs[obj]], color = (1, 0, 0), mode = 'thick')
			outmage[((dpWidth * obj) + int((dpWidth - thisOriShape[0])/2)): \
				((dpWidth * obj) + int((dpWidth - thisOriShape[0])/2) + \
				thisOriShape[0]), ((dpWidth * stn) + \
				int((dpWidth - thisOriShape[1])/2)):((dpWidth * stn) + \
				int((dpWidth - thisOriShape[1])/2) + thisOriShape[1]), :] = \
					thisImage #inmages[stains[stn]][objs[obj]]

	outmage[[x*dpWidth for x in range(nObjs)], :, :] = 1
	outmage[:, [x*dpWidth for x in range(nStains)], :] = 1
	outmage = (outmage*255).astype(np.uint8)

	#if showBoundaries:
	#	bdrimage = showCellImageMatrix(inmages = {'key' : boundaries}, dpWidth = dpWidth)
	#	bdrimage = np.tile(bdrimage, nStains)
	#	outmage = np.dstack((outmage, bdrimage, bdrimage))		

	if showObjIds:
		outmage = Image.fromarray(outmage)
		write = ImageDraw.Draw(outmage)
		myFont = ImageFont.truetype(font, size = fontsize)
		for obj in range(nObjs):
			for stn in range(nStains):
				write.text((dpWidth * stn + hadjust, dpWidth * obj + vadjust), \
					#objs[obj] + stainSuffix[stains[stn]], \
					"", font=myFont, \
					fill = (255, 255, 255) if markBoundaries else (255,))

	return(outmage)	

#Function takes in dict with single image channel per object and single object mask
#Makes all images and masks to dpWidth using np.pad
#Lays them out in a table with specified row and column layout
#Saves them in a specified format
def showUnpooled(gray, red = None, green = None, maskImgs = None, \
	blue = None, dpWidth = 300, ncols = 50, markBoundaries = False, \
	showObjIds = False, fontsize = 20, vadjust = 10, hadjust = 10, \
	font = '/app/UbuntuMono-R.ttf'):
	if ((gray is None) & (red is None) & (green is None) & (blue is None)):
		raise ValueError('No input images.')

	if ((gray is None) & ((red is None) | (green is None) | (blue is None))):
		raise ValueError('Either only gray or all three of RGB channels must be input.')

	if ((gray is not None) & ((red is not None) | (green is not None) | (blue is not None))):
		raise ValueError('Too many input image channels, not sure which ones to plot.')

	rgb = red is not None
	dpWidth = int(dpWidth)
	ncols = int(ncols)
	nrows = int(np.ceil(len(gray)/ncols)) if gray is not None else int(np.ceil(len(red)/ncols))
	outImage = np.zeros((int(nrows * dpWidth), int(ncols * dpWidth), \
		3 if rgb else 1))
	maskImage = np.zeros((int(nrows * dpWidth), int(ncols * dpWidth)), dtype = int)
	objLabels = [x for x in red.keys()] if rgb else [x for x in gray.keys()]
	for ix in range(len(objLabels)):
		thisOriShape = red[objLabels[ix]].shape if rgb else gray[objLabels[ix]].shape
		outImage[(dpWidth * int(np.floor(ix/ncols)) + int((dpWidth - \
			thisOriShape[0])/2)):(dpWidth * int(np.floor(ix/ncols)) + \
			int((dpWidth - thisOriShape[0])/2) + thisOriShape[0]), \
			(dpWidth * int(ix % ncols) + \
			int((dpWidth - thisOriShape[1])/2)):(dpWidth * int(ix % ncols) + \
			int((dpWidth - thisOriShape[1])/2) + \
			thisOriShape[1]), 0] = \
			red[objLabels[ix]] if rgb else gray[objLabels[ix]]
		if markBoundaries:
			maskImage[(dpWidth * int(np.floor(ix/ncols)) + int((dpWidth - \
				thisOriShape[0])/2)):(dpWidth * int(np.floor(ix/ncols)) + \
				int((dpWidth - thisOriShape[0])/2) + thisOriShape[0]), \
				(dpWidth * int(ix % ncols) + \
				int((dpWidth - thisOriShape[1])/2)):(dpWidth * int(ix % ncols) + \
				int((dpWidth - thisOriShape[1])/2) + \
				thisOriShape[1])] = \
				maskImgs[objLabels[ix]]
		if rgb:
			outImage[(dpWidth * int(np.floor(ix/ncols)) + int((dpWidth - \
				thisOriShape[0])/2)):(dpWidth * int(np.floor(ix/ncols)) + \
				int((dpWidth - thisOriShape[0])/2) + thisOriShape[0]), \
				(dpWidth * int(ix % ncols) + \
				int((dpWidth - thisOriShape[1])/2)):(dpWidth * int(ix % ncols) + \
				int((dpWidth - thisOriShape[1])/2) + \
				thisOriShape[1]), 1] = green[objLabels[ix]]

			outImage[(dpWidth * int(np.floor(ix/ncols)) + int((dpWidth - \
				thisOriShape[0])/2)):(dpWidth * int(np.floor(ix/ncols)) + \
				int((dpWidth - thisOriShape[0])/2) + thisOriShape[0]), \
				(dpWidth * int(ix % ncols) + \
				int((dpWidth - thisOriShape[1])/2)):(dpWidth * int(ix % ncols) + \
				int((dpWidth - thisOriShape[1])/2) + \
				thisOriShape[1]), 2] = blue[objLabels[ix]]

	outImage[[x*dpWidth for x in range(nrows)], :, :] = 1
	outImage[:, [x*dpWidth for x in range(ncols)], :] = 1
	outImage = (outImage*255).astype(np.uint8)
	
	if markBoundaries:
		outImage = (mark_boundaries(outImage, \
			maskImage.astype(np.int64), color = (1, 1, 1))*255).astype(np.uint8)

	if showObjIds:
		outImage = Image.fromarray(outImage if rgb else outImage[:, :, 0])
		write = ImageDraw.Draw(outImage)
		myFont = ImageFont.truetype(font, size = fontsize)
		for ix in range(len(objLabels)):
			write.text((dpWidth * int(ix % ncols) + hadjust, \
				dpWidth * int(np.floor(ix/ncols)) + vadjust), \
				objLabels[ix], font=myFont, \
				fill = (255, 255, 255) if rgb else (255,))
	return(outImage)

#Add method based on winsorization
def normalizeForViewing(images, method = "max_scaling", limits = [2, 98]):
	for x in images:
		if (method == "max_scaling"):
			images[x] = images[x]/np.nanmax(images[x])
		elif (method == "rescale_to_range"):
			pctls = np.percentile(images[x][images[x] > 0], limits)
			images[x] = rescale_intensity(images[x], in_range = tuple(pctls), \
						out_range='float32')
		else:
			raise ValueError(method + " is not an available normalization method.")
	return(images)

#Add method to save unpooled images as pngs
#All pngs output will be identical and single cell in all possible formats as in showUnpooled but saved individually
def saveUnpooled(gray, outDir = '.', red = None, green = None, \
	blue = None, dpWidth = 300, ncols = 50, markBoundaries = False, maskImgs = None, \
	showObjIds = False, fontsize = 20, vadjust = 10, hadjust = 10, extension = ".png", \
	font = '/app/UbuntuMono-R.ttf'):
	if ((gray is None) & (red is None) & (green is None) & (blue is None)):
		raise ValueError('No input images.')

	if ((gray is None) & ((red is None) | (green is None) | (blue is None))):
		raise ValueError('Either only gray or all three of RGB channels must be input.')

	if ((gray is not None) & ((red is not None) | (green is not None) | (blue is not None))):
		raise ValueError('Too many input image channels, not sure which ones to plot.')

	if ((maskImgs is None) & markBoundaries):
		raise ValueError('Requested mark boundaries but object masks not provided.')

	rgb = red is not None
	dpWidth = int(dpWidth)
	
	if showObjIds:
		myFont = ImageFont.truetype(font, size = fontsize)

	allImages = OrderedDict()
	objLabels = [x for x in red.keys()] if rgb else [x for x in gray.keys()]
	for ix in range(len(objLabels)):
		thisOriShape = red[objLabels[ix]].shape if rgb else gray[objLabels[ix]].shape
		outImage = np.zeros((dpWidth, dpWidth, 3 if rgb else 1))
		outImage[(int((dpWidth - thisOriShape[0])/2)):(int((dpWidth - \
			thisOriShape[0])/2) + thisOriShape[0]), \
			(int((dpWidth - thisOriShape[1])/2)):(int((dpWidth - \
			thisOriShape[1])/2) + thisOriShape[1]), 0] = \
			red[objLabels[ix]].copy() if rgb else gray[objLabels[ix]].copy()
		if rgb:
			outImage[(int((dpWidth - thisOriShape[0])/2)):(int((dpWidth - \
				thisOriShape[0])/2) + thisOriShape[0]), \
				(int((dpWidth - thisOriShape[1])/2)):(int((dpWidth - \
				thisOriShape[1])/2) + thisOriShape[1]), 1] = \
					green[objLabels[ix]].copy()

			outImage[(int((dpWidth - thisOriShape[0])/2)):(int((dpWidth - \
				thisOriShape[0])/2) + thisOriShape[0]), \
				(int((dpWidth - thisOriShape[1])/2)):(int((dpWidth - \
				thisOriShape[1])/2) + thisOriShape[1]), 2] = \
					blue[objLabels[ix]].copy()

		if markBoundaries:
			outImage = (mark_boundaries(outImage, \
				maskImgs[objLabels[ix]], color = (1, 1, 1))*255).astype(np.uint8)

		if showObjIds:
			outImage = Image.fromarray(outImage if rgb else outImage[:, :, 0])
			write = ImageDraw.Draw(outImage)
			write.text((hadjust, vadjust), \
				objLabels[ix], font=myFont, \
				fill = (255, 255, 255) if rgb else (255,))
			outImage.save(os.path.join(outDir, objLabels[ix] + extension))
		else:
			imsave(os.path.join(outDir, objLabels[ix] + extension), outImage)

	return()









	
