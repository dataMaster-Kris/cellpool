import pandas as pd
import os
import numpy as np
from datetime import datetime
import re
import pandas as pd
import string
import xml.etree.ElementTree as ET
import pandas as pd
import os
import numpy as np

def filterList(l, keepIfPrefix = None, keepIfSuffix = None, keepIfContains = None, \
	throwIfPrefix = None, throwIfSuffix = None, throwIfContains = None):

	if (keepIfPrefix is not None):
		l = [x for x in l if x.startswith(keepIfPrefix)]

	if (keepIfSuffix is not None):
		l = [x for x in l if x.endswith(keepIfSuffix)]

	if (keepIfContains is not None):
		l = [x for x in l if keepIfContains in x]

	if (throwIfPrefix is not None):
		l = [x for x in l if not x.startswith(throwIfPrefix)]

	if (throwIfSuffix is not None):
		l = [x for x in l if not x.endswith(throwIfSuffix)]

	if (throwIfContains is not None):
		l = [x for x in l if not throwIfContains in x]

	return(l)

def getCycleDirMatch(rawDataDir, delims = None, idFormat = None, \
	keepIfPrefix = [], keepIfSuffix = [], keepIfContains = [], \
	throwIfPrefix = [], throwIfSuffix = [], throwIfContains = []):
	"""
	keep and throw arguments must of type list. Example: ["Tom", "U2OS"]
	"""
	dirs = os.listdir(rawDataDir)
	
	for x in keepIfPrefix:
		dirs = filterList(dirs, keepIfPrefix = x)
	for x in keepIfSuffix:
		dirs = filterList(dirs, keepIfSuffix = x)
	for x in keepIfContains:
		dirs = filterList(dirs, keepIfContains = x)
	for x in throwIfPrefix:
		dirs = filterList(dirs, throwIfPrefix = x)
	for x in throwIfSuffix:
		dirs = filterList(dirs, throwIfSuffix = x)
	for x in throwIfContains:
		dirs = filterList(dirs, throwIfContains = x)
	
	cycleInformation = pd.DataFrame(dirs, columns = ['dirs'])	

	#Extract the string with date and time
	dates = [x.split(delims[0])[1].split(delims[1])[0] for x in cycleInformation['dirs']]
	dirAcquisitionOrder = np.argsort([datetime.strptime(x, idFormat) for x in dates])
	
	cycleInformation = cycleInformation.reindex(dirAcquisitionOrder)
	cycleInformation.reset_index(drop = True, inplace = True)
	cycleInformation['cycle'] = np.arange(1, 1 + cycleInformation.shape[0]) 
	return(cycleInformation)


def getImageFileInfo(filenames, pattern):
	r = re.compile(pattern)
	outDat = {}
	for ix in range(len(filenames)):
		outDat[ix] = [m.groupdict() for m in r.finditer(filenames[ix])][0]
		outDat[ix]['filename'] = filenames[ix]
	outDat = pd.DataFrame.from_dict(outDat).T
	outDat[['Row', 'Column', 'Field', 'Z', 'Channel']] = \
		outDat[['Row', 'Column', 'Field', 'Z', 'Channel']].astype('int32', copy = False)
	outDat['Well'] = outDat[['Row', 'Column']].apply( \
		lambda x: list(string.ascii_uppercase)[x[0] - 1] + \
				f"{x[1]:02d}", axis = 1)
	return outDat

#The function works for images exported by the Perkin Elmer's Harmony software.
#Assumes file format according to http://www.perkinelmer.com/PEHH/HarmonyV5
def getChannelInfo(FfcProfileXmlFileDir):
	FfcProfileXmlFilePath = os.path.join(FfcProfileXmlFileDir, \
					os.listdir(FfcProfileXmlFileDir)[0])
	tree = ET.parse(FfcProfileXmlFilePath)
	root = tree.getroot()
	
	#Get the XML namespace
	xmlns = root.tag.split('}')[0].strip('{')
	
	#Subelement of root with flat-field profiles per channel
	#Findall returns a list with all subelements named map. Hence, need to extract index 0.
	mapping = root.findall('{' + xmlns + '}Map')[0] 
	
	entries = [child[0].text for child in mapping]
	
	channelInfo = []
	for entry in entries:
		channel = entry.split('Channel: ')[1].split(',')[0]
		channelName = entry.split('ChannelName: ')[1].split(',')[0]
		dimX, dimY = [int(x) for x in \
			entry.split('Dims: ')[1].split('],')[0].strip('[').split(', ')]
		channelInfo.append({'channel': channel, 'channelName': channelName, \
			'dimX': dimX, 'dimY': dimY})
	
	channelInfo = pd.DataFrame.from_dict(channelInfo)
	return(channelInfo)


#The function works for images exported by the Perkin Elmer's Harmony software.
#Assumes file format according to http://www.perkinelmer.com/PEHH/HarmonyV5
def getMetadata(imageIndexXmlFilePath):
	tree = ET.parse(imageIndexXmlFilePath)
	root = tree.getroot()
	
	#Get the XML namespace
	xmlns = root.tag.split('}')[0].strip('{')
	
	#Subelement of root with image descriptions per channel
	#Findall returns a list with all subelements. Hence, need to extract index 0.
	images = root.findall('{' + xmlns + '}Images')[0] 
	
	metadataTbl = []
	for child in images:
		thisImageTags = []
		thisImageTexts = []
		for gChild in child:
			thisImageTags.append(gChild.tag.split('}')[1])
			thisImageTexts.append(gChild.text)
		metadataTbl.append(dict(zip(thisImageTags, thisImageTexts)))

	metadataTbl = pd.DataFrame.from_dict(metadataTbl)
	numericCols = ['Row', 'Col', 'FieldID', 'PlaneID', 'TimepointID', 'ChannelID', \
		'FlimID', 'ImageResolutionX', 'ImageResolutionY', 'ImageSizeX', \
		'ImageSizeY', 'BinningX', 'BinningY', 'MaxIntensity', 'PositionX', \
		'PositionY', 'PositionZ', 'AbsPositionZ', 'MeasurementTimeOffset', \
		'MainExcitationWavelength', 'MainEmissionWavelength', 'ObjectiveMagnification', \
		'ObjectiveNA', 'ExposureTime']
	metadataTbl[numericCols] = metadataTbl[numericCols].apply(pd.to_numeric, errors='coerce')
	return(metadataTbl)


def distribSummarize(regionmask, intensity):
	this = intensity[regionmask]
	out = dict(zip(['q' + str(x*10) for x in range(11)], \
			np.percentile(this, q=(0, 10, 20, 30, 40, \
							50, 60, 70, 80, 90, 100))))
	out['sum'] = np.sum(this)
	out['mean'] = np.mean(this)
	out['std'] = np.std(this)
	out['median'] = np.median(this)
	return out



