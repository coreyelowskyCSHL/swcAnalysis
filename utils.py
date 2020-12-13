import numpy as np
import copy
import params
import xml.etree.ElementTree as ET

from os.path import join
from os import system
from Node import Node
from sys import exit


def buildNode(swcRow, * ,  parentNode=None, preserveSID=True, preserveRadius=True):

	"""
	Params: swcRow - 1D array, [id, structureID, x, y, z, radius, parentID]
		parentNode - Node
		preserveSID - boolean
		preserveRadius - boolean

	Return: Node object populated with values from swcRow

	"""
	
	sID = swcRow[params.SWC_INDICES['sid']] if preserveSID else 0
	radius = swcRow[params.SWC_INDICES['radius']] if preserveRadius else 0
	
	node = Node(swcRow[params.SWC_INDICES['id']],sID,swcRow[params.SWC_INDICES['x']],swcRow[params.SWC_INDICES['y']],swcRow[params.SWC_INDICES['z']],radius,parentNode=parentNode)

	return node

def fixSWCArrayTypes(swcArray):

	"""
	Params: swcArray, 2D array, swc array
	
	Changes data type of id,structure id,radius,parent id, to int
	Chnages data type of x,y,z to float

	"""

	swcArray[:,[params.SWC_INDICES['id'],params.SWC_INDICES['sid'],params.SWC_INDICES['radius'],params.SWC_INDICES['pid']]] = swcArray[:,[params.SWC_INDICES['id'],params.SWC_INDICES['sid'],params.SWC_INDICES['radius'],params.SWC_INDICES['pid']]].astype(float).astype(int)
	swcArray[:,[params.SWC_INDICES['x'],params.SWC_INDICES['y'],params.SWC_INDICES['z']]] = swcArray[:,[params.SWC_INDICES['x'],params.SWC_INDICES['y'],params.SWC_INDICES['z']]].astype(float)
	
	return swcArray

def getSomaRowFromSWCArray(swcArray):
	
	somaRow = swcArray[swcArray[:,params.SWC_INDICES['pid']]==params.SOMA_PID]

	if len(somaRow) != 1:
		exit('Error: # of somas is ' + str(len(somaRow)) +' and must be 1!')

	return swcArray[swcArray[:,params.SWC_INDICES['pid']]==params.SOMA_PID][0]


def getChildRowsFromSWCArray(swcArray, parentID):

	return swcArray[swcArray[:,params.SWC_INDICES['pid']] == parentID]




###### Utility functions for Oblique -> Coronal Transformation #######

# These functions assume transformations are occuring on images
# therefore, the size of the image must be provided and the
# input coordinates are assumed to be in pixel coordinates of the oblique image

def obliqueToCoronal(coords, fusedDims, shearFactor):


	# convert to numpy arrays
	coords = np.array(coords)
	imageLengths = np.array([fusedDims['x']['length'],fusedDims['y']['length'],fusedDims['z']['length']])

	# shift coordinates 
	coords -= np.array([fusedDims['x']['min'],fusedDims['y']['min'],fusedDims['z']['min']])		

	# perform transformations
	coords, imageLengths = reslice(coords,imageLengths)
	coords = verticalFlip(coords, imageLengths)
	coords, imageLengths = shear(coords,imageLengths,shearFactor)
	coords, imageLengths = reslice(coords,imageLengths)
	coords, imageLengths = rotate(coords,imageLengths)
	
	return list(coords)
		

def reslice(coords, imageLengths):

	return coords[[1,2,0]], imageLengths[[1,2,0]]


def verticalFlip(coords, imageLengths):
	
	coords[1] = imageLengths[1] - coords[1] - 1

	return coords


def shear(coords, imageLengths, shearFactor):

	imageLengthsShear = copy.deepcopy(imageLengths)
	imageLengthsShear[1] = imageLengths[1] + abs(shearFactor)*imageLengthsShear[0]
	imageLengthsShear = np.round(imageLengthsShear).astype(int)

	coords[1] = np.round(coords[1] - imageLengths[1]/2 + shearFactor*(coords[0]-imageLengths[0]/2) + imageLengthsShear[1]/2 )

	return coords, imageLengthsShear

def rotate(coords,imageLengths):

	coords = coords[[1,0,2]]
	coords[0] = imageLengths[1] - coords[0] - 1
	
	return coords, imageLengths[[1,0,2]]

#################################################################

######### Utility Functions for Branch Markup #####################

def parseBranchMarkupFile(path):

	coords, ids = [], [1]

	with open(path,'r') as fp:
		lines = [x for x in fp.readlines() if x != '\n']
		for i,line in enumerate(lines):
			if i == 0:
				coord = [int(x) for x in line.split(': ')[1][:-1].split(',')]
			else:
				coord = [int(x) for x in line.split(': ')[-1][:-1].split(' ')]
			coords.append(coord)
			if i > 0:
				ids.append(line.split(' :')[:2])

	return np.array(coords), ids

def writeBranchMarkupFile(path, coords, ids):

	with open(path, 'w') as fp:
		for i,coord in enumerate(coords):
			if i == 0:
				fp.write('Soma: ' + str(coord[0]) + ',' + str(coord[1]) + ',' + str(coord[2]) + '\n\n')
			else:
				fp.write(ids[i][0] + ' : ' + ids[i][1] +' : ' + str(coord[0]) + ',' + str(coord[1]) + ',' + str(coord[2]) + '\n')

def readBranchMarkups(path):

	labelInfo = []

	with open(path,'r') as fp:
		lines = fp.readlines()
		for line in lines:
			if '\t' in line:
				id_ = int(line.split(':')[1])
				label = int([l for l in line.split('\t') if l != '' and l != '\n'][-1][0])
				labelInfo.append((id_,label))
	
	return labelInfo

#############################################################


def getSWCBrainAndID(swcName):

	# remove extension
	if '.' in swcName:
		swcName = swcName.split('.')[0]

	# special case for 170329_500
	if swcName.startswith('170329'):
		brain = '170329_500'
		swcID = int(swcName.split('_')[-1])
	else:
		brain = swcName.split('_')[0]
		swcID = int(swcName.split('_')[-1])

	brainSWCID = brain + '_' + str(swcID)

	return brain, swcID, brainSWCID

def loadCropInfo():

	cropInfo = np.genfromtxt(params.CROP_INFO_PATH,dtype=object)
	cropInfo[:,0] = cropInfo[:,0].astype(str)
	cropInfo[:,1] = cropInfo[:,1].astype(int)

	cropDict = {}

	for cropping in cropInfo:
		cropDict[cropping[0]] = cropping[1]

	return cropDict


def writeCoordsForRegistration(path, coords):

	with open(path,'w') as fp:
		fp.write('index\n')
		fp.write(str(len(coords))+'\n')
		for i,coord in enumerate(coords):
			fp.write(str(coord[0])+' ')
			fp.write(str(coord[1])+' ')
			fp.write(str(coord[2]))
			if i != len(coords)-1:
				fp.write('\n')

def parseTransformixOutput(path):

	registeredCoords = np.zeros(shape=(0,3))
	with open(path,'r') as fp:
		lines = fp.readlines()
		for line in lines:
			coords = np.array([x for x in line.split(';') if 'OutputPoint' in x][0].split(' ')[4:7],dtype=np.float)
			registeredCoords = np.vstack((registeredCoords,coords))

	return registeredCoords

def runTransformix(outPath, tpPath, inPath):
	
	system('transformix -out ' + outPath  + ' -tp ' + tpPath + ' -def ' + inPath)




