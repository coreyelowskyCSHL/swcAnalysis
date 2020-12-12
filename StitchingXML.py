from sys import exit
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt


"""
Stitching XML 
	- provides functions to parse Stitching XML file from Big Stitcher
"""


class StitchingXML:

	def __init__(self, xmlPath, numYVolumes=None, numZVolumes=None):
		

		self.xmlPath = xmlPath
		self.root =  ET.parse(xmlPath).getroot()
	
		# parse XML
		self.setupsAndRegistrations = self.parseSetups()
		self.parseRegistrations()
		self.pairwiseStitchings = self.parsePairwiseStitchings()

		if numYVolumes == None or numZVolumes == None:
			self.numYVolumes, self.numZVolumes = self.inferYZVolumes()
		else:
			self.numYVolumes = numYVolumes
			self.numZVolumes = numZVolumes

		self.numVolumes = self.numYVolumes*self.numZVolumes
		self.anisotropyFactor = self.calculateAnisotropyFactor()



	def calculateAnisotropyFactor(self):
		
		voxelSize = self.setupsAndRegistrations[0]['voxel size']
		return voxelSize[2]/voxelSize[1]


	def setupIDtoY(self, setupID):
		return (setupID % self.numYVolumes) + 1

	def setupIDtoZ(self, setupID):
		return (setupID // self.numYVolumes) + 1

	def setupIDtoVolume(self, setupID):

		y = '{:02}'.format(self.setupIDtoY(setupID))
		z = '{:02}'.format(self.setupIDtoZ(setupID))

		return 'Z' + z + '_Y' + y

	def volumeToSetupID(self, volume):

		z = int(volume.split('_')[0][1:])
		y = int(volume.split('_')[1][1:])

		return (z-1)*self.numYVolumes + (y-1)

	def parseSetups(self):

		"""
		Parses setup information into dictionary where setup ids are keys
		"""

		setups = {}

		for setup in self.root.iter('ViewSetup'):

			setupID = int(setup.find('id').text)
			size = setup.find('size').text.split(' ')
			voxelSize = setup.find('voxelSize').find('size').text.split(' ')
			
			setups[setupID] = {'size':[int(x) for x in size], 'voxel size':[float(x) for x in voxelSize]}

		return setups

	def parseRegistrations(self):

		"""	
		Parses registration information which are three matrices
		and stores them in setupsAndRegistrations dictionary 
		"""
		
		for registration in self.root.iter('ViewRegistration'):

			setupID = int(registration.attrib['setup'])

			for transform in registration.iter('ViewTransform'):

				name = transform.find('Name').text
				matrix = np.array([float(x) for x in transform.find('affine').text.split()]).reshape(3,4)
				self.setupsAndRegistrations[setupID][name] = np.vstack((matrix,np.array([0,0,0,1])))

	def parsePairwiseStitchings(self):
		
		"""
		Parses pairwise stitchings and stores in list of dictionaries
		"""


		pairwiseStitchings = []

		for pairwise in self.root.iter('PairwiseResult'):

			setupID_a = int(pairwise.attrib['view_setup_a'])
			setupID_b = int(pairwise.attrib['view_setup_b'])
			shift = np.array([float(x) for x in pairwise.find('shift').text.split()]).reshape(3,4)
			correlation = np.round(float(pairwise.find('correlation').text),2)
			bbox = [float(x) for x in pairwise.find('overlap_boundingbox').text.split()]

			pairwiseStitchings.append({'setupID_a':setupID_a,'setupID_b':setupID_b,'shift':shift,'correlation':correlation,'bbox':bbox})
			

		return pairwiseStitchings

	
	def inferYZVolumes(self):
		
		"""
		Calculates number of Y and Z volumes based on translsation to grid matrices
		"""

		x1 = self.setupsAndRegistrations[0]['Translation to Regular Grid'][0,-1]

		for i,setupID in enumerate(sorted(self.setupsAndRegistrations.keys())):
			if i > 0:			
				x2 = self.setupsAndRegistrations[setupID]['Translation to Regular Grid'][0,-1]

				if x1 == x2:
					numYVolumes = i
					break

		numZVolumes = int(len(self.setupsAndRegistrations)/numYVolumes)

		if numZVolumes*numYVolumes != len(self.setupsAndRegistrations):
			exit('Error: Number of y and z volumes was not infered correctly!')	

		return numYVolumes, numZVolumes


	def analyzeCorrelations(self, * , zPair = None, yPair = None):
		
		print()
		print('Analyzing Correlations...')
		print()

		yCorrs, zCorrs, yzCorrs = [], [], []

		for pairwise in self.pairwiseStitchings:

			setupID_a = pairwise['setupID_a']
			setupID_b = pairwise['setupID_b']

			if zPair:
				if type(zPair) != set:
					exit('Error: zPair must be of type set!')
				if {self.setupIDtoZ(setupID_a), self.setupIDtoZ(setupID_b)} != zPair:
					continue

			if yPair:
				if type(yPair) != set:
					exit('Error: zPair must be of type set!')
				if {self.setupIDtoY(setupID_a), self.setupIDtoY(setupID_b)} != zPair:
					continue
			
			correlation = pairwise['correlation']

			diff = abs(setupID_a - setupID_b)
	
			if self.setupIDtoZ(setupID_a) == self.setupIDtoZ(setupID_b):
				yCorrs.append(correlation)
			elif self.setupIDtoY(setupID_a) == self.setupIDtoY(setupID_b):
				zCorrs.append(correlation)
			else:
				yzCorrs.append(correlation)

		yCorrs.sort()
		zCorrs.sort()
		yzCorrs.sort()

		yCorrs = yCorrs[::-1]
		zCorrs = zCorrs[::-1]
		yzCorrs = yzCorrs[::-1]
		
		print('Y Corrs:',yCorrs)
		print('Z Corrs:',zCorrs)
		print('YZ Corrs:',yzCorrs)


	def getStitchingMatrices(self, volume):

		"""
		Gets transformation matrices for given volume

		Params: volume, e.g. Z12_Y06
		Return: translation to grid matrix, stitching matrix, calibration matrix stored in dictionary 

		"""

		stitchingDict = self.setupsAndRegistrations[self.volumeToSetupID(volume)]

		stitchingMatrices = {}
		stitchingMatrices['Translation to Regular Grid'] = stitchingDict['Translation to Regular Grid']
		stitchingMatrices['Stitching Transform'] = stitchingDict['Stitching Transform']
		stitchingMatrices['calibration'] = stitchingDict['calibration']

		return stitchingMatrices

	def getFusedDimensions(self, * , downsampling=1):
		
		"""
		Calculated dimensions of big stitcher image after fusing

		Params: volumeSize, pixel dimensions of volumes [xLength, yLength, zLength]
		Return: 
		"""

		xCoords, yCoords, zCoords = [], [], []

		for setupID in self.setupsAndRegistrations:

			volumeSize = self.setupsAndRegistrations[setupID]['size']

			volume = self.setupIDtoVolume(setupID)
			stitchingMatrices = self.getStitchingMatrices(volume)
			translationToGridMatrix = stitchingMatrices['Translation to Regular Grid']
			stitchingMatrix = stitchingMatrices['Stitching Transform']

			# append min and max bounds
			xCoords.append(translationToGridMatrix[0,3] + stitchingMatrix[0,3])
			xCoords.append(translationToGridMatrix[0,3] + stitchingMatrix[0,3] + volumeSize[0])
			yCoords.append(translationToGridMatrix[1,3] + stitchingMatrix[1,3])
			yCoords.append(translationToGridMatrix[1,3] + stitchingMatrix[1,3] + volumeSize[1])
			zCoords.append((translationToGridMatrix[2,3] + stitchingMatrix[2,3])/self.anisotropyFactor)
			zCoords.append((translationToGridMatrix[2,3] + stitchingMatrix[2,3])/self.anisotropyFactor + volumeSize[2])

		minX, maxX = min(xCoords)/downsampling, max(xCoords)/downsampling
		minY, maxY = min(yCoords)/downsampling, max(yCoords)/downsampling
		minZ, maxZ = min(zCoords)/downsampling, max(zCoords)/downsampling
		lengthX = int(np.round(maxX - minX))
		lengthY = int(np.round(maxY - minY))
		lengthZ = int(np.round(maxZ - minZ))

		dimensions = {'anisotropy factor':self.anisotropyFactor ,'x':{'min':minX,'max':maxX,'length':lengthX},'y':{'min':minY,'max':maxY,'length':lengthY},'z':{'min':minZ,'max':maxZ,'length':lengthZ}}
		
		return dimensions


					


if __name__ == '__main__':

	xml = StitchingXML('/data/elowsky/OLST_2/data/Vglut_stitching/N5_gzip_36Z/dataset.xml~2')

	xml.analyzeCorrelations(zPair={3,4})



	


