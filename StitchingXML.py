import xml.etree.ElementTree as ET
import numpy as np
from sys import exit
import matplotlib.pyplot as plt


"""
Stitching XML 
	- provides functions to parse Stitching XML file from Big Stitcher
"""


class StitchingXML:

	def __init__(self, xmlPath, yVolumes=None, zVolumes=None):
		
		self.xmlPath = xmlPath
		self.root =  ET.parse(xmlPath).getroot()

		# parse XML
		self.stitchings = self.parseSetups()
		self.parseRegistrations()
		self.pairwiseStitchings = self.parsePairwiseStitchings()

		if yVolumes == None and zVolumes == None:
			print()
			print('Infer # stacks from XML...')
			self.yVolumes, self.zVolumes = self.inferYZVolumes()
			print('# Y Volumes:',self.yVolumes)
			print('# Z Volumes:',self.zVolumes)
	
		else:
			self.yVolumes = yVolumes
			self.zVolumes = zVolumes

		self.numVolumes = self.yVolumes*self.zVolumes



	def parseSetups(self):

		setups = {}

		for setup in self.root.iter('ViewSetup'):

			id_ = int(setup.find('id').text)
			size = setup.find('size').text.split(' ')
			voxelSize = setup.find('voxelSize').find('size').text.split(' ')
			
			setups[id_] = {'size':size, 'voxel size':voxelSize}

		return setups

	def parseRegistrations(self):
		
		for registration in self.root.iter('ViewRegistration'):

			id_ = int(registration.attrib['setup'])

			for transform in registration.iter('ViewTransform'):

				name = transform.find('Name').text
				self.stitchings[id_][name] = np.array([float(x) for x in transform.find('affine').text.split()]).reshape(3,4)

	def parsePairwiseStitchings(self):
		

		pairwiseStitchings = []

		for pairwise in self.root.iter('PairwiseResult'):

			id_a = int(pairwise.attrib['view_setup_a'])
			id_b = int(pairwise.attrib['view_setup_b'])
			shift = np.array([float(x) for x in pairwise.find('shift').text.split()]).reshape(3,4)
			correlation = np.round(float(pairwise.find('correlation').text),2)
			bbox = [float(x) for x in pairwise.find('overlap_boundingbox').text.split()]

			pairwiseStitchings.append({'id_a':id_a,'id_b':id_b,'shift':shift,'correlation':correlation,'bbox':bbox})
			

		return pairwiseStitchings

	
	def inferYZVolumes(self):
		
		x1 = self.stitchings[0]['Translation to Regular Grid'][0,-1]

		
		for i,id_ in enumerate(sorted(self.stitchings.keys())):
			if i > 0:			
				x2 = self.stitchings[id_]['Translation to Regular Grid'][0,-1]

				if x1 == x2:
					yVolumes = i
					break

		zVolumes = int(len(self.stitchings)/yVolumes)

		if zVolumes*yVolumes != len(self.stitchings):
			exit('Error: Number of y and z volumes was not infered correctly!')	

		return yVolumes, zVolumes


	### allow to ask for specific y or z
	def analyzeCorrelations(self, zPair = None, yPair = None):
		
		print()
		print('Analyzing Correlations...')
		print()

		yCorrs, zCorrs, yzCorrs = [], [], []

		for pairwise in self.pairwiseStitchings:

			id_a = pairwise['id_a']
			id_b = pairwise['id_b']

			if zPair:
				if type(zPair) != set:
					exit('Error: zPair must be of type set!')
				if {self.getZ(id_a), self.getZ(id_b)} != zPair:
					continue

			if yPair:
				if type(yPair) != set:
					exit('Error: zPair must be of type set!')
				if {self.getY(id_a), self.getY(id_b)} != zPair:
					continue
			

			correlation = pairwise['correlation']

			diff = abs(id_a - id_b)
	
			if self.getZ(id_a) == self.getZ(id_b):
				yCorrs.append(correlation)
			elif self.getY(id_a) == self.getY(id_b):
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


	def getY(self, id_):
		return (id_ % self.yVolumes) + 1

	def getZ(self, id_):
		return (id_ // self.yVolumes) + 1

	def idToVolumeID(self, id_):
		y = '{:02}'.format(self.getY(id_))
		z = '{:02}'.format(self.getZ(id_))

		return 'Z' + z + '_Y' + y


					


if __name__ == '__main__':

	xml = StitchingXML('/data/elowsky/OLST_2/data/Vglut_stitching/N5_gzip_36Z/dataset.xml~2')

	xml.analyzeCorrelations(zPair={3,4})



	


