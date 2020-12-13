import numpy as np
import math
import tifffile as tif
import copy
import matplotlib.pyplot as plt
import params
import utils
from Node import Node
from os.path import exists, join
from sys import exit
from mpl_toolkits.mplot3d import Axes3D
from skimage.io import imread


"""
SWC object
	- swc is internally represented as a tree of nodes
	- implements functions for processing of swcs
"""

class SWC:

	def __init__(self, swcPath, * , fromPath=True, swcArray=None, preserveSID=True, preserveRadius=True):

		"""
		Params:	swcPath - string, absolute path of swc file 
			fromPath - boolean, if true reads swc from file, else reads from array
			swcArray - 2d numpy arraym 
			preserveSID - boolean, if false will set all structure id values to 0
			preserveRadius - boolean, if false will set all radius values to 0
		
		This constructor only calls buildSWCTree() to create the tree structure
		that will be used to represent the SWC

		"""		

		self.swcPath = swcPath

		# build SWC tree sutructure
		self.buildSWCTree(fromPath=fromPath, swcArray=swcArray, preserveSID=preserveSID, preserveRadius=preserveRadius)


	def buildSWCTree(self, * , fromPath=True, swcArray=None, preserveSID=True, preserveRadius=True):

		"""
		loads SWC file and builds a tree structure where each coordinate in the SWC is
		represented by a Node in the tree

		Params: fromPath, boolean if true will read swc from file, otherwise will read from swcArray parameter
			swcArray, 2d array representation of swc 

		"""

		# load swc

		if fromPath: 
			if exists(self.swcPath):
				print('Building SWC Tree:', self.swcPath.split('/')[-1])
				swcArray = np.genfromtxt(self.swcPath,delimiter=params.SWC_DELIMITER,dtype=object)
			else:
				exit('ERROR: SWC Path Does Not Exist!!!')
		else:
			if swcArray is None:
				exit('ERROR: If not reading from path must provide swcArray!')

			if not isinstance(swcArray,np.ndarray):
				exit('ERROR: swcArray must be of type numpy array!')


		# dimension errors
		if swcArray.ndim != 2:
			exit('Error: # Dimensions of SWC Array is ' + str(swcArray.ndim) + ' and must be 2!')

		if swcArray.shape[1] != len(params.SWC_INDICES):
			exit('Error: # Columns of SWC Array is '+ str(swcArray.shape[1]) + ' and must be '+str(len(params.SWC_INDICES)))


		# modify data types
		swcArray = utils.fixSWCArrayTypes(swcArray)
		
		# build node for soma
		somaRow = utils.getSomaRowFromSWCArray(swcArray)
		self.nodeSoma = utils.buildNode(somaRow, preserveSID=preserveSID, preserveRadius=preserveRadius)

		stack = [self.nodeSoma]

		while len(stack):

			node = stack.pop()

			# get all children
			childRows = utils.getChildRowsFromSWCArray(swcArray, node.id_)

			# build nodes and push child on stack
			for childRow in childRows:
				nodeChild = utils.buildNode(childRow, parentNode=node, preserveSID=preserveSID, preserveRadius=preserveRadius)
				node.childNodes.insert(0,nodeChild)
				stack.append(nodeChild)


	def getNodesInTree(self, nodeRoot=None):

		"""
		Return: list of  all nodes in tree
		"""
			
		if nodeRoot == None:
			nodeRoot = self.nodeSoma

		nodes = []

		stack = [nodeRoot]

		while len(stack):

			node = stack.pop()
			nodes.append(node)

			for childNode in node.childNodes:
				stack.append(childNode)

		return nodes

	def translateCoords(self, shifts, nodeRoot=None):

		"""
		Params: shifts, list of 3 floats [x shift, y shift, z shift]

		Translates (shifts) all coordinates
		"""

		if nodeRoot == None:
			nodeRoot = self.nodeSoma
		
		xShift, yShift, zShift = shifts[0], shifts[1] ,shifts[2]

		stack = [nodeRoot]

		while len(stack):

			node = stack.pop()

			node.x += xShift
			node.y += yShift
			node.z += zShift

			for childNode in node.childNodes:
				stack.append(childNode)

	def scaleCoords(self, scaleFactors, nodeRoot=None):

		"""
		Params: scaleFactors, list of 3 floats [x scale, y scale, z scale] or just scalar 

		Scales all coordinates
		"""
		
		if isinstance(scaleFactors,float) or isinstance(scaleFactors,int):
			xScale, yScale, zScale = scaleFactors, scaleFactors, scaleFactors
		elif isinstance(scaleFactors,list):

			if len(scaleFactors) == 3:
				xScale, yScale, zScale = scaleFactors[0], scaleFactors[1] ,scaleFactors[2]
			else:
				exit('ERROR: if scaleFactors is a list it must be of length 3!')
		else:
			exit('ERROR: scaleFactors must be of type int, float, or list!')

		if nodeRoot == None:
			nodeRoot = self.nodeSoma

		stack = [nodeRoot]

		while len(stack):

			node = stack.pop()
			
			node.x *= xScale
			node.y *= yScale
			node.z *= zScale

			for childNode in node.childNodes:
				stack.append(childNode)

	def invertY(self, nodeRoot=None):

		"""
		Inverts all y coordinates
		Assumes soma is centered already
		
		"""

		self.scaleCoords(scaleFactors=[1,-1,1], nodeRoot=nodeRoot)

	def pixelsToMicrons(self,orientation, *, centerAroundSoma=False,  nodeRoot=None):

		"""
		Params: orientation - string, orientation of coordinates (coronal or oblique)

		Converts all coordinates from pixels to microns
	
		"""
		
		if orientation == 'oblique':
			resFactors = params.RES_OBLIQUE
		elif orientation == 'coronal':
			resFactors = params.RES_CORONAL
		else:
			exit('Orientation Variable Not Valid: ' + orientation)

		if centerAroundSoma:
			self.centerAroundSoma()

		self.scaleCoords(scaleFactors=resFactors, nodeRoot=nodeRoot)

	def micronsToPixels(self,orientation, somaCenter=None, nodeRoot=None):

		"""
		Params: orientation - string, orientation of coordinates (coronal or oblique)

		Converts all coordinates from microns to pixels
	
		"""
		
		if orientation == 'oblique':
			resFactors = params.RES_OBLIQUE
		elif orientation == 'coronal':
			resFactors = params.RES_CORONAL
		else:
			exit('Orientation Variable Not Valid: ' + orientation)

		self.scaleCoords(scaleFactors=list(1/np.array(resFactors)), nodeRoot=nodeRoot)

		if somaCenter != None:
			self.translateCoords(somaCenter, nodeRoot=nodeRoot)

	def applyStitchingMatrices(self, stitchingMatrices):

		self.transformCoords(stitchingMatrices['calibration'])
		self.transformCoords(stitchingMatrices['Translation to Regular Grid'])
		self.transformCoords(stitchingMatrices['Stitching Transform'])
		self.transformCoords(np.linalg.inv(stitchingMatrices['calibration']))


	def transformCoords(self, transformationMatrix):

		stack = [self.nodeSoma]

		while len(stack):

			node = stack.pop()
			nodeCoords = transformationMatrix @ np.array([node.x,node.y,node.z,1]).reshape(-1,1)
			node.x = nodeCoords.flatten()[0]
			node.y = nodeCoords.flatten()[1]
			node.z = nodeCoords.flatten()[2]


			for childNode in node.childNodes:
				stack.append(childNode)

		

	def extractTypeTree(self, structure):

		"""
		Params: structure - string, name of structure (basal dendrite, apical dendrite, etc...)
		
		Removes all branches that are not of the same type as structure. For example, if only 
		want to do analysis on basal tree, make struture = 'basal dendrite'

		"""

		if not structure in params.SID_DICT:
			exit('Structure Not In Structure Dictionary:' + structure)

		sID = params.SID_DICT[structure]

		newChildArray = [node for node in self.nodeSoma.childNodes if node.structureID == sID]
		self.nodeSoma.childNodes = newChildArray


	def obliqueToCoronal(self, fusedDimensions, shearFactor):

		stack = [self.nodeSoma]

		while len(stack):

			node = stack.pop()

			coronalCoords = utils.obliqueToCoronal([node.x,node.y,node.z], fusedDimensions, shearFactor)
			node.x = coronalCoords[0]
			node.y = coronalCoords[1]
			node.z = coronalCoords[2]

			for childNode in node.childNodes:
				stack.append(childNode)

		

	def centerAroundSoma(self):

		"""
		Translates all coordinates so soma is at (0,0,0)

		"""

		self.translateCoords(shifts=[-self.nodeSoma.x, -self.nodeSoma.y ,-self.nodeSoma.z])

	def numNodes(self):
		
		"""
		Return: int, number of nodes in SWC 		
	
		"""

		count = 0
		stack = [self.nodeSoma]

		while len(stack):

			node = stack.pop()
			count += 1

			for childNode in node.childNodes:
				stack.append(childNode)

		return count


	def numStems(self):

		"""
		Return: int, number of stems (number of children of soma)		
	
		"""
		
		return len(self.nodeSoma.childNodes)
		


	def setRadii(self, somaRadius, radius):

		"""
		Params: somaRadius - radius of soma node 
			radius - radius of all other node besides soma

		Sets the radius of all nodes in the tree

		"""

		stack = [self.nodeSoma]

		while len(stack):

			node = stack.pop()

			if node == self.nodeSoma:
				node.radius = somaRadius
			else:
				node.radius = radius

			for childNode in node.childNodes:
				stack.append(childNode)

	def setBranchIdsSame(self):

		"""
		Sets all ids to ids of branch (id of child of soma
		"""
		
		for child in self.nodeSoma.childNodes:

			branchID = child.id_
			
			stack = [child]

			while len(stack):

				node = stack.pop()

				node.id_ = branchID

				for childNode in node.childNodes:
					stack.append(childNode)

	def saveSkeleton(self, path, imageShape,branchIdsSame=False):

		""" 
		Params:
		

		Saves tiff file that is skeleton
		"""

		if branchIdsSame:
			skeleton = np.zeros(shape=imageShape, dtype=np.uint8)
			print('Saving Tiff as 8 bit (Branch Ids Same)...')
		else:
			# calculate how many nodes
			numNodes = self.numNodes()

			if numNodes <= 2**8 - 1:
				print('Saving Tiff as 8 bit...')
				skeleton = np.zeros(shape=imageShape, dtype=np.uint8)
			elif numNodes <= 2**16 - 1:
				print('Saving Tiff as 16 bit...')
				skeleton = np.zeros(shape=imageShape, dtype=np.uint16)
			elif numNodes <= 2**32 - 1:
				print('Saving Tiff as 32 bit...')
				skeleton = np.zeros(shape=imageShape, dtype=np.uint32)
			elif numNodes <= 2**64 - 1:
				print('Saving Tiff as 64 bit...')
				skeleton = np.zeros(shape=imageShape, dtype=np.uint64)

		# set soma intensity to 1
		skeleton[int(self.nodeSoma.z),int(self.nodeSoma.y),int(self.nodeSoma.x)] = 1
		 

		id_ = 2
		for child in self.nodeSoma.childNodes:

			stack = [child]

			while len(stack):

				node = stack.pop()

				if branchIdsSame:
					skeleton[int(node.z),int(node.y),int(node.x)] = id_
				else:
					skeleton[int(node.z),int(node.y),int(node.x)] = node.id_
	
				for childNode in node.childNodes:
					stack.append(childNode)	
			id_ += 1

		tif.imsave(path,skeleton)

	def saveBranchInfo(self, path):

		"""
		Params:

		Saves branch info to text file for label markup

		"""

		with open(path, 'w') as fp:
			fp.write('Soma: ' + str(int(self.nodeSoma.x)) + ',' + str(int(self.nodeSoma.y)) + ',' + str(int(self.nodeSoma.z)) + '\n\n')
			id_ = 2
			for i,child in enumerate(self.nodeSoma.childNodes):

				
				if i < len(self.nodeSoma.childNodes) - 1:
					fp.write(str(id_) + ' : ' + str(child.id_) + ' : ' + str(int(child.x)) + ' ' + str(int(child.y)) + ' ' + str(int(child.z)) + '\n')
				else:
					fp.write(str(id_) + ' : ' + str(child.id_) + ' : ' + str(int(child.x)) + ' ' + str(int(child.y)) + ' ' + str(int(child.z)))
				id_ += 1


	def generateSWCArray(self, onlyCoords=False):

		"""
		Return: 2D array, SWC in 2D array
	
		"""

		outSWC = []

		id_ = 1

		stack = [(self.nodeSoma,params.SOMA_PID)]

		while len(stack):
		
			nodeTup = stack.pop()
			node, pid = nodeTup[0], nodeTup[1]
		
			outSWC.append([id_,node.structureID,node.x,node.y,node.z,node.radius,pid])			

			for c in node.childNodes:
				stack.append((c,id_))

			id_ += 1

		if onlyCoords:
			return np.array(outSWC)[:,[params.SWC_INDICES['x'], params.SWC_INDICES['y'], params.SWC_INDICES['z']]]
		else:
			return np.array(outSWC)

	def saveSWC(self, outPath, fmt=params.SWC_FORMAT_FLOAT):

		"""
		Params: outPath - string, path to save SWC
			fmt - array of strings, datatype specifier for each column in SWC
			      ex: ['%d','%d','%d','%d','%d','%d','%d']

		Generates and Saves SWC

		"""

		print('Saving SWC:',outPath)

		# remove duplicates
		self.removeDuplicates()

		swcArray = self.generateSWCArray()
		np.savetxt(outPath,swcArray,delimiter=params.SWC_DELIMITER,fmt=fmt)





	def heightWidthDepth(self):

		"""
		Return: float,float,float

		Calculates range of coordinates in all dimensions
			
		"""

		minX, maxX = float('inf'), float('-inf')
		minY, maxY = float('inf'), float('-inf')
		minZ, maxZ = float('inf'), float('-inf')

		stack = [self.nodeSoma]

		while len(stack):

			node = stack.pop()

			if node.x < minX:
				minX = node.x
			if node.x > maxX:
				maxX = node.x
			if node.y < minY:
				minY = node.y
			if node.y > maxY:
				maxY = node.y
			if node.z < minZ:
				minZ = node.z
			if node.z > maxZ:
				maxZ = node.z
		
			for childNode in node.childNodes:
				stack.append(childNode)

		return maxX-minX,maxY-minY,maxZ-minZ


	def distance(self, a, b):

		"""
		Params: a - Node
			b - Node

		Return: float, euclidean distance between two nodes


		"""

		return math.sqrt((a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2)


	def Lmeasure_EucDistance(self):
		
		"""
		Return: float, float, float, float, float

		Recreation of EucDistance function from 'lmeasure'
		Distances are calculated from soma to every node

		"""

		distances = []

		stack = [self.nodeSoma]

		while len(stack):

			node = stack.pop()
			distances.append(self.distance(node,self.nodeSoma))

			for childNode in node.childNodes:
				stack.append(childNode)
		
		totalDistance = np.sum(distances)
		minDistance = np.min(distances)
		avgDistance = np.mean(distances)
		maxDistance = np.max(distances)
		stdDistance = np.std(distances)	

		return totalDistance, minDistance, avgDistance, maxDistance, stdDistance

	def Lmeasure_Length(self):

		"""
		Return: float, float, float, float, float

		Recreation of Length function from 'lmeasure'
		Distances are calculated from soma to every node

		"""

		distances = []

		stack = [self.nodeSoma]

		while len(stack):

			node = stack.pop()

			for childNode in node.childNodes:

				distances.append(self.distance(childNode,node))
				stack.append(childNode)
		
		totalDistance = np.sum(distances)

		minDistance = np.min(distances)
		avgDistance = np.mean(distances)
		maxDistance = np.max(distances)
		stdDistance = np.std(distances)	

		return totalDistance, minDistance, avgDistance, maxDistance, stdDistance

	def Lmeasure_PathDistance(self):

		"""
		Return: float, float, float, float, float

		Recreation of PathDistance function from 'lmeasure'
		Distances are calculated from soma to every node

		"""

		distances = []

		stack = [(self.nodeSoma,0)]

		while len(stack):

			nodeTup = stack.pop()
			node = nodeTup[0]
			dist = nodeTup[1]

			distances.append(dist)

			for childNode in node.childNodes:
				stack.append((childNode,dist+self.distance(childNode,node)))
		
		totalDistance = np.sum(distances)
		minDistance = np.min(distances)
		avgDistance = np.mean(distances)
		maxDistance = np.max(distances)
		stdDistance = np.std(distances)	

		return totalDistance, minDistance, avgDistance, maxDistance, stdDistance


	
	def terminalNodes(self,nodeRoot=None):

		if nodeRoot == None:
			nodeRoot = self.nodeSoma

		terminalNodes = []
		stack = [self.nodeSoma]

		while len(stack):

			node = stack.pop()
			
			if len(node.childNodes) == 0:
				terminalNodes.append(node)

			for childNode in node.childNodes:
				stack.append(childNode)

		return terminalNodes
	
	def numTerminalNodes(self, nodeRoot=None):

		"""
		Return: int, number of terminal nodes

		"""

		return len(self.terminalNodes(nodeRoot=nodeRoot))

	def numBifurcationNodes(self):

		"""
		Return: int, number of bifucrcation nodes (not including soma)

		"""

		numBifurcationNodes = 0
		stack = [self.nodeSoma]

		while len(stack):

			node = stack.pop()
			
			if len(node.childNodes) > 1 and node != self.nodeSoma:
				numBifurcationNodes += 1

			for childNode in node.childNodes:
				stack.append(childNode)

		return numBifurcationNodes

	def collapseSoma(self,collapseRadius=None):

		"""
		Params: collapseRadius, float - remove all nodes within collapseRadius from soma

		"""
		

		if collapseRadius == None:
			collapseRadius = params.SOMA_COLLAPSE_RADIUS
		
		print('Collapse Soma: ' + str(collapseRadius) + ' (um)')

		innerNodes,outerNodes = self.getInnerOuterNodes(collapseRadius)
		innerNodes = [innerNode for innerNode in innerNodes if innerNode != self.nodeSoma]
		outer_to_inner_nodes = [innerNode for innerNode in innerNodes if innerNode.parentNode in outerNodes]
		print('# Nodes (before):',self.numNodes())
		print('# Inner Nodes:',len(innerNodes))
		print('# Outer Nodes:',len(outerNodes))
		print('# Outer to Inner Nodes:',len(outer_to_inner_nodes))


		# get all terminal nodes that are outside radius
		terminalNodes = [terminalNode for terminalNode in self.terminalNodes() if self.distance(terminalNode, self.nodeSoma) >= collapseRadius]

		# track backwards from all terminal nodes until reaches node in radius
		newSomaChildNodes = []

		for terminalNode in terminalNodes:
			
			pointer = terminalNode

			while pointer:
		
				# parent reference might already be removed
				if pointer.parentNode == None:
					break

				if self.distance(self.nodeSoma, pointer.parentNode) < collapseRadius:
					pointer.parentNode = None
					newSomaChildNodes.append(pointer)
					break
		
				pointer = pointer.parentNode

		# reset somas child nodes
		self.nodeSoma.childNodes = newSomaChildNodes
		for child in self.nodeSoma.childNodes:
			child.parentNode = self.nodeSoma

				
		# Error Checks
		innerNodes = [node for node in self.getNodesInTree(self.nodeSoma) if self.distance(node,self.nodeSoma) < collapseRadius and node != self.nodeSoma]
		if len(innerNodes) > 0:				
			exit('Error: Node found < ' + str(collapseRadius) + ' microns from soma')

		print('# Nodes (after):', self.numNodes())


	def removeDuplicates(self):

		"""
		Removes duplicate nodes in treee
		Nodes are considered duplicates if same coordinates and are in parent-child relationship
		

		"""

		stack = [self.nodeSoma]

		while len(stack):

			node = stack.pop()

			if node != self.nodeSoma:
				if self.areCoordsSame(node, node.parentNode):
					
					print('REMOVING DUPLICATE...')

					node.parentNode.childNodes = node.childNodes

					# set childrens parent pointers to grandparent
					for childNode in node.childNodes:
						childNode.parentNode = node.parentNode

			for childNode in node.childNodes:
				stack.append(childNode)


	def areCoordsSame(self, a, b):

		"""
		Params: a - Node
			b - Node

		Return: boolean, True if x,y,z coordinates of nodes are the same

		"""

		return a.x == b.x and a.y == b.y and a.z == b.z 

	def checkForDuplicates(self):

		"""
		Return: list, contains ids of duplcate nodes

		"""

		coordsVisited, duplicateIDs = [], []

		stack = [self.nodeSoma]

		while len(stack):

			node = stack.pop()

			if node != self.nodeSoma:
				if self.areCoordsSame(node, node.parentNode):
					duplicateIDs.append(node.id_)

			for childNode in node.childNodes:
				stack.append(childNode)

		return duplicateIDs


	def smooth(self, orientation):

		"""
		Params: orientation, string, 'corona'l or 'oblique', oblique smooths z, coronal smooths y
			

		- smooths coordinates to  (non biological) sharp edges
		- doesnt change soma, children of soma, terminal nodes, bifurcation nodes, or children of bifurcation nodes
		
		"""

		if orientation != 'coronal' and orientation != 'oblique':
			exit('Error: Invalid Orientation!')

		diff = self.smoothingDiffs(orientation)

		diffChange, iters = np.Inf, 0

		while diffChange > params.SMOOTHING_DIFF_THRESHOLD or iters < params.SMOOTHING_ITERS_THRESHOLD:

			# smooth
			stack = [self.nodeSoma]

			while len(stack):

				node = stack.pop()
	
				# make sure not soma and has 1 child (discludes bifurcation and terminal nodes
				if node != self.nodeSoma and len(node.childNodes) == 1 and len(node.parentNode.childNodes) == 1:

					if orientation == 'coronal':
						node.tmp = (node.parentNode.y + node.y + node.childNodes[0].y)/3
					elif orientation == 'oblique':
						node.tmp = (node.parentNode.z + node.z + node.childNodes[0].z)/3

				for childNode in node.childNodes:
					stack.append(childNode)

			# update coordinates
			for node in self.getNodesInTree():
				if node.tmp:
					if orientation == 'coronal':
						node.y = node.tmp
					elif orientation == 'oblique':
						node.z = node.tmp


			# recaulate differences
			diffSmoothed = self.smoothingDiffs(orientation)
			diffChange = abs(diffSmoothed - diff)
			diff = diffSmoothed
			iters += 1
			print(diffChange)
		print('# iters:',iters)


	def smoothingDiffs(self,orientation):

		stack = [self.nodeSoma]

		diff = 0

		while len(stack):

			node = stack.pop()

			for childNode in node.childNodes:
				if node != self.nodeSoma:
					if orientation == 'coronal':
						diff += abs(node.y - childNode.y)
					elif orientation == 'oblique':
						diff += abs(node.z - childNode.z)
				stack.append(childNode)

		return diff

	def sholl(self, radialPercentIncrement = 5):

		maxRadialDistance = self.maxRadialDistance()
		
		for radiusPercentage in range(0,101,radialPercentIncrement):
			if radiusPercentage == 0:
				continue
			radius = radiusPercentage*maxRadialDistance/100

			print('Radius (%):',radiusPercentage)
			print('Radius (um):',radius)
	
			globalLengthInner, globalLengthTotal = 0,0
			

			# for all branches (children of soma)
			for stemNode in self.nodeSoma.childNodes:
				
				# accumulate path distance and only add to sum if intersects sphere
				intersect = False
				length = self.distance(stemNode, self.nodeSoma)
				totalLength = length
				stack = [stemNode]

				while len(stack):

					node = stack.pop()

					if self.distance(node,self.nodeSoma) > radius:
						intersect = True

					for childNode in node.childNodes:
						
						if self.distance(node,self.nodeSoma) < radius:
							length += self.distance(childNode,node)
						totalLength += self.distance(childNode,node)
						stack.append(childNode)

				if intersect:
					globalLengthInner += length
				globalLengthTotal += totalLength

			print(globalLengthInner / globalLengthTotal)
				

			print('##########')


	def maxRadialDistance(self):

		"""
		Return: float, maximal radial distance froma soma 		
	
		"""
		maxRadialDist = -1
		stack = [self.nodeSoma]

		while len(stack):

			node = stack.pop()
			
			dist = self.distance(self.nodeSoma, node)
			if dist > maxRadialDist:
				maxRadialDist = dist

			for childNode in node.childNodes:
				stack.append(childNode)

		return maxRadialDist
		

	def getInnerOuterNodes(self, collapseRadius):

		innerNodes = []
		outerNodes = [] 

		stack = [self.nodeSoma]

		while len(stack):

			node = stack.pop()
			
			if self.distance(node, self.nodeSoma) < collapseRadius:
				innerNodes.append(node)
			else:
				outerNodes.append(node)

			for childNode in node.childNodes:
				stack.append(childNode)

		if len(innerNodes) + len(outerNodes) != self.numNodes():
			sys.exit('Error: Inner nodes + Outer nodes != All nodes')

		return innerNodes, outerNodes
		
	def setStructureLabels(self, labelInfo):
		
		self.nodeSoma.structureID = 1

		for childNode in self.nodeSoma.childNodes:
			
			for labelTup in labelInfo:
				sID = None
				if labelTup[0] == childNode.id_:
					sID = labelTup[1]
					break

			if sID == None:
				exit('Error: Soma Child with Label Not Found!')

			stack = [childNode]

			while len(stack):

				node = stack.pop()
			
				node.structureID = sID

				for childNode in node.childNodes:
					stack.append(childNode)
		
		# Error Check
		nodes = self.getNodesInTree()
		if sum([node for node in nodes if node.structureID == 0]) > 0:
			exit('Some nodes do not have label!')

	def register(self,brain):

		"""
		Assumes SWC is in coronal orientation
		"""

		cropDict = utils.loadCropInfo()

		# apply cropping shift
		shift = cropDict[brain]
		self.translateCoords([0,-shift*params.CROP_FACTOR,0])
				
		# get coordinates
		swcCoords = self.generateSWCArray(onlyCoords=True)

		brainPath = join(params.REGISTRATION_BASE_PATH, brain)

		# write coords to file
		outPath = join(brainPath,'swc_in.txt')
		utils.writeCoordsForRegistration(outPath, swcCoords)

		# run transformix
		tpPath = join(brainPath,'TransformParameters_labels.1_full_scale.txt')
		inPath = join(brainPath,'swc_in.txt')
		utils.runTransformix(brainPath, tpPath, inPath)

		# read in registered coordinates
		transformedCoordsPath = join(brainPath,'outputpoints.txt')
		registeredCoords = utils.parseTransformixOutput(transformedCoordsPath)

		swcArray = self.generateSWCArray()
		swcArray[:,[params.X_INDEX, params.Y_INDEX, params.Z_INDEX]] = registeredCoords

		# re build SWC Tree
		self.buildSWCTree(fromPath=False, swcArray=swcArray)

		self.scaleCoords([params.REGISTRATION_RES_MICRONS])


	def overlayOnReferenceSpace(self, outPath, * , downsampling=25):

		referenceImage = imread(params.MOP_REFERENCE_PATHS[downsampling])

		swc.scaleCoords([1/downsampling])

		swcArray = self.generateSWCArray(onlyCoords=True).astype(int)

		referenceImage[swcArray[:,2],swcArray[:,1],swcArray[:,0]] = params.OVERLAY_INTENSITY

		tif.imsave(outPath, referenceImage)

		

		
				
			
			

	def plot(self,dim=2):
	
		swcArray = self.generateSWCArray()		
		uniqueSID = np.unique(swcArray[:,params.SID_INDEX]).astype(int)		

		if dim == 2:
			for sid in uniqueSID:
				sidCoords = swcArray[swcArray[:,params.SID_INDEX] == sid]
				plt.scatter(sidCoords[:,params.X_INDEX],-sidCoords[:,params.Y_INDEX], color=params.SID_COLORS[sid], s=params.SID_PLOTTING_RADII[sid])

			plt.show()
		elif dim == 3:

			fig = plt.figure()
			ax = fig.add_subplot(111,projection='3d')
			for sid in uniqueSID:
				sidCoords = swcArray[swcArray[:,params.SID_INDEX] == sid]
				ax.scatter(sidCoords[:,params.X_INDEX],-sidCoords[:,params.Y_INDEX],sidCoords[:,params.Z_INDEX], color=params.SID_COLORS[sid], s=params.SID_PLOTTING_RADII[sid])
			plt.show()			

		else:
			exit('Error: # Dimensions must be 2 or 3')
	



	def __str__(self):
		return self.swcPath


if __name__ == '__main__':

	swcPath = '/data/elowsky/OLST/swc_analysis/automatically_traced/flagship/layer_6/190416/registered/190416_23.swc'
	swc = SWC(swcPath)
	swc.scaleCoords(2)
	

	







