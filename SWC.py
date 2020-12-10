import numpy as np
import math
import params
from Node import Node
from os.path import exists
from sys import exit
import tifffile as tif

"""
SWC object
	- swc is internally represented as a tree of nodes
	- implements functions for processing of swcs
"""

class SWC:

	def __init__(self, swcPath, preserveSID=True, preserveRadius=True):

		"""
		Params:	swcPath - string, absolute path of swc file 
			preserveSID - boolean, if false will set all structure id values to 0
			preserveRadius - boolean, if false will set all radius values to 0
		
		This constructor only calls buildSWCTree() to create the tree structure
		that will be used to represent the SWC

		"""		

		self.swcPath = swcPath
		self.preserveSID = preserveSID
		self.preserveRadius = preserveRadius

		self.buildSWCTree()


	def buildSWCTree(self):

		"""
		loads SWC file and builds a tree structure where each coordinate in the SWC is
		represented by a Node in the tree

		"""

		print('Building SWC Tree:', self.swcPath.split('/')[-1])

		# load swc 
		if exists(self.swcPath):
			swcArray = np.genfromtxt(self.swcPath,delimiter=params.SWC_DELIMITER,dtype=object)
		else:
			exit('ERROR: SWC Path Does Not Exist!!!')

		# check for dimension error
		if swcArray.ndim != 2:
			exit('Error: # Dimensions of SWC Array is ' + str(swcArray.ndim))

		self.fixTypes(swcArray)
		
		# build node for soma
		self.nodeSoma = self.buildNode(swcArray[swcArray[:,params.PID_INDEX]==params.SOMA_PID][0])

		stack = [self.nodeSoma]

		while len(stack):

			node = stack.pop()

			# get all children
			childRows = swcArray[swcArray[:,params.PID_INDEX] == node.id_]

			# build nodes and push child on stack
			for childRow in childRows:
				nodeChild = self.buildNode(childRow, parentNode=node)
				node.childNodes.insert(0,nodeChild)
				stack.append(nodeChild)

	def buildNode(self,swcRow,parentNode=None):

		"""
		Params: swcRow - 1D array, [id, structureID, x, y, z, radius, parentID]
			parentNode - Node

		Return: Node object populated with values from swcRow

		"""
		
		sID = swcRow[params.SID_INDEX] if self.preserveSID else 0
		radius = swcRow[params.RADIUS_INDEX] if self.preserveRadius else 0
		
		node = Node(swcRow[params.ID_INDEX],sID,swcRow[params.X_INDEX],swcRow[params.Y_INDEX],swcRow[params.Z_INDEX],radius,parentNode=parentNode)

		return node


	def fixTypes(self, swcArray):
	
		"""
		Params: swcArray, 2D array, swc array
		
		Changes data type of id,structure id,radius,parent id, to int
		Chnages data type of x,y,z to float

		"""

		swcArray[:,[params.ID_INDEX,params.SID_INDEX,params.RADIUS_INDEX,params.PID_INDEX]] = swcArray[:,[params.ID_INDEX,params.SID_INDEX,params.RADIUS_INDEX,params.PID_INDEX]].astype(float).astype(int)
		swcArray[:,[params.X_INDEX,params.Y_INDEX,params.Z_INDEX]] = swcArray[:,[params.X_INDEX,params.Y_INDEX,params.Z_INDEX]].astype(float)	


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
		Params: scaleFactors, list of 3 floats [x scale, y scale, z scale]

		Scales all coordinates
		"""
		
		if nodeRoot == None:
			nodeRoot = self.nodeSoma
		
		xScale, yScale, zScale = scaleFactors[0], scaleFactors[1] ,scaleFactors[2]

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

		self.scaleCoords(scaleFactors=[0,-1,0], nodeRoot=nodeRoot)

	def pixelsToMicrons(self,orientation, nodeRoot=None):

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

		self.scaleCoords(scaleFactors=resFactors, nodeRoot=nodeRoot)	
	


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
			numNodes = self.numNodes

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


	def generateSWCArray(self):

		"""
		Return: 2D array, SWC in 2D array
	
		"""

		out_swc = []

		id_ = 1

		stack = [(self.nodeSoma,params.SOMA_PID)]

		while len(stack):
		
			nodeTup = stack.pop()
			node, pid = nodeTup[0], nodeTup[1]
		
			out_swc.append([id_,node.structureID,node.x,node.y,node.z,node.radius,pid])			

			for c in node.childNodes:
				stack.append((c,id_))

			id_ += 1

		return np.array(out_swc)

	def save(self, outPath, fmt=params.SWC_FORMAT_FLOAT):

		"""
		Params: outPath - string, path to save SWC
			fmt - array of strings, datatype specifier for each column in SWC
			      ex: ['%d','%d','%d','%d','%d','%d','%d']

		Generates and Saves SWC

		"""

		print('Saving SWC:',outPath)
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





	

		

	def smooth(self, coord='z'):
	# (n+2) moving average where n is number of children
	### CAREFUL WITH CORONAL OR OBLIQUE ######
	# if smoothing in coronal, must smooth y (index 1)

		stack = [node_root]

		while len(stack) > 0:

			node = stack.pop(0)

			# dont smooth soma or children of soma
			if node != node_root:
				child_nodes = node.childNodes
				# dont change terminal nodes or bifurcation nodes or children of bifurcation nodes
				if len(child_nodes) == 1 and len(node.parentNode.childNodes) == 1:
					child_y_sum = sum([n.temp_y for n in child_nodes])
					updated_y = (child_y_sum + node.temp_y + node.parentNode.temp_y)/(len(child_nodes) + 2)
					node.y = updated_y

			for c in node.childNodes:
				stack.append(c)

		# set z
		for node in nodes_in_tree(node_root):
			node.temp_y = node.y
	

		return node_root

	def get_differences(node_root):
		stack = [node_root]

		dx,dy,dz = 0,0,0
		while len(stack) > 0:

			node = stack.pop(0)

			for c in node.childNodes:
				if node != node_root:
					dx += abs(node.x - c.x)
					dy += abs(node.y - c.y)
					dz += abs(node.z - c.z)
		
				stack.append(c)

		return [dx,dy,dz]

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
		


	
			


	def __str__(self):
		return self.swcPath


if __name__ == '__main__':
	swcPath = '/data/elowsky/OLST/swc_analysis/automatically_traced/flagship/morphometric_analysis/swcs/2_171012_13.swc'

	swc = SWC(swcPath)
	
	#print(swc.checkForDuplicates())
	print('# Stems:',swc.numStems())
	print('# Nodes:',swc.numNodes())
	#print('# Terminal Nodes:',swc.numTerminalNodes())
	#print('# Bifurcation Nodes:',swc.numBifurcationNodes())
	#print('EucDistances:',swc.Lmeasure_EucDistance())
	#print('Length:',swc.Lmeasure_Length())
	#print('PathDistance:',swc.Lmeasure_PathDistance())
	
	print('Maximum Radial Distance:', swc.maxRadialDistance())
	swc.sholl()






