class Node:

	def __init__(self,id_,structureID,x,y,z,radius,parentNode=None):

		self.id_ = id_
		self.structureID = structureID
		self.x, self.y, self.z = x,y,z
		self.radius = radius

		self.parentNode = parentNode
			

		self.childNodes = []

		# temporary coord for smoothing
		self.tmp = None


	def getCoords(self):
		return [self.x, self.y, self.z]

	def __str__(self):
		return 'id: '+str(self.id_)+'\ncoords: ('+str(self.x)+','+str(self.y)+','+str(self.z)+')'
