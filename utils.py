import numpy as np
import copy
import params
import xml.etree.ElementTree as ET


###### Utility functions for Oblique -> Coronal Transformation #######

# These functions assume transformations are occuring on images
# therefore, the size of the image must be provided and the
# input coordinates are assumed to be in pixel coordinates of the oblique image

def obliqueToCoronal(coords, imageLengths):

	# convert to numpy arrays
	coords = np.array(coords)
	imageLengths = np.array(imageLengths)		

	# perform transformations
	coords, imageLengths = reslice(coords,imageLengths)
	coords = verticalFlip(coords, imageLengths)
	coords, imageLengths = shear(coords,imageLengths)
	coords, imageLengths = reslice(coords,imageLengths)
	coords, imageLengths = rotate(coords,imageLengths)
	
	return list(coords)
		

def reslice(coords, imageLengths):

	return coords[[1,2,0]], imageLengths[[1,2,0]]


def verticalFlip(coords, imageLengths):
	
	coords[1] = imageLengths[1] - coords[1] - 1

	return coords


def shear(coords, imageLengths):

	imageLengthsShear = copy.deepcopy(imageLengths)
	imageLengthsShear[1] = imageLengths[1] + abs(params.SHEAR_FACTOR)*imageLengthsShear[0]
	imageLengthsShear = np.round(imageLengthsShear).astype(int)

	coords[1] = np.round(coords[1] - imageLengths[1]/2 + params.SHEAR_FACTOR*(coords[0]-imageLengths[0]/2) + imageLengthsShear[1]/2 )

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



# Extracts stitching information from xml file generated from BigStitcher
def extract_stitching_parameters(xml_file):

	# read in xml
	tree = ET.parse(xml_file)
	root = tree.getroot() 

	# output lists
	files = []
	registrations = []
	stitchings = []

	# find nodes for files, registrations, stitchings
	for child in root.iter():
		if child.tag == 'files':
			files_node = child
		elif child.tag == 'ViewRegistrations':
			registrations_node = child
		elif child.tag == 'StitchingResults':
			stitchings_node = child

	for child in files_node:
		setup_number = child.attrib['view_setup']
		file_name = child[0].text
		dict_data = {"setup number":setup_number,"file name":file_name}
		files.append(dict_data)
	
	for child in registrations_node:
		setup_number = child.attrib['setup']
		stitching_transform = np.fromstring(child[0][1].text, sep=' ')
		translation_regular_grid = np.fromstring(child[1][1].text, sep=' ')
		calibration = np.fromstring(child[2][1].text, sep=' ')
		dict_data = {"setup number":setup_number,"stitching transform":stitching_transform,"translation to regular grid":translation_regular_grid,"calibration":calibration}
		registrations.append(dict_data)

	for child in stitchings_node:
		setup_a_number = child.attrib['view_setup_a']
		setup_b_number = child.attrib['view_setup_b']	
		shift = np.fromstring(child[0].text, sep=' ')
		bounding_box = np.fromstring(child[3].text, sep=' ')
		dict_data = {"setup number a":setup_a_number,"setup number b":setup_b_number,"shift":shift,"bounding box":bounding_box}
		stitchings.append(dict_data)

	return files, registrations, stitchings

def get_stitching_matrices(files,registrations,volume):

	# Get associated setup id
	setup_id = -1
	for f in files:
		if volume in f['file name']:
			setup_id = f['setup number']
			break

	# Make sure volume was found in XML
	if setup_id == -1:
		print("Error: Volume name not found in XML")
		return -1

	# Get associated stitching parameters
	stitch_params = -1
	for r in registrations:
		if r['setup number'] == setup_id:
			stitch_params = r
			break

	# Make sure volume was found in XML
	if stitch_params == -1:
		print("Error: Volume stitching parameters name not found in XML")
		return -1

	# Extract affine transformation matrices and add fourth row (matrices are 3x4)
	translation_to_grid_matrix = np.vstack([np.reshape(stitch_params['translation to regular grid'], (3,4)),np.array([0,0,0,1])])
	stitching_matrix = np.vstack([np.reshape(stitch_params['stitching transform'], (3,4)),np.array([0,0,0,1])])
	calibration_matrix = np.vstack([np.reshape(stitch_params['calibration'], (3,4)),np.array([0,0,0,1])])

	return translation_to_grid_matrix, stitching_matrix, calibration_matrix
