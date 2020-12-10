# Oblique Full Scale Resolution
X_RES_OBLIQUE = .406
Y_RES_OBLIQUE  = .406
Z_RES_OBLIQUE  = 2.5
RES_OBLIQUE = [X_RES_OBLIQUE,Y_RES_OBLIQUE,Z_RES_OBLIQUE]

# Coronal Full Scale Resolution
X_RES_CORONAL = .406
Y_RES_CORONAL = 2.5
Z_RES_CORONAL = .406
RES_CORONAL = [X_RES_CORONAL,Y_RES_CORONAL,Z_RES_CORONAL]

# Anisotropy factor
ANISOTROPY_FACTOR = Z_RES_OBLIQUE/X_RES_OBLIQUE

# All brains, Brains with 4000 z size, Brains with 3600 z size
BRAINS = ['170329_500','171012','180206','180523','180606','180614','180926','181004','181115','190123','190306','190327','190416','190522']
BRAINS_4000 = ['170329_500','171012', '180206', '180523', '180606', '180614', '180926', '181004', '190306', '190327', '190416', '190522']
BRAINS_3600 = ['181115', '190123']

# SWC Data Type Formats
SWC_FORMAT_INT = ['%d','%d','%d','%d','%d','%d','%d']
SWC_FORMAT_FLOAT = ['%d','%d','%g','%g','%g','%d','%d']

# SWC Indices
ID_INDEX = 0
SID_INDEX = 1
X_INDEX = 2
Y_INDEX = 3
Z_INDEX = 4
RADIUS_INDEX = 5
PID_INDEX = 6

# Soma Parent ID
SOMA_PID = -1

# SWC delimiter
SWC_DELIMITER = ' '

# Structure ID Dictionary
SID_DICT = {'soma': 1, 'axon': 2, 'basal dendrite': 3, 'apical dendrite': 4}

# Structure Color Dictionary
SID_COLORS = {0: 'green',1: 'black', 2: 'yellow', 3: 'red', 4: 'blue',5:'purple'}
SID_PLOTTING_RADII = {0:1,1:100,2:1,3:1,4:1,5:1}

# soma collapse radius
SOMA_COLLAPSE_RADIUS = 30

# Lmeasure Indicies Dict
LMEASURE_INDICES = {'Total_Sum':2,'Minimum':5,'Average':6,'Maximum':7,'S.D.':8}

# Order of Lmeasure Metrics
METRIC_ORDER = ['Bif_ampl_remote', 'Bif_tilt_remote', 'Branch_Order', 'Branch_pathlength', 'Contraction', 'Depth', 'EucDistance', 'Height', 'Length', 'N_stems', 'Partition_asymmetry', 'PathDistance', 'TerminalSegment', 'Terminal_degree', 'Width']

METRIC_DISPLAY_DICT = {'PathDistance':'Path Distance', 'Branch_pathlength':'Branch Path Length', 'N_stems':'N Stems', 'Length':'Length', 'EucDistance':'EucDistance', 'Height':'Height'}

SHEAR_FACTOR = -.7 / ANISOTROPY_FACTOR

# Volume Sizes
X_SIZE = 2048
Y_SIZE = 1024

