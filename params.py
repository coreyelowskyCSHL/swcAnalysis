# Paths
INITIAL_SWC_BASE_PATH = '/data/palmer/data/reconstructed_brains/'
SOMA_BASE_PATH = '/data/elowsky/OLST/soma_detection/'
SWC_BASE_PATH = '/data/elowsky/OLST/swc_analysis/automatically_traced/flagship/'
REGISTRATION_BASE_PATH = '/data/elowsky/OLST/registration/'

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

# Shearing Factor for Oblique -> Coronal
SHEAR_FACTOR_ISOTROPIC = -.7 
SHEAR_FACTOR_ANISOTROPIC = SHEAR_FACTOR_ISOTROPIC / ANISOTROPY_FACTOR

# All Brains
BRAINS = ['170329_500','171012','180206','180523','180606','180614','180926','181004','181115','190123','190306','190327','190416','190522']

# Volume Sizes
X_SIZE = 2048
Y_SIZE = 1024
Z_SIZE = {'170329_500':4000,'171012':4000, '180206':4000, '180523':4000, '180606':4000, '180614':4000, '180926':4000, '181004':4000, '190306':4000, '190327':4000, '190416':4000, '190522':4000,'181115':3600, '190123':3600}

# SWC Data Type Formats
SWC_FORMAT_INT = ['%d','%d','%d','%d','%d','%d','%d']
SWC_FORMAT_FLOAT = ['%d','%d','%g','%g','%g','%d','%d']

# SWC Indices
SWC_INDICES = {'id':0 , 'sid':1 , 'x':2 , 'y':3 ,'z':4 ,'radius':5 ,'pid':6 }

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

# Smoothing Parameters
SMOOTHING_DIFF_THRESHOLD = .075
SMOOTHING_ITERS_THRESHOLD = 10

# Reigstration Parameters
CROP_INFO_PATH = '/data/elowsky/OLST/registration/crop_info.txt'
CROP_FACTOR = 10
REGISTRATION_RES_MICRONS = 25
MOP_REFERENCE_PATHS = {10:'/data/elowsky/OLST/registration/MOpul_Layers_528_10x10x10.tif',25:'/data/elowsky/OLST/registration/MOpul_Layers_528.tif'}
OVERLAY_INTENSITY = 100






