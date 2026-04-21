from brian2 import *

# --- Simulation Control ---
SIM_DURATION = 120 * second
NUM_CELLS = 10
TAU_CLOCK = 1 * msecond
DT_SCALING = 20  # defaultclock.dt = TAU_CLOCK / DT_SCALING
TRANSIENT = 10

# --- Coupling & Global Logic ---
# ISOLATE = 0: decoupled, 1: coupled
ISOLATE = 1
COUPLING_STRENGTH = 0.1 # - Timed array
W_MAX = 0.006

# --- Population 1: Hindmarsh-Rose (HR) ---
HR_A = 1.0
HR_B = 3.0
HR_C = 1.0
HR_D = 5.0
HR_S = 4.0
HR_I_APP = 6
HR_X_NAUGHT = -3.5 # - timed array
HR_R = 0.00002 / msecond
HR_SIGMA = 1/50
HR_THRESHOLD = 'x > 1.5'
HR_REFRACTORY_CONDITION = 'x >= 0' 
I_SCALE = 0.05 * uamp

# --- Population 2: Morris-Lecar (ML) ---
ML_CM = 20 * ufarad
ML_I_APP = 37 * uamp
ML_V1 = -1.2 * mvolt
ML_V2 = 18 * mvolt
ML_V3 = 12 * mvolt
ML_V4 = 17.4 * mvolt
ML_PHI = 0.067 / msecond
ML_E_CA = 120 * mvolt
ML_E_K = -84 * mvolt
ML_E_L = -60 * mvolt
ML_GL = 2 * msiemens
ML_GCA = 4.0 * msiemens
ML_GK = 8.0 * msiemens
ML_SIGMA = 50 * uA     # Corresponding to sigma_2
ML_THRESHOLD = 'x > 0.95'
ML_REFRACTORY_CONDITION = 'x >= 0'

# --- Synaptic Parameters ---
SYN_VT = 2 * mV
SYN_KP = 5 * mV
SYN_TMAX = 1 * mmolar

# Excitatory (Pop 1)
SYN_ALPHA_EXC = 1.1 / (mmolar * msecond)
SYN_BETA_EXC = 0.19 / msecond
SYN_E_EXC = 0 * mV

# Inhibitory (Pop 2)
SYN_ALPHA_INH = 5 / (mmolar * msecond)
SYN_BETA_INH = 0.18 / msecond
SYN_E_INH = -80 * mV

# Conductances - Replaced with Timed arrays
# G_INTRA = 0.1 * uS  
# G_INTER = 0.2 * uS

# Plasticity
THETA_LTD_START = 0.25 # 0.05
THETA_LTD_END = 0.5 # 0.1
THETA_LTP_START = 0.75 # 0.15
A_LTD = 0.5
A_LTP = 2

TAU_WPRE = 5 * second
TAU_CA = 200 * msecond

# Timed array schedules
X_NAUGHT_VALS = [-2.5]
COUPLING_VALS = [0]
# G_INTER_VALS = [1, 1, 4, 4, 1] * uS
# G_INTRA_VALS = [0.5, 2, 2, 0.5, 0.5] * uS
G_INTER_VALS = [1] * uS
G_INTRA_VALS = [1] * uS

ML_Z_BAR_SCALE = 0.3
ML_Z_BAR_OFFSET = 6

# sigmoid calcium function parameters
CA_SIGMOID_SHIFT = 0.8
CA_SIGMOID_SLOPE = 0.2

