from brian2 import *

# Simulation Control
SIM_DURATION = 60 * second
NUM_CELLS = 40
DT_SCALING = 20 # Used to scale defaultclock.dt = tau / DT_SCALING
TAU_CLOCK = 1 * msecond


# Coupling & Noise
# ISOLATE = 0: decoupled, 1: coupled
ISOLATE = 1
COUPLING_STRENGTH = 1
# Note: Paper discrepancy: Table says [1, 20], Text says [2.5, 20]mV
W_MAX = 0.006          


# Population 1: Hindmarsh-Rose
HR_A = 1.0
HR_B = 3.0
HR_C = 1.0
HR_D = 5.0
HR_S = 8.0
HR_I_APP = 3.1
HR_X_NAUGHT = -3
HR_R = 0.00002 / msecond
HR_SIGMA = 1/50
HR_THRESHOLD = 'x > 1.5'


# Population 2: Morris-Lecar
ML_CM = 20 * ufarad
ML_I_APP = 40 * uamp
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
ML_SIGMA = 50 * uA
ML_THRESHOLD = 'x > 0.95'


# Synaptic Parameters
SYN_VT = 2 * mV
SYN_KP = 5 * mV
SYN_TMAX = 1 * mmolar

# Excitatory
SYN_ALPHA_EXC = 1.1 / (mmolar * msecond)
SYN_BETA_EXC = 0.19 / msecond
SYN_E_EXC = 0 * mV

# Inhibitory
SYN_ALPHA_INH = 5 / (mmolar * msecond)
SYN_BETA_INH = 0.18 / msecond
SYN_E_INH = -80 * mV

# Conductances
G_INTRA = 0.1 * uS  # Within population
G_INTER = 0.2 * uS  # Between populations