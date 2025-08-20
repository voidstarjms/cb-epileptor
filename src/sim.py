from brian2 import *
from matplotlib import pyplot as plt
import argparse

def main():
    pop1_eqs = '''
    dx/dt = y - a * x ** 3 + b * x ** 2 - z + I_app + C_E * (x_bar - x) +
        sigma * I_syn_rflx + sigma * I_syn_ext + W(t) : 1
    dy/dt = c - d * x ** 2 - y : 1
    dz/dt = r * (s * (x + x2_bar - x_naught)) : 1

    a : 1
    b : 2
    c : 1
    d : 1
    s : 1
    I : 1
    x_naught : 1

    # This needs to be passed in throughout runtime
    x2_bar : 1
    '''

    pop2_eqs = '''
    dv/dt = (I_app - gL * (v - E_L) - gK * n * (V - E_K) - gCa * m_inf *
        (v - E_Ca) + sigma * C_E * (x_bar - x) + I_syn_ext + I syn_rflx -
        sigma * 0.3 * (z_bar - 3) + sigma * W(t)) / Cm : volt
    dn/dt = phi * (n_inf - n) / tau_n : 1

    m_inf = 0.5 * (1 + tanh((v - v1) / v2) : 1
    tau_n = 1 / cosh((v - v3) / (2 * v4)) : 1
    n_inf = 0.5 * (1 + tanh((v - v3) / v4)) : 1

    x = v/20 : 1

    Cm : ufarad/cm2
    I : amp
    v1 : volt
    v2 : volt
    v3 : volt
    v4 : volt
    phi : Hz
    E_Ca : volt
    E_K : volt
    E_L : volt
    '''
