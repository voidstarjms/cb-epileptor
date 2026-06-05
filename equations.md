# This file details all the equations and parameters used in the model

## Global Simulation Control
| Parameter | Value | Description |
|-----------|-------|-------------|
| `SIM_DURATION` | $120\mathrm{s}$ | Total simulated time |
| `NUM_CELLS` | $10$ | Cells per population |
| `TAU_CLOCK` | $1\mathrm{msec}$ | Base clock step |
| `DT_SCALING` | $20$ | Integration step multiplier |
| `TRANSIENT` | $10\mathrm{s}$ | Discarded transient |

## Coupling and Global Logic
| Parameter | Value | Description |
|-----------|-------|-------------|
| `ISOLATE` | $1$ | $0$: decoupled, $1$: coupled |
| $J$ | $0.1$ | Intrapopulation coupling gain |
| $W_\text{max}$ | $0.006$ | Maximum noise amplitude |

## Neuron Populations

### Excitatory (Hindmarsh-Rose)
$${x}' = y - a x^3 + b x^2 - z + I_{\text{app,1}} + J (\bar{x} - x) + \sigma_1 I_\text{syn,tot} + W(t)$$

Unitless membrane potential. Evolves according to cubic derivative modulated by y and z variables. Stimulus comes from a constant current, a stochastic current, and synaptic currents. Intrapopulation gap junction coupling simulated by partially coupling neuron potential to population mean.

$$y' = c - d x^2 - y$$

Spiking variable. Models ion transport by voltage-gated sodium and potassium ion channels.

$$z' = r(s(x + \bar{V}^\* - \eta) - \bar{z})$$

Slow adaptation variable. Introduces chaotic dynamics influenced by individual excitatory neuron potential and mean inhibitory population potential.

### Excitatory Neuron Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| $a$ | $1.0$ | x nullcline cubic scalar |
| $b$ | $3.0$ | x nullcline quadratic scalar |
| $c$ | $1.0$ | y nullcline constant element |
| $d$ | $5.0$ | y nullcline quadratic scalar |
| $s$ | $4.0$ | z excitability scalar |
| $I_{\text{app,1}}$ | $6$ | Applied current |
| $\eta$ | $-3.5$ | Resting potential |
| $r$ | $2\times10^{-5}\,(\mathrm{ms})^{-1}$ | Slow adaptation rate |
| $\sigma_1$ | $0.02$ | Noise amplitude |
| Threshold | $x > 1.5$ | Spike condition |
| Refractory | $x \geq 0$ | Reset condition |
| $\psi$ | $0.05\,\mu\mathrm{A}$ | Current scaling |

### Inhibitory (Morris-Lecar)
$$C_mV' = I_{\text{app,2}} - g_L(V - E_L) - g_K n(V - E_K) - g_{Ca} m_\infty(V)(V - E_{Ca})
    + \sigma_2((\bar{V}^\* - V^{\*}) + I_\text{syn,tot} - 0.15(\bar{z} - 6) + W(t))$$

Membrane voltage. Evolves according to ion channel currents similar to the Hodgkin-Huxley model. Intrapopulation gap junction coupling simulated by partially coupling neuron potential to population mean. Input from excitatory population slow variable drives slow dynamics.

$$n = \phi(n_\infty(V) - n)/\tau_n(V)$$

Fraction of open potassium channels.

$$m_\infty(V) = \frac{1}{2}[1 + \tanh((V-h_\text{Ca})/\lambda_\text{Ca})]$$

Fraction of open calcium channels.

$$\tau_n(V) = \frac{1}{\cosh((V-h_\text{K})/ (2 \lambda_\text{K}))}$$

Potassium channel opening time scale.

$$V^{*} = \frac{V}{20}$$

Scaled unitless voltage. For use in the synapses and excitatory population.

### Inhibitory Neuron Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| $C_m$ | $20\mu\mathrm{F}$ | Membrane capacitance |
| $I_{\text{app,2}}$ | $37\mu\mathrm{A}$ | Applied current |
| $h_{\text{Ca}}$ | $-1.2\mathrm{mV}$ | Ca activation midpoint |
| $\lambda_{\text{Ca}}$ | $18\mathrm{mV}$ | Ca activation slope |
| $h_{\text{K}}$ | $12\mathrm{mV}$ | K activation midpoint |
| $\lambda_{\text{K}}$ | $17.4\mathrm{mV}$ | K activation slope |
| $\phi$ | $0.067\,(\mathrm{ms})^{-1}$ | K rate constant |
| $E_{\text{Ca}}$ | $120\mathrm{mV}$ | Ca reversal potential |
| $E_{\text{K}}$ | $-84\mathrm{mV}$ | K reversal potential |
| $E_{\text{L}}$ | $-60\mathrm{mV}$ | Leak reversal potential |
| $g_{\text{L}}$ | $2\mathrm{mS}$ | Leak conductance |
| $g_{\text{Ca}}$ | $4.0\mathrm{mS}$ | Ca conductance |
| $g_{\text{K}}$ | $8.0\mathrm{mS}$ | K conductance |
| $\sigma_2$ | $50\mu\mathrm{A}$ | Noise amplitude |
| Threshold | $x > 0.95$ | Spike condition |
| Refractory | $x \geq 0$ | Reset condition |

## Synapse Dynamics

### Presynaptic Plasticity Fixed Point
$$\Omega(\mathrm{[Ca]}) = \begin{cases}
    A_{\mathrm{LTP}} & [\mathrm{Ca}] > \theta_{\mathrm{LTP, start}} \\
    A_{\mathrm{LTD}} & \theta_{\mathrm{LTD, start}} < [\mathrm{Ca}] < \theta_{\mathrm{LTD, end}} \\
    1 & \mathrm{otherwise}
\end{cases}$$

Piecewise function to model bidirectional plasticity.

### Presynaptic Synapse Weight
$$\frac{dW_{\mathrm{pre}}}{dt} = \frac{\Omega((1-\alpha_Wq_\text{CBD})[\mathrm{Ca}]) - W_{\mathrm{pre}}}{\tau_{W_{\mathrm{pre}}}}$$

Scaling factor on strength of synapse as determined by endocannabinoid binding.

### Postsynaptic Calcium Concentration
$$\frac{d[\mathrm{Ca}]}{dt} = \frac{\sigma(x_\mathrm{post}) - (\beta_\text{Ca}-\gamma_{Ca}q_\text{CBD})[\mathrm{Ca}]}{\tau_{\mathrm{Ca}}}$$

Used as a direct proxy for endocannabinoid release amount.

### Postsynaptic Voltage-to-Calcium Mapping
$$\sigma_\text{Ca}(x_\mathrm{post}) = \frac{1}{1 + e^{-(x_\mathrm{post} + 0.8) / 0.2}}$$

Maps postsynaptic voltage to a unitless concentration of calcium between 0 and 1.

### Synaptic Channel Opening Rate
$$T(x_\mathrm{pre} ) = \frac{T_{\mathrm{max}}}{1 + e^{-\frac{x_\mathrm{pre} - V_t}{K_p}}}$$

### Fraction of Open Synaptic Channels
$$\frac{du}{dt} = \alpha_u T(x_\mathrm{pre}) (1 - u) - \beta_u u$$

### Current from Neuron i (pre) to Neuron j (post)
$$I_{\mathrm{syn}}(i,j) = -G^{i,j}\_sW_{\mathrm{pre}}u(x_\mathrm{post}-E_{\mathrm{syn}})$$

### Total Current into Neuron j
$$I_\mathrm{syn,tot} = \sum_{i \to j} I_\mathrm{syn}(i,j)$$

## Synaptic Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| $V_T$ | $2\mathrm{mV}$ | Synaptic cleft opening rate half-activation |
| $K_p$ | $5\mathrm{mV}$ | Synaptic cleft opening rate slope |
| $T_{\max}$ | $1\mathrm{mM}$ | Max transmitter conc. |
| | | |
| $\alpha_{\text{exc}}$ | $1.1(\mathrm{mM\cdot ms})^{-1}$ | Exc. rise rate |
| $\beta_{\text{exc}}$ | $0.19(\mathrm{ms})^{-1}$ | Exc. decay rate |
| $E_{\text{exc}}$ | $0\mathrm{mV}$ | Exc. reversal potential |
| | | |
| $\alpha_{\text{inh}}$ | $5(\mathrm{mM\cdot ms})^{-1}$ | Inh. rise rate |
| $\beta_{\text{inh}}$ | $0.18(\mathrm{ms})^{-1}$ | Inh. decay rate |
| $E_{\text{inh}}$ | $-80\mathrm{mV}$ | Inh. reversal potential |

## Plasticity Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| $\theta_{\text{LTD,start}}$ | $0.25$ | LTD window start |
| $\theta_{\text{LTD,end}}$ | $0.5$ | LTD window end |
| $\theta_{\text{LTP,start}}$ | $0.75$ | LTP threshold |
| $A_{\text{LTD}}$ | $0.5$ | LTD amplitude |
| $A_{\text{LTP}}$ | $2$ | LTP amplitude |
| $\tau_W$ | $5\mathrm{s}$ | Presynaptic weight time constant |
| $\tau_{\text{Ca}}$ | $200\mathrm{ms}$ | Postsynaptic calcium time constant |
| $\alpha_W$ | $0.5$ | Maximum effective CB receptor activation reduction |
| $\beta_\text{Ca}$ | $1$ | Base dendritic calcium loss rate |
| $\gamma_\text{Ca}$ | $0.5$ | Maximum dendritic calcium loss reduction $^1$ |
| $q_\text{CBD}$ | $0$ | CBD concentration |

$^1$ In our model, dendritic calcium loss is a direct proxy for degradation of eCB by FAAH

## Sources
- Cui, Y., Ilya Prokin, Xu, H., Delord, B., Genet, S., Laurent Venance, & Berry, H. (2016). Endocannabinoid dynamics gate spike-timing dependent depression and potentiation. ELife, 5. https://doi.org/10.7554/elife.13185
- Naze, S., Bernard, C., & Jirsa, V. K. (2015). Computational Modeling of Seizure Dynamics Using Coupled Neuronal Networks: Factors Shaping Epileptiform Activity. PLOS Computational Biology, 11(5), e1004209–e1004209. https://doi.org/10.1371/journal.pcbi.1004209
