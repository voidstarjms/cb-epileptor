# This file details all the equations and parameters used in the model

Values listed are the defaults in `parameters/params.yaml`. Where a parameter is swept, the sweep range is given in the description.

## Global Simulation Control
| Parameter | Value | Description |
|-----------|-------|-------------|
| `SIM_DURATION` | $40\mathrm{s}$ | Total simulated time (recorded, after transient) |
| `NUM_CELLS` | $5$ | Cells per population |
| `TAU_CLOCK` | $1\mathrm{msec}$ | Base clock step ($\tau$) |
| `DT_SCALING` | $20$ | Integration step is `TAU_CLOCK / DT_SCALING` |
| `TRANSIENT` | $10\mathrm{s}$ | Discarded transient, run before recording |
| `NOISE_INIT_OFFSET` | $1.5$ | Offset of initial excitatory $x$ above $\eta$ |

## Coupling and Global Logic
| Parameter | Value | Description |
|-----------|-------|-------------|
| $\iota$ (`ISOLATE`) | $1$ | $0$: decoupled, $1$: coupled $^1$ |
| $J$ | $0.1$ | Intrapopulation coupling gain (swept over $[0,\,0.2]$) |
| $W_\text{max}$ | $0.012$ | Noise amplitude |

$^1$ $\iota$ gates gap junction coupling, noise, synaptic input, and the interpopulation mean-field terms. Setting $\iota = 0$ leaves each neuron fully isolated and noiseless.

Noise is Gaussian white noise, $W(t) = W_\text{max}\,\xi(t)$.

## Time-Varying Schedules
By default, $\eta$, $J$, $G_\text{intra}$, and $G_\text{inter}$ are held constant at their scalar values. With `--timed-array` (`-t`), each is instead stepped through its `*_VALS` list, advancing one entry every `*_DT`, and holding the last value once the list is exhausted.

| Parameter | Schedule (`*_VALS`) | Step (`*_DT`) |
|-----------|---------------------|---------------|
| $\eta$ | `BASE_EXCITE_VALS`: $[-2.5]$ | $20\mathrm{s}$ |
| $J$ | `COUPLING_VALS`: $[0.2]$ | $20\mathrm{s}$ |
| $G_\text{inter}$ | `G_INTER_VALS`: $[1]\mu\mathrm{S}$ | $20\mathrm{s}$ |
| $G_\text{intra}$ | `G_INTRA_VALS`: $[0.5]\mu\mathrm{S}$ | $20\mathrm{s}$ |

Schedule time starts at $t = 0$, including the transient.

## Initial Conditions
| Variable | Initial value |
|----------|---------------|
| $x$ (excitatory) | $\eta(0) + $ `NOISE_INIT_OFFSET` $ + W_\text{max}\,\mathcal{N}(0,1)$ |
| $y$ | $c - d x^2$ |
| $z$ | $s(x - \eta(0))$ |
| $V$ (inhibitory) | $E_L + W_\text{max}\,\mathcal{N}(0,1)\,\mathrm{V}$ |
| $n$ | $n_\infty(V)$ |
| $W_\text{pre}$ | $1$ |
| $[\mathrm{Ca}]$ | $\frac{\theta_\text{LTD,start}}{2}(1 + \mathcal{N}(0,1))$ |
| $u$ | $0$ |

## Neuron Populations

### Excitatory (Hindmarsh-Rose)
$$\tau x' = y - a x^3 + b x^2 - z + I_{\text{app,1}} + \iota\left[J (\bar{x} - x) + W(t) + \sigma_1 \frac{I_\text{syn,tot}}{\psi}\right]$$

Unitless membrane potential. Evolves according to cubic derivative modulated by y and z variables. Stimulus comes from a constant current, a stochastic current, and synaptic currents. Intrapopulation gap junction coupling simulated by partially coupling neuron potential to population mean.

$$\tau y' = c - d x^2 - y$$

Spiking variable. Models ion transport by voltage-gated sodium and potassium ion channels.

$$z' = r(s(x + \iota\bar{V}^\* - \eta) - \bar{z})$$

Slow adaptation variable. Introduces chaotic dynamics influenced by individual excitatory neuron potential and mean inhibitory population potential.

$\bar{x}$ and $\bar{z}$ are excitatory population means; $\bar{V}^\*$ is the inhibitory population mean of $V^\*$.

### Excitatory Neuron Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| $a$ | $1.0$ | x nullcline cubic scalar |
| $b$ | $3.0$ | x nullcline quadratic scalar |
| $c$ | $1.0$ | y nullcline constant element |
| $d$ | $5.0$ | y nullcline quadratic scalar |
| $s$ | $4.0$ | z excitability scalar |
| $I_{\text{app,1}}$ | $6$ | Applied current |
| $\eta$ | $-3.5$ | Resting potential (swept over $[-4.5,\,-2.5]$) |
| $r$ | $2\times10^{-5}(\mathrm{ms})^{-1}$ | Slow adaptation rate |
| $\sigma_1$ | $0.02$ | Synaptic input gain |
| $\psi$ | $0.05\mu\mathrm{A}$ | Current scaling |
| Threshold | $x > 1.5$ | Spike condition |
| Refractory | $x \geq 0$ | Reset condition |

### Inhibitory (Morris-Lecar)
$$C_mV' = I_{\text{app,2}} - g_L(V - E_L) - g_K n(V - E_K) - g_{Ca} m_\infty(V)(V - E_{Ca})
    + \iota\left[\sigma_2\left(W(t) + J(\bar{V}^\* - V^{\*}) - \kappa_z(\bar{z} - z_0)\right) + I_\text{syn,tot}\right]$$

Membrane voltage. Evolves according to ion channel currents similar to the Hodgkin-Huxley model. Intrapopulation gap junction coupling simulated by partially coupling neuron potential to population mean. Input from excitatory population slow variable $\bar{z}$ drives slow dynamics.

$$n' = \phi\frac{n_\infty(V) - n}{\tau_n(V)}$$

Fraction of open potassium channels.

$$m_\infty(V) = \frac{1}{2}[1 + \tanh((V-h_\text{Ca})/\lambda_\text{Ca})]$$

Fraction of open calcium channels.

$$n_\infty(V) = \frac{1}{2}[1 + \tanh((V-h_\text{K})/\lambda_\text{K})]$$

Steady-state fraction of open potassium channels.

$$\tau_n(V) = \frac{1}{\cosh((V-h_\text{K})/ (2 \lambda_\text{K}))}$$

Potassium channel opening time scale.

$$V^{*} = \frac{V}{20\mathrm{mV}}$$

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
| $\phi$ | $0.067(\mathrm{ms})^{-1}$ | K rate constant |
| $E_{\text{Ca}}$ | $120\mathrm{mV}$ | Ca reversal potential |
| $E_{\text{K}}$ | $-84\mathrm{mV}$ | K reversal potential |
| $E_{\text{L}}$ | $-60\mathrm{mV}$ | Leak reversal potential |
| $g_{\text{L}}$ | $2\mathrm{mS}$ | Leak conductance |
| $g_{\text{Ca}}$ | $4.0\mathrm{mS}$ | Ca conductance |
| $g_{\text{K}}$ | $8.0\mathrm{mS}$ | K conductance |
| $\sigma_2$ | $50\mu\mathrm{A}$ | Noise and coupling scaling |
| $\kappa_z$ | $0.15$ | Excitatory slow variable drive gain (`IPOP_Z_BAR_SCALE`) |
| $z_0$ | $6$ | Excitatory slow variable drive offset (`IPOP_Z_BAR_OFFSET`) |
| Threshold | $V^\* > 0.95$ | Spike condition |
| Refractory | $V^\* \geq 0$ | Reset condition |

## Synapse Dynamics

Every synapse group (E→E, E→I, I→E, I→I) uses the same equations. Each ordered pair of neurons, including self-connections within a population, is connected with probability $p$.

### Presynaptic Plasticity Fixed Point
$$\Omega(\mathrm{[Ca]}) = \begin{cases}
    1 + A_{\mathrm{LTP}} & [\mathrm{Ca}] > \theta_{\mathrm{LTP, start}} \\
    1 - A_{\mathrm{LTD}} & \theta_{\mathrm{LTD, start}} < [\mathrm{Ca}] < \theta_{\mathrm{LTD, end}} \\
    1 & \mathrm{otherwise}
\end{cases}$$

Piecewise function to model bidirectional plasticity.

### Presynaptic Synapse Weight
$$\frac{dW_{\mathrm{pre}}}{dt} = \frac{\Omega(\rho_\text{CB}(1-\alpha_Wq_\text{CBD})[\mathrm{Ca}]) - W_{\mathrm{pre}}}{\tau_{W_{\mathrm{pre}}}}$$

Scaling factor on strength of synapse as determined by endocannabinoid binding. $\rho_\text{CB}$ is the CB receptor prevalence of the presynaptic population. When CB synapses are disabled (`--cb` not passed), $W_\text{pre}$ still evolves but is held at $1$ in the synaptic current.

### Postsynaptic Calcium Concentration $^2$
$$\frac{d[\mathrm{Ca}]}{dt} = (1 - [\mathrm{Ca}])\frac{\sigma_\text{Ca}(x_\mathrm{post}) - (\beta_\text{Ca}-\gamma_{Ca}q_\text{CBD})[\mathrm{Ca}]}{\tau_{\mathrm{Ca}}}$$

Used as a direct proxy for endocannabinoid release amount. The $(1 - [\mathrm{Ca}])$ factor keeps $[\mathrm{Ca}] \leq 1$.

### Postsynaptic Voltage-to-Calcium Mapping $^2$
$$\sigma_\text{Ca}(x_\mathrm{post}) = \frac{1}{1 + e^{-(x_\mathrm{post} + \theta_\sigma) / k_\sigma}}$$

Maps postsynaptic voltage to a unitless concentration of calcium between 0 and 1.

### Synaptic Voltage Scaling $^2$
$$\tilde{x} = \frac{x}{\sigma_1}\,\mathrm{mV}$$

Converts unitless potential to millivolts for use in the transmitter and current equations ($\tilde{x} = 50x\,\mathrm{mV}$ at default $\sigma_1$).

### Synaptic Channel Opening Rate $^2$
$$T(x_\mathrm{pre} ) = \frac{T_{\mathrm{max}}}{1 + e^{-\frac{\tilde{x}_\mathrm{pre} - V_t}{K_p}}}$$

### Fraction of Open Synaptic Channels $^2$
$$\frac{du}{dt} = \alpha_u T(x_\mathrm{pre}) (1 - u) - \beta_u u$$

### Current from Neuron i (pre) to Neuron j (post) $^2$
$$I_{\mathrm{syn}}(i,j) = -G^{i,j}\_sW_{\mathrm{pre}}u(\tilde{x}_\mathrm{post}-E_{\mathrm{syn}})$$

$G^{i,j}_s$ is $G_\text{intra}$ for synapses within a population and $G_\text{inter}$ for synapses between populations. $\alpha_u$, $\beta_u$, and $E_\text{syn}$ take the excitatory or inhibitory values according to the presynaptic population.

$^2$ In these equations, $x$ refers to $V^\*$ if the neuron in question is inhibitory

### Total Current into Neuron j
$$I_\mathrm{syn,tot} = \sum_{i \to j} I_\mathrm{syn}(i,j)$$

## Synaptic Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| $V_T$ | $2\mathrm{mV}$ | Synaptic cleft opening rate half-activation |
| $K_p$ | $5\mathrm{mV}$ | Synaptic cleft opening rate slope |
| $T_{\max}$ | $1\mathrm{mM}$ | Max transmitter conc. |
| $p$ | $1$ | Connection probability (`PCT_CONNECT`) |
| | | |
| $\alpha_{\text{exc}}$ | $1.1(\mathrm{mM\cdot ms})^{-1}$ | Exc. rise rate |
| $\beta_{\text{exc}}$ | $0.19(\mathrm{ms})^{-1}$ | Exc. decay rate |
| $E_{\text{exc}}$ | $0\mathrm{mV}$ | Exc. reversal potential |
| | | |
| $\alpha_{\text{inh}}$ | $5(\mathrm{mM\cdot ms})^{-1}$ | Inh. rise rate |
| $\beta_{\text{inh}}$ | $0.18(\mathrm{ms})^{-1}$ | Inh. decay rate |
| $E_{\text{inh}}$ | $-80\mathrm{mV}$ | Inh. reversal potential |
| | | |
| $G_{\text{intra}}$ | $0.5\mu\mathrm{S}$ | Intrapopulation synapse strength (swept over $[0.5,\,2]\mu\mathrm{S}$) |
| $G_{\text{inter}}$ | $1\mu\mathrm{S}$ | Interpopulation synapse strength (swept over $[1,\,4]\mu\mathrm{S}$) |

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
| $\theta_\sigma$ | $0.8$ | Voltage-to-calcium sigmoid shift (`CA_SIGMOID_SHIFT`) |
| $k_\sigma$ | $0.2$ | Voltage-to-calcium sigmoid slope (`CA_SIGMOID_SLOPE`) |
| $\alpha_W$ | $0.5$ | Maximum effective CB receptor activation reduction |
| $\beta_\text{Ca}$ | $1$ | Base dendritic calcium loss rate |
| $\gamma_\text{Ca}$ | $0$ | Maximum dendritic calcium loss reduction $^3$ |
| $q_\text{CBD}$ | $0$ | CBD concentration |
| $\rho_\text{CB,exc}$ | $1$ | CB receptor prevalence on excitatory presynaptic terminals |
| $\rho_\text{CB,inh}$ | $1$ | CB receptor prevalence on inhibitory presynaptic terminals |

$^3$ In our model, dendritic calcium loss is a direct proxy for degradation of eCB by FAAH

## Sources
- Cui, Y., Ilya Prokin, Xu, H., Delord, B., Genet, S., Laurent Venance, & Berry, H. (2016). Endocannabinoid dynamics gate spike-timing dependent depression and potentiation. ELife, 5. https://doi.org/10.7554/elife.13185
- Naze, S., Bernard, C., & Jirsa, V. K. (2015). Computational Modeling of Seizure Dynamics Using Coupled Neuronal Networks: Factors Shaping Epileptiform Activity. PLOS Computational Biology, 11(5), e1004209–e1004209. https://doi.org/10.1371/journal.pcbi.1004209
