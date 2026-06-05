# This file details all the equations used in the model

## Neuron Populations

### Excitatory (Hindmarsh-Rose)
$${x}' = y - a x^3 + b x^2 - z + I_{\text{app,1}} + J (\bar{x} - x) + \sigma_1 I_\text{syn,tot} + W(t)$$

Unitless membrane potential. Evolves according to cubic derivative modulated by y and z variables. Stimulus comes from a constant current, a stochastic current, and synaptic currents. Intrapopulation gap junction coupling simulated by partially coupling neuron potential to population mean.

$$y' = c - d x^2 - y$$

Spiking variable. Models ion transport by voltage-gated sodium and potassium ion channels.

$$z' = r(s(x + \bar{V}^\* - \eta) - \bar{z})$$

Slow adaptation variable. Introduces chaotic dynamics influenced by individual excitatory neuron potential and mean inhibitory population potential.

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

## Sources
- Computational Modeling of Seizure Dynamics Using Coupled Neuronal Networks: Factors Shaping Epileptiform Activity
    - Naze et al, 2015
- Endocannabinoid dynamics gate spike-timing dependent depression and potentiation
    - Cui et al, 2016