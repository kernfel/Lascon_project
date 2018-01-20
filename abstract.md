#The acute effects of dopamine on the cortico-basal ganglia-thalamic circuit
Augusto T. Figueiredo, Physiology Department, Universidade Federal de Sergipe, São Cristóvão, Sergipe, Brazil

Felix B. Kern, Sussex Neuroscience, University of Sussex, Brighton, UK

##Abstract
In this project, we aim to explore how the presence or absence of dopaminergic substantia nigra pars compacta (SNc) neurons affects the activity of a cortex-basal ganglia-thalamus circuit. The absence of dopaminergic SNc neurons is the pathophysiological marker of Parkinson's Disease (PD), yet most models of PD focus on the disease's impact on long-term synaptic plasticity, ignoring the short-term effects of dopamine on excitability and synaptic function. To address this gap, we will use a circuit model composed of Izhikevich neurons fitted to the activity patterns of the more biophysically detailed neurons used in [1], aiming to maintain the overall properties of the circuit. To model the effects of dopamine, we will replace the striatal medium spiny neurons (MSN) with a dopamine-sensitive Izhikevich model [2] and add a population of SNc neurons modelled after [3]. Synaptic dopamine concentration at the MSN will be modelled using Michaelis-Menten kinetics [4]. The circuit will be simulated in NEST using custom models written in NESTML. We hope to see effects related to the therapeutic effects (amelioration of bradykinesia) or side effects (dyskinesia) of levodopa treatment, measured in terms of thalamocortical activity patterns as described in [5].

1. [Kumaravelu, K., Brocker, D. T., & Grill, W. M. (2016). A biophysical model of the cortex-basal ganglia-thalamus network in the 6-OHDA lesioned rat model of Parkinson’s disease. Journal of computational neuroscience, 40(2), 207-229.](http://dx.doi.org/10.1007/s10827-016-0593-9)
2. [Humphries, M. D., Lepora, N., Wood, R., & Gurney, K. (2009). Capturing Dopaminergic Modulation and Bimodal Membrane Behaviour of Striatal Medium Spiny Neurons in Accurate, Reduced Models. Frontiers in Computational Neuroscience, 3, 26.](http://doi.org/10.3389/neuro.10.026.2009)
3. [Kuznetsova, A. Y., Huertas, M. A., Kuznetsov, A. S., Paladini, C. A., & Canavier, C. C. (2010). Regulation of firing frequency in a computational model of a midbrain dopaminergic neuron. Journal of computational neuroscience, 28(3), 389-403.](https://doi.org/10.1007/s10827-010-0222-y)
4. [Kawagoe, K. T., Garris, P. A., Wiedemann, D. J., & Wightman, R. M. (1992). Regulation of transient dopamine concentration gradients in the microenvironment surrounding nerve terminals in the rat striatum. Neuroscience, 51(1), 55-64.](https://doi.org/10.1016/0306-4522(92)90470-M)
5. [Kerr, C. C., Van Albada, S. J., Neymotin, S. A., Chadderdon, G. L., Robinson, P. A., & Lytton, W. W. (2013). Cortical information flow in Parkinson’s disease: a composite network/field model. Frontiers in Computational Neuroscience, 7, 39.](https://dx.doi.org/10.3389%2Ffncom.2013.00039)