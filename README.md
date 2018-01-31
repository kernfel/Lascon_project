# LASCON project:
## The acute effects of dopamine on the cortico-basal ganglia-thalamic circuit
### Augusto T. Figueiredo, Felix B. Kern
### VII LASCON, 7 January - 3 February 2018, São Paulo, Brasil

To run the reference model as we did, open Matlab in the KumaraveluEtAl2016 subdirectory and run `simulate_network_model(0, n, 0, 1)` with n = 0 (healthy) or n = 1 (Parkinson's)

To run our implementation of the circuit, which includes a population of dopaminergic SNc neurons as an input to the striatal neurons, run `python simulate.py`.

To see the effects of dopamine modelled in the striatal medium spiny neurons, run `python dopa_demo.py`.
