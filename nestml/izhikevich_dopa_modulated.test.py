import nest
import numpy as np
import pylab

# NOTE: setting theta currently only available via C++ hacking:
# Find the functions `set_status` (2 blocks to add) and `get_status` (1 block to add) in <modelname>.h; add entries for theta

nest.SetKernelStatus({"resolution":0.1})
simulation_time = 100000.0

nest.ResetKernel()

nest.Install("models")
model = 'izhikevich_dopa_modulated'
#model = 'izhikevich_dopa_modulated_immediate'

D1_cell_params = {
	'Kappa': 0.0289,
	'Lambda': 0.331,
	'alpha': 0.,
	'beta1': 6.3,
	'beta2': 0.,
	'nmda_ratio': 0.5,
	'V_gaba': -80.,
	'theta': 0.5
}

D2_cell_params = {
	'Kappa': 0.,
	'Lambda': 0.,
	'alpha': 0.032,
	'beta1': 0.,
	'beta2': 0.215,
	'nmda_ratio': 0.5,
	'V_gaba': -80.,
	'theta': 0.5
}

n_probe_spikes = 51 * int(simulation_time/1000)
probe_spikes = [round(i) for i in np.linspace(1, int(simulation_time)+1, n_probe_spikes)]

probe_exc_spikegen_params = {
	'spike_times': probe_spikes[::2],
	'spike_weights': [1.]*int(np.ceil(n_probe_spikes/2.))
}

probe_inh_spikegen_params = {
	'spike_times': probe_spikes[1::2],
	'spike_weights': [-1.]*int(np.floor(n_probe_spikes/2.))
}

dopa_spikes = [39., 59.] \
	+ np.linspace(62,80,10).tolist() \
	+ np.linspace(201,220,20).tolist() \
	+ np.linspace(501,600,100).tolist()

dopa_spikegen_params = {
	'spike_times': dopa_spikes,
	'spike_weights': np.ones(len(dopa_spikes)) * 45
}

syn_exc = {'receptor_type': 1}
syn_inh = {'receptor_type': 2}
syn_dopa = {'receptor_type': 3}

D1_neuron = nest.Create( model, 1, params = D1_cell_params )
D1_base = nest.Create( model, 1, params = D1_cell_params )
D2_neuron = nest.Create( model, 1, params = D2_cell_params )
D2_base = nest.Create( model, 1, params = D2_cell_params )

probe_exc_spikegen = nest.Create('spike_generator', 1, params = probe_exc_spikegen_params)
probe_inh_spikegen = nest.Create('spike_generator', 1, params = probe_inh_spikegen_params)
nest.Connect(probe_exc_spikegen, D1_neuron, syn_spec = syn_exc)
nest.Connect(probe_inh_spikegen, D1_neuron, syn_spec = syn_inh)
nest.Connect(probe_exc_spikegen, D2_neuron, syn_spec = syn_exc)
nest.Connect(probe_inh_spikegen, D2_neuron, syn_spec = syn_inh)

dopa_spikegen = nest.Create('spike_generator', 1, params = dopa_spikegen_params)
nest.Connect(dopa_spikegen, D1_neuron + D1_base + D2_neuron + D2_base, syn_spec = syn_dopa)

multimeter = nest.Create( "multimeter" )
nest.SetStatus( multimeter, { "withtime" : True, "record_from" : ["V_m"] })
nest.Connect( multimeter, D1_neuron + D1_base + D2_neuron + D2_base )
meter_shape = (int(simulation_time)-1, 4)

nest.Simulate( simulation_time )


multimeter_events = nest.GetStatus( multimeter )[0]['events']
voltages = np.reshape(multimeter_events["V_m"], meter_shape)
times = np.reshape(multimeter_events["times"], meter_shape).T[0]
voltages_nrnD1 = voltages.T[0]
voltages_baseD1 = voltages.T[1]
voltages_nrnD2 = voltages.T[2]
voltages_baseD2 = voltages.T[3]
epspsD1 = voltages_nrnD1 - voltages_baseD1
epspsD2 = voltages_nrnD2 - voltages_baseD2

pylab.subplot(2,1,1)
pylab.plot( times, epspsD1 )

pylab.subplot(2,1,2)
pylab.plot( times, voltages_baseD1 )


pylab.figure()
pylab.subplot(2,1,1)
pylab.plot( times, epspsD2 )

pylab.subplot(2,1,2)
pylab.plot( times, voltages_baseD2 )

pylab.show()
