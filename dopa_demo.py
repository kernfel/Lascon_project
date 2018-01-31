import nest
import numpy as np
import pylab as pl
from nodes import cell_params

nest.SetKernelStatus({"resolution":0.1})
simtime = 80000.0

chunks = [10000., 20000., 30000., 40000., simtime]
weights = [2500., 5000., 7500., 10000., 0.]

nest.ResetKernel()

nest.Install("models")

# Create MSNs
D1_neuron = nest.Create(cell_params['MSN_D1']['model'], 1, params = cell_params['MSN_D1']['params'])
D2_neuron = nest.Create(cell_params['MSN_D2']['model'], 1, params = cell_params['MSN_D2']['params'])
nest.SetStatus(D1_neuron + D2_neuron, {'theta': 0., 'dopa': 0., 'V_m': -80.})

# Create an SNc neuron
SNc = nest.Create(cell_params['SNc']['model'], 1, params = cell_params['SNc']['params'])
nest.SetStatus(SNc, {'I_e': 10.})

# Create a stimulator
probe_spikes = np.arange(10, simtime, 1000) # 1 spike / second
probe = nest.Create('spike_generator', 1,
	params = {'spike_times': probe_spikes, 'spike_weights': np.ones(probe_spikes.shape)})


# Connections
nest.Connect(probe, D1_neuron + D2_neuron,
	syn_spec = {'receptor_type': nest.GetStatus(D1_neuron, 'receptor_types')[0]['SPIKESGLU']})

nest.Connect(SNc, D1_neuron + D2_neuron,
	syn_spec = {'receptor_type': nest.GetStatus(D1_neuron, 'receptor_types')[0]['SPIKESDOPA'],
				'weight': 0.})
dopa_syn = nest.GetConnections(source=SNc)


# Recordings
meter = nest.Create('multimeter')
nest.SetStatus(meter, {"withtime" : True, "record_from" : ['dopa', 'theta', 'V_m']})
nest.Connect(meter, D1_neuron + D2_neuron)
meter_shape = (int(simtime)-1, 2)


# Run in chunks
t = 0.
for weight, tNext in zip(weights, chunks):
	nest.SetStatus(dopa_syn, {'weight': weight})
	print t, tNext, weight
	nest.Simulate(tNext - t)
	t = tNext


# Get data from the meter
events = nest.GetStatus(meter, 'events')[0]
times = np.reshape(events["times"], meter_shape).T[0]
dopa = np.reshape(events["dopa"], meter_shape)
theta = np.reshape(events["theta"], meter_shape)
V_m = np.reshape(events["V_m"], meter_shape)

# Plot it!
pl.subplot(3,1,1)
pl.plot(times, V_m.T[0], 'g')
pl.ylabel('V_m')

pl.subplot(3,1,2)
pl.plot(times, V_m.T[1], 'r')
pl.hlines(max(V_m.T[1]), 0, times[-1], linestyles='dashed')
pl.ylim([-80.1, -79.2])
pl.ylabel('V_m')

ax1 = pl.subplot(3,1,3)
pl.xlabel('Time (ms)')
ax1.plot(times, theta.T[0])
ax1.set_ylabel('theta')

ax2 = ax1.twinx()
ax2.plot(times, dopa.T[0], 'black')
ax2.set_ylabel('Synaptic dopamine (nM)')


# Plot EPSP overlays
pl.figure()
tObs = 1000, 20000, 40000
for i in 1, 2:
	pl.subplot(2,1,i)
	pl.title('D%d' % i)
	pl.ylabel('V_m')
	V0 = V_m.T[i-1][tObs[0]]
	for t in tObs:
		pl.plot(times[t:t+40] - t, V_m.T[i-1][t:t+40] - (V_m.T[i-1][t] - V0), label = '%d ms' % t)
pl.xlabel('Time (ms)')
pl.legend()


pl.show()
