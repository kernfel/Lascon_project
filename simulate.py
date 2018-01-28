import nest
import network
import nodes
import nest.raster_plot
import nest.voltage_trace
import pylab as pl
from tools import *

# simulation parameters
simtime = 2000.            # simulation time (ms)
dt = 0.1                   # simulation resolution (ms)

nest.ResetKernel()
nest.SetKernelStatus({
    'resolution': dt,    # set simulation resolution
    'print_time': True   # enable printing of simulation progress
})


cell_params = nodes.cell_params
cell_params['MSN_D1']['params']['theta'] = cell_params['MSN_D2']['theta'] = 0.8
pop = network.create_populations(cell_params, scale = 1)
network.create_network(pop)

spikes = {}
voltages = {}
for key in pop:
	spikes[key] = nest.Create('spike_detector')
	nest.SetStatus(spikes[key], [{'withtime': True,
                          'withgid': True,
                          'to_file': False}])
	nest.Connect(pop[key], spikes[key])

	voltages[key] = nest.Create('voltmeter')
	nest.SetStatus(voltages[key], {'withgid': True, 'withtime': True})
	nest.Connect(voltages[key], pop[key])

thetameter = nest.Create('multimeter')
nest.SetStatus(thetameter, {'withtime': True, 'record_from': ['theta']})
nest.Connect(thetameter, [pop['MSN_D1'][0]])



nest.Simulate(simtime)

#for key in voltages:
#	nest.voltage_trace.from_device(voltages[key])
#pl.figure()

raster(pop, spikes, simtime)

#pl.figure()
#thetaEvents = nest.GetStatus(thetameter)[0]['events']
#pl.plot(thetaEvents['times'], thetaEvents['theta'])

pl.show()
