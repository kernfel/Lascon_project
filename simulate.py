import nest
import network
import nodes
import nest.raster_plot
import nest.voltage_trace
import pylab as pl

# simulation parameters
simtime = 20000.            # simulation time (ms)
dt = 0.1                   # simulation resolution (ms)

nest.ResetKernel()
nest.SetKernelStatus({
    'resolution': dt,    # set simulation resolution
    'print_time': True   # enable printing of simulation progress
})


cell_params = nodes.cell_params
#cell_params['MSN_D1']['params']['theta'] = cell_params['MSN_D2']['theta'] = 0.8
#for pa in cell_params.values():
#	pa['n'] = 50
pop = network.create_populations(cell_params)
network.create_network(pop, w = 1)

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

#thetameter = nest.Create('multimeter')
#nest.SetStatus(thetameter, {'withtime': True, 'record_from': ['theta']})
#nest.Connect(thetameter, pop['MSN_D1'])



nest.Simulate(simtime)

#for key in voltages:
#	nest.voltage_trace.from_device(voltages[key])
#pl.figure()

k = 0;
for key in pop:
	k = k+1
	senders = nest.GetStatus(spikes[key], 'events')[0]['senders'] - pop[key][0]
	times = nest.GetStatus(spikes[key], 'events')[0]['times']
	pl.subplot(len(pop), 1, k)
	pl.scatter(times, senders, marker='.')
	pl.xlim([0,simtime])
	pl.ylim([0,len(pop[key])])
	pl.ylabel(key)

	rate = nest.GetStatus(spikes[key], 'n_events')[0] / (simtime/1000. * len(pop[key]))
	print key, rate
#	if rate > 0 and key in rasters:
#		nest.raster_plot.from_device(spikes[key], title=key)

#thetaEvents = nest.GetStatus(thetameter)[0]['events']
#pl.plot(thetaEvents['times'], thetaEvents['theta'])

pl.show()
