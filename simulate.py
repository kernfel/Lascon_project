import nest
import network
import nodes
import nest.raster_plot
import nest.voltage_trace
import pylab as pl
from tools import *

# simulation parameters
simtime = 60000.            # simulation time (ms)
dt = 0.1                   # simulation resolution (ms)

nest.Install('models')

def SNc_experiment(frac, theta, weight = 1.5, outdeg = 0.6, title = ''):
	nest.ResetKernel()
	nest.SetKernelStatus({
		'resolution': dt,    # set simulation resolution
		'print_time': True   # enable printing of simulation progress
	})

	cell_params = nodes.cell_params
	cell_params['MSN_D1']['params']['theta'] = theta
	cell_params['MSN_D2']['params']['theta'] = theta
	pop = network.create_populations(cell_params, scale = 1)
	network.create_network(pop)
	network.connect_SNc(pop, frac = frac, weight = weight, outdeg = 0.6)

	#stim_times = [(1000,1010)]
	#network.add_stims(pop['SNc'], stim_times, amp = 500)

	spikes, voltages, thetameter = setup_recordings(pop)

	nest.Simulate(simtime)

	pl.figure()
	pl.title(title)
	rates = plot_raster(pop, spikes, simtime)

	pl.figure()
	pl.title(title)
	plot_theta(pop, thetameter)

	return rates

if __name__ == '__main__':
	rates_healthy = SNc_experiment(1, 0.8, title = 'Healthy')
	rates_untreated = SNc_experiment(0.1, 0.2, title = 'untreated PD')
	rates_ldopa = SNc_experiment(0.1, 0.3, weight = 6., title = 'PD + levodopa')
	
	pl.show()
