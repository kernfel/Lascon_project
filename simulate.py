import nest
import network
import nodes
import nest.raster_plot
import nest.voltage_trace
import pylab as pl
from tools import *

# simulation parameters
simtime = 20000.            # simulation time (ms)
dt = 0.1                   # simulation resolution (ms)

nest.ResetKernel()
nest.SetKernelStatus({
    'resolution': dt,    # set simulation resolution
    'print_time': True   # enable printing of simulation progress
})


cell_params = nodes.cell_params
cell_params['MSN_D1']['params']['theta'] = 0.8
cell_params['MSN_D2']['params']['theta'] = 0.8
pop = network.create_populations(cell_params, scale = 1)
network.create_network(pop)
network.connect_SNc(pop, frac = 1, weight = 1.5, outdeg = 0.6)

stim_times = [(10000,15000)]
network.add_stims(pop['SNc'], stim_times, rate = 100, weight = 20)

spikes, voltages, thetameter = setup_recordings(pop)



nest.Simulate(simtime)

plot_raster(pop, spikes, simtime)

pl.figure()
plot_theta(pop, thetameter)

pl.show()
