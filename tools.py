import pylab as pl
import nest

def raster(pop, spikes, simtime, showRate=True):
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

		if showRate:
			rate = nest.GetStatus(spikes[key], 'n_events')[0] / (simtime/1000. * len(pop[key]))
			print key, rate
