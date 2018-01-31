import pylab as pl
import nest
import numpy as np

def plot_raster(pop, spikes, simtime, showRate=True):
	k = 0;
	rates = {}
	for key in pop:
		k = k+1
		senders = nest.GetStatus(spikes[key], 'events')[0]['senders'] - pop[key][0]
		times = nest.GetStatus(spikes[key], 'events')[0]['times']
		pl.subplot(len(pop), 1, k)
		pl.scatter(times, senders, marker='.')
		pl.xlim([0,simtime])
		pl.ylim([0,len(pop[key])])
		pl.ylabel(key)
		pl.yticks([])

		if k < len(pop):		
			pl.xticks([])
		else:
			pl.xlabel('Time (ms)')

		rates[key] = nest.GetStatus(spikes[key], 'n_events')[0] / (simtime/1000. * len(pop[key]))
		if showRate:
			print key, rates[key]

	return rates

def plot_theta(pop, thetameter):
	nev = nest.GetStatus(thetameter, 'n_events')[0]
	nD1, nD2 = len(pop['MSN_D1']), len(pop['MSN_D2'])
	shape = (nev/(nD1+nD2), nD1+nD2)
	events = nest.GetStatus(thetameter, 'events')[0]
	theta = np.reshape(events['theta'], shape)
	time = np.reshape(events['times'], shape)
	pl.plot(time.T[0], theta.T[0], 'b', label='D1')
	pl.plot(time.T[1:nD1].T, theta.T[1:nD1].T, 'b')
	pl.plot(time.T[-1].T, theta.T[-1].T, 'r', label='D2')
	pl.plot(time.T[nD1:-1].T, theta.T[nD1:-1].T, 'r')
	pl.legend()
	pl.ylabel('theta')
	pl.xlabel('Time (ms)')

def setup_recordings(pop):
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
	nest.Connect(thetameter, pop['MSN_D1'] + pop['MSN_D2'])

	return spikes, voltages, thetameter

def plot_rate_bars(rates, titles, colours):
	nGroups, nBars = len(rates), len(rates[0])
	k = 0
	w = 1. / (nGroups+1)
	rects = []
	ind = np.arange(nBars)
	for ratedict in rates:
		rects.append( pl.bar(ind + w*k, ratedict.values(), w, color = colours[k]) )
		k += 1
	pl.legend(rects, titles)
	pl.xticks(ind + w*(nGroups-1)/2., rates[0].keys())
	pl.ylabel('Spike rate (Hz)')
