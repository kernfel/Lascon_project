import nest
from numpy.random import uniform, normal, choice

def create_populations(cell_params, scale = 1):
	populations = {}
	for pop_name in cell_params:
		p = cell_params[pop_name]
		populations[pop_name] = nest.Create(p['model'], scale * p['n'], params = p['params'])

	for i in populations['SNc']:
		nest.SetStatus([i], {'V_m': uniform(-100.,-40.), 'U_m': uniform(-100.,100.)})

	return populations

def create_network(pop):
	nest.Connect(
		pop['Thal'], pop['Pyr'],
		syn_spec = {'weight': 5., 'delay': 5.6},
		conn_spec = {'rule': 'one_to_one'}
	)
	nest.Connect(
		pop['Inh'], pop['Pyr'],
		syn_spec = {'weight': -18., 'delay': 0.1},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 4}
	)
	PyrInput = nest.Create('poisson_generator')
	nest.SetStatus(PyrInput, {'rate': 20.})
	nest.Connect(PyrInput, pop['Pyr'], syn_spec={'weight': 3.})

	nest.Connect(
		pop['Pyr'], pop['Inh'],
		syn_spec = {'weight': 3.3, 'delay': 0.1},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 4}
	)

	nest.Connect(
		pop['Pyr'], pop['MSN_D1'],
		syn_spec = {'weight': 43., 'delay': 5.1, 'receptor_type':
			nest.GetStatus(pop['MSN_D1'])[0]['receptor_types']['SPIKESGLU']},
		conn_spec = {'rule': 'one_to_one'}
	)
	nest.Connect(
		pop['MSN_D1'], pop['MSN_D1'],
		syn_spec = {'weight': -8., 'delay': 0.1, 'receptor_type':
			nest.GetStatus(pop['MSN_D1'])[0]['receptor_types']['SPIKESGABA']},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 3, 'autapses': False},
	)

	nest.Connect(
		pop['Pyr'], pop['MSN_D2'],
		syn_spec = {'weight': 43., 'delay': 5.1, 'receptor_type':
			nest.GetStatus(pop['MSN_D2'])[0]['receptor_types']['SPIKESGLU']},
		conn_spec = {'rule': 'one_to_one'}
	)
	nest.Connect(
		pop['MSN_D2'], pop['MSN_D2'],
		syn_spec = {'weight': -8., 'delay': 0.1, 'receptor_type':
			nest.GetStatus(pop['MSN_D2'])[0]['receptor_types']['SPIKESGABA']},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 3, 'autapses': False}
	)

	nest.Connect(
		pop['GPe'], pop['STN'],
		syn_spec = {'weight': -4., 'delay': 4.0},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)
	nest.Connect(
		pop['Pyr'], pop['STN'],
		syn_spec = {'weight': 7., 'delay': 5.9},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)

	nest.Connect(
		pop['GPe'], pop['GPe'],
		syn_spec = {'weight': -2., 'delay': 0.1},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)
	nest.Connect(
		pop['MSN_D2'], pop['GPe'],
		syn_spec = {'weight': -50. / len(pop['MSN_D2']), 'delay': 5.0},
		conn_spec = {'rule': 'all_to_all'}
	)
	nest.Connect(
		pop['STN'], pop['GPe'][::2],
		syn_spec = {'weight': 5., 'delay': 2.0},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)

	nest.Connect(
		pop['GPe'], pop['GPi'],
		syn_spec = {'weight': -5., 'delay': 3.0},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)
	nest.Connect(
		pop['MSN_D1'], pop['GPi'],
		syn_spec = {'weight': -50. / len(pop['MSN_D1']), 'delay': 4.0},
		conn_spec = {'rule': 'all_to_all'}
	)
	nest.Connect(
		pop['STN'], pop['GPi'][::2],
		syn_spec = {'weight': 5., 'delay': 1.5},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)

	nest.Connect(
		pop['GPi'], pop['Thal'],
		syn_spec = {'weight': -25., 'delay': 5.0},
		conn_spec = {'rule': 'one_to_one'}
	)

def connect_SNc(pop, frac = 1., weight = 12, outdeg = 0.6):
	if frac >= 1.:
		# Use the whole population
		sn = pop['SNc']
	else:
		# Pick a fractional part of the population
		n = int(round(len(pop['SNc']) * frac))
		idx = choice(len(pop['SNc']), n, replace=False)
		sn = [pop['SNc'][i] for i in idx]
	
	nest.Connect(
		sn, pop['MSN_D1'],
		syn_spec = {'weight': 100.0 * weight, 'delay': 3.0, 'receptor_type':
			nest.GetStatus(pop['MSN_D1'])[0]['receptor_types']['SPIKESDOPA']},
		conn_spec = {'rule': 'fixed_outdegree', 'outdegree': int(round(outdeg*len(pop['MSN_D1'])))}
	)
	nest.Connect(
		sn, pop['MSN_D2'],
		syn_spec = {'weight': 100.0 * weight, 'delay': 3.0, 'receptor_type':
			nest.GetStatus(pop['MSN_D2'])[0]['receptor_types']['SPIKESDOPA']},
		conn_spec = {'rule': 'fixed_outdegree', 'outdegree': int(round(outdeg*len(pop['MSN_D2'])))}
	)

def add_stims(pop, times, amp = 100):
	''' Creates a poisson input for each (tStart, tEnd) in times and connects it to the given population. '''
	for start, stop in times:
		stim = nest.Create('dc_generator')
		nest.SetStatus(stim, {'start': float(start), 'stop': float(stop), 'amplitude': float(amp)})
		nest.Connect(stim, pop)


if __name__ == '__main__':
	pop = create_populations()
	create_network(pop)

