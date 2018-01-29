import nest
from numpy.random import uniform, normal

def create_populations(cell_params, scale = 1):
	nest.Install('models')
	populations = {}
	for pop_name in cell_params:
		p = cell_params[pop_name]
		populations[pop_name] = nest.Create(p['model'], scale * p['n'], params = p['params'])
#		i_e = nest.GetStatus(populations[pop_name])[0]['I_e']
#		for i in populations[pop_name]:
#			nest.SetStatus([i], {'V_m': uniform(-100., -40.), 'U_m': uniform(0.,1.), 'I_e': i_e + normal()})
	return populations

def create_network(pop, wdopa = 150):
	nest.Connect(
		pop['Thal'], pop['Pyr'],
		syn_spec = {'weight': 0.5, 'delay': 5.6},
		conn_spec = {'rule': 'one_to_one'}
	)
	nest.Connect(
		pop['Inh'], pop['Pyr'],
		syn_spec = {'weight': -10., 'delay': 0.1},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 4}
	)
	PyrInput = nest.Create('poisson_generator')
	nest.SetStatus(PyrInput, {'rate': 60.})
	nest.Connect(PyrInput, pop['Pyr'], syn_spec={'weight': 10.})

	nest.Connect(
		pop['Pyr'], pop['Inh'],
		syn_spec = {'weight': 10., 'delay': 0.1},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 4}
	)

	nest.Connect(
		pop['Pyr'], pop['MSN_D1'],
		syn_spec = {'weight': 13., 'delay': 5.1, 'receptor_type':
			nest.GetStatus(pop['MSN_D1'])[0]['receptor_types']['SPIKESGLU']},
		conn_spec = {'rule': 'one_to_one'}
	)
	nest.Connect(
		pop['MSN_D1'], pop['MSN_D1'],
		syn_spec = {'weight': -5.5, 'delay': 0.1, 'receptor_type':
			nest.GetStatus(pop['MSN_D1'])[0]['receptor_types']['SPIKESGABA']},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 3, 'autapses': False},
	)

	nest.Connect(
		pop['Pyr'], pop['MSN_D2'],
		syn_spec = {'weight': 13., 'delay': 5.1, 'receptor_type':
			nest.GetStatus(pop['MSN_D2'])[0]['receptor_types']['SPIKESGLU']},
		conn_spec = {'rule': 'one_to_one'}
	)
	nest.Connect(
		pop['MSN_D2'], pop['MSN_D2'],
		syn_spec = {'weight': -5.5, 'delay': 0.1, 'receptor_type':
			nest.GetStatus(pop['MSN_D2'])[0]['receptor_types']['SPIKESGABA']},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 3, 'autapses': False}
	)

	nest.Connect(
		pop['GPe'], pop['STN'],
		syn_spec = {'weight': -20., 'delay': 4.0},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)
	nest.Connect(
		pop['Pyr'], pop['STN'],
		syn_spec = {'weight': 10., 'delay': 5.9},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)

	nest.Connect(
		pop['GPe'], pop['GPe'],
		syn_spec = {'weight': -3., 'delay': 0.1},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)
	nest.Connect(
		pop['MSN_D2'], pop['GPe'],
		syn_spec = {'weight': -15., 'delay': 5.0},
		conn_spec = {'rule': 'all_to_all'}
	)
	nest.Connect(
		pop['STN'], pop['GPe'][::2],
		syn_spec = {'weight': 8., 'delay': 2.0},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)

	nest.Connect(
		pop['GPe'], pop['GPi'],
		syn_spec = {'weight': -15., 'delay': 3.0},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)
	nest.Connect(
		pop['MSN_D1'], pop['GPi'],
		syn_spec = {'weight': -15., 'delay': 4.0},
		conn_spec = {'rule': 'all_to_all'}
	)
	nest.Connect(
		pop['STN'], pop['GPi'][::2],
		syn_spec = {'weight': 8., 'delay': 1.5},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)

	nest.Connect(
		pop['GPi'], pop['Thal'],
		syn_spec = {'weight': -10., 'delay': 5.0},
		conn_spec = {'rule': 'one_to_one'}
	)

	nest.Connect(
		pop['SNc'], pop['MSN_D1'],
		syn_spec = {'weight': 100.0*wdopa, 'delay': 3.0, 'receptor_type':
			nest.GetStatus(pop['MSN_D1'])[0]['receptor_types']['SPIKESDOPA']},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)
	nest.Connect(
		pop['SNc'], pop['MSN_D2'],
		syn_spec = {'weight': 100.0*wdopa, 'delay': 3.0, 'receptor_type':
			nest.GetStatus(pop['MSN_D2'])[0]['receptor_types']['SPIKESDOPA']},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)


if __name__ == '__main__':
	pop = create_populations()
	create_network(pop)

