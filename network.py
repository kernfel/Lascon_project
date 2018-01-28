import nest

def create_populations(cell_params):
	nest.Install('models')
	populations = {}
	for pop_name in cell_params:
		p = cell_params[pop_name]
		populations[pop_name] = nest.Create(p['model'], p['n'], params = p['params'])
	return populations

def create_network(pop, w = 1, wdopa = 150):
	nest.Connect(
		pop['Thal'], pop['Pyr'],
		syn_spec = {'weight': 0.05*w, 'delay': 5.6},
		conn_spec = {'rule': 'one_to_one'}
	)
	nest.Connect(
		pop['FSI'], pop['Pyr'],
		syn_spec = {'weight': -.01*w, 'delay': 0.1},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 4}
	)
	PyrInput = nest.Create('poisson_generator')
	nest.SetStatus(PyrInput, {'rate': 120.})
	nest.Connect(PyrInput, pop['Pyr'], syn_spec={'weight': 20.})

	nest.Connect(
		pop['Pyr'], pop['FSI'],
		syn_spec = {'weight': .03*w, 'delay': 0.1},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 4}
	)

	nest.Connect(
		pop['Pyr'], pop['MSN_D1'],
		syn_spec = {'weight': .07*w, 'delay': 5.1, 'receptor_type':
			nest.GetStatus(pop['MSN_D1'])[0]['receptor_types']['SPIKESGLU']},
		conn_spec = {'rule': 'one_to_one'}
	)
	nest.Connect(
		pop['MSN_D1'], pop['MSN_D1'],
		syn_spec = {'weight': -0.033*w, 'delay': 0.1, 'receptor_type':
			nest.GetStatus(pop['MSN_D1'])[0]['receptor_types']['SPIKESGABA']},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 3, 'autapses': False},
	)

	nest.Connect(
		pop['Pyr'], pop['MSN_D2'],
		syn_spec = {'weight': .07*w, 'delay': 5.1, 'receptor_type':
			nest.GetStatus(pop['MSN_D2'])[0]['receptor_types']['SPIKESGLU']},
		conn_spec = {'rule': 'one_to_one'}
	)
	nest.Connect(
		pop['MSN_D2'], pop['MSN_D2'],
		syn_spec = {'weight': -0.025*w, 'delay': 0.1, 'receptor_type':
			nest.GetStatus(pop['MSN_D2'])[0]['receptor_types']['SPIKESGABA']},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 3, 'autapses': False}
	)

	nest.Connect(
		pop['GPe'], pop['STN'],
		syn_spec = {'weight': -.04*w, 'delay': 4.0},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)
	nest.Connect(
		pop['Pyr'], pop['STN'],
		syn_spec = {'weight': .04*w, 'delay': 5.9},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)

	nest.Connect(
		pop['GPe'], pop['GPe'],
		syn_spec = {'weight': -.03*w, 'delay': 0.1},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)
	nest.Connect(
		pop['MSN_D2'], pop['GPe'],
		syn_spec = {'weight': -.03*w, 'delay': 5.0},
		conn_spec = {'rule': 'all_to_all'}
	)
	nest.Connect(
		pop['STN'], pop['GPe'][::2],
		syn_spec = {'weight': .04*w, 'delay': 2.0},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)

	nest.Connect(
		pop['GPe'], pop['GPi'],
		syn_spec = {'weight': -.03*w, 'delay': 3.0},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)
	nest.Connect(
		pop['MSN_D1'], pop['GPi'],
		syn_spec = {'weight': -.03*w, 'delay': 4.0},
		conn_spec = {'rule': 'all_to_all'}
	)
	nest.Connect(
		pop['STN'], pop['GPi'][::2],
		syn_spec = {'weight': .04*w, 'delay': 1.5},
		conn_spec = {'rule': 'fixed_indegree', 'indegree': 2}
	)

	nest.Connect(
		pop['GPi'], pop['Thal'],
		syn_spec = {'weight': -0.01*w, 'delay': 5.0},
		conn_spec = {'rule': 'one_to_one'}
	)
	ThalInput = nest.Create('poisson_generator')
	nest.SetStatus(ThalInput, {'rate': 120.})
	nest.Connect(ThalInput, pop['Thal'], syn_spec={'weight': 100.})

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
