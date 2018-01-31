''' Neuron model parameters '''

cell_params = {}
cell_order = ['Pyr', 'Inh', 'MSN_D1', 'MSN_D2', 'STN', 'GPe', 'GPi', 'Thal', 'SNc']

tau_in = 5.
tau_ex = 2.

cell_params['SNc'] = {'model': 'izhikevich', 'params':
	{# Source: Hand tuned using KuznetsovaEtAl2010/izhikevich.py
		'V_m': -62.6,
		'a': .0024,
		'b': .42,
		'c': -60.,
		'd': 50.,
		'V_th': 10.
	},
	'n': 100
}

cell_params['MSN_D1'] = {'model': 'izhikevich_dopa_modulated', 'params':
	{# Source: Fountas & Shanahan 2017
		"a" :          0.01,
		"b" :          -20.,
		"c" :          -55.,
		"d" :          91.,
		"V_peak" :      40.,
		"V_r" :         -80.,
		"V_t" :         -29.7,
		"k" :          1.0,
		"C_m" :          15.2,

		"Kappa" :      0.0289,
		"Lambda" :     0.331,
		#"alpha" :      0.032,

		"beta1" :      0.5,
		#"beta2" :      0.3,

		#'I_e': 230., # tuned

		# Values from Kumaravelu:
		"V_gaba" : -80.,
		"tau_gaba": 13.,

		# ???
		'nmda_ratio': 0.5,
		'tau_ampa': 3.0,
		'tau_nmda': 30.

		# Values from Fountas & Shanahan:
#		"V_gaba" : -60.,
#		"tau_gaba": 4.,
#		'nmda_ratio': 0.33,
#		'tau_ampa': 6.0,
#		'tau_nmda': 160.
	},
	'n': 10
}

cell_params['MSN_D2'] = {'model': 'izhikevich_dopa_modulated', 'params':
	{# Source: Fountas & Shanahan 2017
		"a" :          0.01,
		"b" :          -20.,
		"c" :          -55.,
		"d" :          91.,
		"V_peak" :      40.,
		"V_r" :         -80.,
		"V_t" :         -29.7,
		"k" :          1.0,
		"C_m" :          15.2,

		#"Kappa" :      0.0289,
		#"Lambda" :     0.331,
		"alpha" :      0.032,

		#"beta1" :      0.5,
		"beta2" :      0.3,

		#'I_e': 230., # tuned

		# Values from Kumaravelu:
		"V_gaba" : -80.,
		"tau_gaba": 13.,

		# ???
		'nmda_ratio': 0.5,
		'tau_ampa': 3.0,
		'tau_nmda': 30.

		# Values from Fountas & Shanahan:
#		"V_gaba" : -60.,
#		"tau_gaba": 4.,
#		'nmda_ratio': 0.33,
#		'tau_ampa': 6.0,
#		'tau_nmda': 160.
	},
	'n': 10
}

cell_params['STN'] = {'model': 'izhikevich_psc_alpha', 'params':
	{# Source: Mandali et al. 2015 + translation to Izhi2007b
		'C_m': 1.,
		'k': 0.04,
		'V_t': -42.3,
		'V_r': -82.7,
		'tau_syn_in': tau_in,
		'tau_syn_ex': tau_ex,
		'a': .005,
		'b': .265,
		'c': -65.,
		'd': 1.5,
		'I_e': 16.5 + 6. # 20.
	},
	'n': 10
}

cell_params['GPe'] = {'model': 'izhikevich_psc_alpha', 'params':
	{# Source: Mandali et al. 2015 + translation to Izhi2007b
		'C_m': 1.,
		'k': 0.04,
		'V_t': -42.3,
		'V_r': -82.7,
		'tau_syn_in': tau_in,
		'tau_syn_ex': tau_ex,
		'a': .1,
		'b': .2,
		'c': -65.,
		'd': 2.,
		'I_e': 16.5 + 7. # 10.
	},
	'n': 10
}

cell_params['GPi'] = {'model': 'izhikevich_psc_alpha', 'params':
	{# Source: Mandali et al. 2015 + translation to Izhi2007b
		'C_m': 1.,
		'k': 0.04,
		'V_t': -42.3,
		'V_r': -82.7,
		'tau_syn_in': tau_in,
		'tau_syn_ex': tau_ex,
		'a': .1,
		'b': .2,
		'c': -65.,
		'd': 2.,
		'I_e': 16.5 + 7. # 10.
	},
	'n': 10
}

cell_params['Pyr'] = {'model': 'izhikevich_psc_alpha', 'params':
	{# Source: Kumaravelu et al. 2016 + translation to Izhi2007b
		'C_m': 1.,
		'k': 0.04,
		'V_t': -42.3,
		'V_r': -82.7,
		'tau_syn_in': tau_in,
		'tau_syn_ex': tau_ex,
		'a': .02,
		'b': .2,
		'c': -65.,
		'd': 8.,
		'I_e': 16.5
	},
	'n': 10
}

cell_params['Inh'] = {'model': 'izhikevich_psc_alpha', 'params':
	{# Source: Kumaravelu et al. 2016 + translation to Izhi2007b
		'C_m': 1.,
		'k': 0.04,
		'V_t': -42.3,
		'V_r': -82.7,
		'tau_syn_in': tau_in,
		'tau_syn_ex': tau_ex,
		'a': .1,
		'b': .2,
		'c': -65.,
		'd': 2.,
		'I_e': 16.5 + 2. # Tuned
	},
	'n': 10
}

cell_params['Thal'] = {'model': 'izhikevich_psc_alpha', 'params':
	{# Source: Izhikevich & Edelman 2008 + synaptic adjustment
		'tau_syn_in': tau_in,
		'tau_syn_ex': tau_ex,
		'C_m': 200.,
		'k': 1.6,
		'V_r': -60.,
		'V_t': -50.,
		'V_peak': 40.,
		'a': .1,
		'b': 15.,
		'c': -60.,
		'd': 10.,
		't_ref': 1., # Tuned
		'I_e': 160. # 0.
	},
	'n': 10
}
