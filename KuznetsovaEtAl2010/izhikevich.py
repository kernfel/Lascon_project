import nest
import numpy as np
import pylab

t_m20, V_m20 = np.loadtxt('step_-10pA.dat', unpack=True, skiprows=2)
t_p20, V_p20 = np.loadtxt('step_50pA.dat', unpack=True, skiprows=2)
t_100, V_100 = np.loadtxt('step_400pA.dat', unpack=True, skiprows=2)
t_rec = (t_m20, t_p20, t_100)
V_rec = (V_m20, V_p20, V_100)
nStim = 3

nest.SetKernelStatus({"resolution":0.01})
simulation_time = 10000.0

nest.ResetKernel()

nProbe = 1
n = nest.Create( 'izhikevich', nStim*nProbe, {'V_m': -62.6,
									'a': .0024, # higher: Squash top, squash bottom
									'b': .42, # Squash top, invariant bottom
									'c': -60.,
									'd': 50.,
									'V_th': 10.} )

ps = [0.]#np.linspace(0., 10., nProbe)
for (i,p) in enumerate(ps):
	nest.SetStatus([n[0] + i + j*nProbe for j in range(nStim)], {'I_e': p})

stim = [nest.Create('step_current_generator',
	params = {'amplitude_times': [2000., 6000.], 'amplitude_values': [i, 0.]}) for i in [-20., 20., 100.]]
for i in range(nStim):
	nest.Connect(stim[i], range(n[0] + i*nProbe, n[0] + (i+1)*nProbe))

multimeter = nest.Create( "multimeter" )
nest.SetStatus( multimeter, { "withtime" : True, "record_from" : ["V_m", 'U_m'] })
nest.Connect( multimeter, n )
meter_shape = (9999, len(n))

nest.Simulate( simulation_time )


multimeter_events = nest.GetStatus( multimeter )[0]['events']
voltages = np.reshape(multimeter_events["V_m"], meter_shape)
u = np.reshape(multimeter_events['U_m'], meter_shape)
times = np.reshape(multimeter_events["times"], meter_shape).T[0]

for row in range(nStim):
	for col in range(nProbe):
		i = row*nProbe + col
		pylab.subplot(nStim, nProbe, i + 1)
		pylab.plot(t_rec[row], V_rec[row])
		pylab.plot(times, u.T[i])
		pylab.plot(times, voltages.T[i])
#		pylab.ylim(-70,20)
		pylab.title(ps[col])

pylab.show()
