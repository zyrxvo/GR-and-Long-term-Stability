import rebound
import reboundx # !!! Uses a modified version of REBOUNDx, https://github.com/zyrxvo/reboundx/tree/GR_Sweep !!!
import numpy as np
from sys import argv

twopi = 2.*np.pi
tmax = twopi * 12.5e9  # Maximum integration time.
dt = np.sqrt(11)*twopi/365.25 # Timestep of about 3.3 days.
save_interval = int(2e5*twopi/dt) # save snapshots every 200,000 years
n = int(argv[1]) # Simulation number.

delta = 0.00163677264386899684
taus = np.insert(np.array([10**(8 + n*delta) for n in range(1279)]), 0,0) # Assuming 1280 simulations, with taus that span from 1e8 yrs to 2.5e11 yrs.
tau = twopi * taus[n] # Rate of change for the strength of the General Relativistic correction term.

folder = 'archives/'
filename = folder + 'gr_{0:04d}.bin'.format(n)

try:
    sim = rebound.Simulation(filename)
    sim.automateSimulationArchive(filename, step=save_interval, deletefile=False)
except:
    try:
        sim = rebound.Simulation('ss.bin')
        sim.automateSimulationArchive(filename, step=save_interval, deletefile=True)
    except:
        sim = rebound.Simulation()
        sim.add(['sun', 'mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune'], date='2000-01-01 12:00')
        sim.save('ss.bin')
        sim.automateSimulationArchive(filename, step=save_interval, deletefile=True)

sim.move_to_com()
sim.integrator = 'whckl'
sim.dt = dt
sim.ri_whfast.safe_mode = 0
sim.ri_whfast.keep_unsynchronized = True
sim.exit_min_distance = 1e-2 # ~4x the distance between the Earth and Moon.
sim.exit_max_distance = 1000. # 1000 AU, Unlikely that a planet is bound.

rebx = reboundx.Extras(sim)
gr = rebx.load_force('gr_potential')
rebx.add_force(gr)
gr.params['c'] = tau # Set the rate of change for the strength of the General Relativistic correction term.
# !!! Uses a modified version of REBOUNDx, https://github.com/zyrxvo/reboundx/tree/GR_Sweep !!!

fate = ''
try: sim.integrate(tmax, exact_finish_time=False)
except rebound.Escape as esc: fate = 'esc'
except rebound.Encounter as enc: fate = 'enc'

endname = 'completed/gr_{0:04d}_{1}.bin'.format(n, fate)
sim.save(endname)