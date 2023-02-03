import numpy as np
import multiprocessing as mp

def model(tau, tmax=12.5, dt=0.001, withoutg1GR=False):
    '''
        Diffusion model for imitating the instability time of Mercury.
        Time is in units of Gyr.
        Use tau=1e-9 to imitate No GR and tau=1e1000 to imitate constant GR.
    '''
    # Define values for Mercury and Jupiter perihelion precession frequencies
    g1 = 5.577
    g5 = 4.257
    g1max = 5.72367
    g_GR = 0.4298
    if withoutg1GR: g1 -= g_GR
    
    # Define the diffusion constants based on the timescales for instability given by Mogavero & Laskar (2021) and Laskar & Gastineau (2009).
    tscale_GR = 27.6 # Gyr
    D_GR = (( (g1max - g5)**2.0 ) / (4.0 * tscale_GR))
    kick_GR = np.sqrt(2 * D_GR * dt)

    tscale_noGR = 3.216507262547807 # Gyr
    D_noGR = (( (g1max - g5)**2.0 ) / (4.0 * tscale_noGR))
    kick_noGR = np.sqrt(2 * D_noGR * dt)
    
    # Solve the advection-diffusion equation numerically.
    t = 0
    while g1 > g5 and t < tmax:
        # While GR decreases, remove additional precession from GR (advection).
        if (t < tau and (not withoutg1GR)): g1 -= g_GR * min(1, dt/tau)
        
        # Linearly interpolate between GR and no GR regimes.
        alpha = max(1. - t/tau,0)
        kick = kick_GR*alpha + (1.-alpha)*kick_noGR
        
        # Kick g1 frequency in random direction (diffusion).
        g1 += kick*np.random.normal()

        # Upper reflecting barrier at g1max.
        if g1 > g1max: g1 = g1max - (g1-g1max)
        t += dt


    return t

def run_model(taus, nsamples=1, tmax=12.5, dt=0.001, withoutg1GR=False, processes=None):
    '''
        Run the advection-diffusion model in parallel over the given taus.
        nsamples determines many times to run the model for each tau.
        Time is in units of Gyr.
    '''
    model_data = []
    results = lambda res: model_data.append(res)
    
    p = mp.Pool(processes=processes)
    for i in range(nsamples):
        for tau in taus:
            p.apply_async(model, args=(tau, tmax, dt, withoutg1GR,), callback=results)
    p.close()
    p.join()
    
    return np.array(model_data)
