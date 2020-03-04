'''Get observables analytically'''

import sys, time, warnings
from math import exp, log
import numpy as np
from scipy.integrate import ode
import scipy.special as sp
import scipy.optimize as opt
import matplotlib.pyplot as plt
import Run_Sobol as rs

def compute_I(R_p, R_m, c_1, c_2, c_3):
    '''Compute I integrals. Reworked a bit to match Andy's MATLAB.
    Note, the Python def of Ein is different than in MATLAB.
    '''

    return c_1*(log(R_p)-log(R_m)) + (c_1-c_2)*(sp.expi(-c_3*R_m)-sp.expi(-c_3*R_p))


def rootfind(rat, Rp, N, v, lam, greeks):
    Ival1 = compute_I(Rp, Rp/rat, greeks['theta'], greeks['beta'], greeks['delta'])
    Ival2 = compute_I(Rp, Rp/rat, greeks['eta'], greeks['alpha'], greeks['gamma'])
    return (N*lam/v)*Ival1 - log(rat)*Ival2


def compute_c(v, greeks, R_m, R_p):
    '''Get travel wave speed, reworked a bit to match Andy's MATLAB'''

    I_1 = compute_I(R_p, R_m, greeks['eta'], greeks['alpha'], greeks['gamma'])
    I_2 = compute_I(R_p, R_m, greeks['theta'], greeks['beta'], greeks['delta'])
    rat = I_1/I_2

    return v*rat/(1+rat)


def Lasymp(Rp, N, v, lam, greeks):
    a = greeks['alpha']
    b = (greeks['eta']-greeks['alpha'])*Fasymp(Rp,greeks['gamma']) - (N*lam/v)*greeks['beta']
    c = -(N*lam/v)*(greeks['theta']-greeks['beta'])*Fasymp(Rp,greeks['delta'])
    return (-b + np.sqrt(b**2-4*a*c))/(2*a)


def Fasymp(Rp, gamma):
    return np.euler_gamma + log(Rp*gamma) - sp.expi(-Rp*gamma)


def casymp(logRpRm, Rp, v, greeks):
    I1overI2 = (
               greeks['alpha']*logRpRm +
               (greeks['eta']-greeks['alpha'])*Fasymp(Rp,greeks['gamma'])
               )/(
               greeks['beta']*logRpRm +
               (greeks['theta']-greeks['beta'])*Fasymp(Rp,greeks['delta'])
               )
    return v*I1overI2/(1+I1overI2)


def find_c_and_R_minus(v, lam, greeks, N, R_p):
    ratmax = 10**20 # switchover point for Rp/Rm from root finding to asymptotic

    fun = lambda x_rat: rootfind(x_rat, R_p, N, v, lam, greeks)

    if fun(ratmax) < 0:
        ratout = opt.bisect(fun, 1+10**-8, ratmax, maxiter=300)
        Rm = R_p/ratout
        c = compute_c(v, greeks, Rm, R_p)
    else:
        logRpRm = Lasymp(R_p, N, v, lam, greeks)
        c = casymp(logRpRm, R_p, v, greeks)
        try:
            Rm = R_p/exp(logRpRm)
        except OverflowError:
            Rm = 0

    return (c, Rm)

###                                 ###
#####       START OLD STUFF       #####
###                                 ###

def n_func(v, lam, greeks, R_m, R_p):
    '''Ask Andy...'''

    I_1 = compute_I(greeks['eta'], greeks['alpha'], greeks['gamma'], R_m, R_p)
    I_2 = compute_I(greeks['theta'], greeks['beta'], greeks['delta'], R_m, R_p)

    return v/lam*I_1/I_2*log(R_p/R_m)


def find_R_minus_old(v, lam, greeks, N, R_p):
    '''Find it!'''

    f = lambda x, v, lam, greeks, R_p : n_func(v, lam, greeks, x, R_p) - N

    R_m = opt.bisect(f, 1e-300, R_p*(1-1e-15), args=(v, lam, greeks, R_p))

    return R_m


def solve_odes(x0, tstart, v, lam, greeks, c):
    '''Solve the ODEs with given initial condition, start time,
    end time, and parameters.
    We will actually solve twice: once to get the length of the relevant
    spatial domain, and then again to get the solution on a reasonable mesh.'''

    x_sol = [x0] # create a list to hold the solution. it starts with the IC
    tlist = [tstart] # create a list to hold the times which correspond to the solution points

    ###### First, solve to get approximate width of the domain ######
    ### setup solver
    solver = ode(ode_eqns).set_integrator('dopri5', nsteps=1000)
    # Pass in the IC, initial time, and parameter dictionary to the solver object
    solver.set_initial_value(x0, tstart).set_f_params(v, lam, greeks, c)
    ### solve with very coarse step size until end condition is met

    stime = time.perf_counter()
    while solver.successful and solver.y[1] >= 1e-8:
        solver.integrate(solver.t+100)
        if time.perf_counter() - stime > 5:
            raise RuntimeError('Coarse solver timed out.')
    # if this was way too coarse, resolve a bit within 100
    stime = time.perf_counter()
    if solver.t < 150:
        solver = ode(ode_eqns).set_integrator('dopri5')
        solver.set_initial_value(x0, tstart).set_f_params(v, lam, greeks, c)
        while solver.successful and solver.y[1] >= 1e-8:
            solver.integrate(solver.t+1)
            if time.perf_counter() - stime > 5:
                raise RuntimeError('Middle solver timed out.')
    ### use final t value to define coarsity of our solve routine
    dx = solver.t/2000

    ###### Solve the second time, on the mesh we will use ######
    ### setup solver again
    solver = ode(ode_eqns).set_integrator('dopri5')
    # Pass in the IC, initial time, and parameter dictionary to the solver object
    solver.set_initial_value(x0, tstart).set_f_params(v, lam, greeks, c)

    stime = time.perf_counter()
    while solver.successful and solver.y[1] >= 1e-8:
        solver.integrate(solver.t+dx) #integrate up to next time we want to record
        # record solution at that time
        x_sol.append(solver.y)
        tlist.append(solver.t)
        if time.perf_counter() - stime > 5:
            raise RuntimeError('Main solver timed out.')

    # After we've finished, return all the values for plotting, recording, etc.
    return (x_sol, tlist)

###                               ###
#####       END OLD STUFF       #####
###                               ###

def ode_eqns(t, x, v, lam, greeks, c):
    '''Define the ODEs for R and rho'''
    R_t = lam/v*((v-c)/c)*x[0]*x[1]
    rho_t = (
        k(greeks['eta'],greeks['alpha'],greeks['gamma'],x[0])/c -
        k(greeks['theta'],greeks['beta'],greeks['delta'],x[0])/(v-c)
        )*x[1]
    return np.array([-R_t, -rho_t])
    
def k(A,B,C,R):
    return A - (A-B)*exp(-C*R)


def solve_odes_hardstop(x0, tstart, v, lam, greeks, c):
    '''Solve the ODEs with given initial condition, start time,
    end time, and parameters.
    We will actually solve twice: once to get the length of the relevant
    spatial domain, and then again to get the solution on a reasonable mesh.'''

    x_sol = [x0] # create a list to hold the solution. it starts with the IC
    tlist = [tstart] # create a list to hold the times which correspond to the solution points

    ###### First, solve to get approximate width of the domain ######
    ### setup solver
    solver = ode(ode_eqns).set_integrator('dopri5', nsteps=1000)
    # Pass in the IC, initial time, and parameter dictionary to the solver object
    solver.set_initial_value(x0, tstart).set_f_params(v, lam, greeks, c)
    ### solve with very coarse step size until end condition is met or
    ###   we have gone 10 km
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        while solver.successful and solver.y[1] >= 1e-8 and solver.t < 10000:
            solver.integrate(solver.t+10)
            # record solution at that time
            x_sol.append(solver.y)
            tlist.append(solver.t)
    
    ### use final t value to define coarsity of our solve routine
    dx = solver.t/1000

    ###### Solve the second time, on the mesh we will actually use... ######
    ######   ...unless this is a waste of time.                       ######
    if solver.t < 10000:
        ### setup solver again
        x_sol = [x0]
        tlist = [tstart]
        solver = ode(ode_eqns).set_integrator('dopri5')
        # Pass in the IC, initial time, and parameter dictionary to the solver object
        solver.set_initial_value(x0, tstart).set_f_params(v, lam, greeks, c)

        while solver.successful and solver.y[1] >= 1e-8 and solver.t < 10000:
            solver.integrate(solver.t+dx) #integrate up to next time we want to record
            # record solution at that time
            x_sol.append(solver.y)
            tlist.append(solver.t)
    
    # After we've finished, return all the values for plotting, recording, etc.
    return (x_sol, tlist)

def ratio(v,c):
    '''Stationary over Moving'''
    return v/c - 1


def get_observables(speed, R_minus, ratio, wave_form, mesh):
    ### Peak size
    peak_size = np.max(wave_form)

    ### Width, as std
    mass = np.trapz(wave_form, mesh)
    mean = np.trapz(wave_form*mesh, mesh)/mass

    # integrate to get variance, take sqrt
    width = np.sqrt(np.trapz(wave_form*(mesh-mean)**2, mesh)/mass)

    # Traveling wave speed
    # Peak size
    # Width (btwn 12 locusts/m^2)
    # Ratio of total stationary locusts to total moving
    # Resource density behind locusts, asymptotic
    return (speed, peak_size, width, ratio, R_minus)


def main(lam, v, N, Rplus, greeks):
    '''Grab info about the asymptotic solution and return'''

    #import pdb; pdb.set_trace()
    tstart = 0

    # R_m = find_R_minus(v, lam, greeks, N, Rplus)
    # c = compute_c(v, lam, greeks, R_m, Rplus)
    c, R_m = find_c_and_R_minus(v, lam, greeks, N, Rplus)
    x0 = np.array([Rplus, 1e-6])

    try:
        x_sol, mesh = solve_odes_hardstop(x0, tstart, v, lam, greeks, c)
    except:
        # If the ODEs fail for any reason, still return c and R_m
        return (sys.exc_info()[1], 'c = {}'.format(c), 'R_m = {}'.format(R_m), 
                None, None)

    # return wave speed, skew, R_minus, traveling wave profile, spatial mesh
    #import pdb; pdb.set_trace()
    return (c, R_m, np.array(x_sol, order='F')[:,1], np.array(mesh))


def test_me():

    ### Biological Parameters ###
    # lam = 10**(-5.084778)
    # v = 0.041886
    # N = 9502.868652
    # R_p = 209.351196
    # beta = 0.020937
    # theta_beta = 9.762344
    # delta = 0.029516
    # Delta = 0.111419
    # eta_alpha = 0.309109
    # gamma = 0.037436
    # greeks = {}
    # greeks['beta'] = beta
    # greeks['theta'] = beta*theta_beta
    # greeks['delta'] = delta

    # greeks['alpha'] = Delta*beta*greeks['theta']/(greeks['theta']-beta*eta_alpha)
    # greeks['eta'] = greeks['alpha']*eta_alpha
    # greeks['gamma'] = gamma

    ### Test case 1 ###
    # greeks = {}
    # greeks['alpha'] = 2.0
    # greeks['eta'] = 1.0
    # greeks['gamma'] = .01
    # greeks['beta'] = 1
    # greeks['theta'] = 2
    # greeks['delta'] = .01
    # v = 0.45
    # lam = 1e-5
    # N = 10000
    # R_p = 200

    ### Test case 2 ###
    greeks = {}
    greeks['alpha'] = 2.0
    greeks['eta'] = 1.0
    greeks['gamma'] = .04
    greeks['beta'] = 4.0
    greeks['theta'] = 5
    greeks['delta'] = .01
    v = 0.45
    lam = 1e-5
    N = 10000
    R_p = 200

    tstart = 0

    # R_m = find_R_minus(v, lam, greeks, N, R_p)
    # c = compute_c(v, lam, greeks, R_m, R_p)
    c, R_m = find_c_and_R_minus(v, lam, greeks, N, R_p)
    x0 = np.array([R_p, 1e-6])

    x_sol, mesh = solve_odes_hardstop(x0, tstart, v, lam, greeks, c)

    mesh = np.array(mesh)
    x_sol = np.array(x_sol)

    ratio_res = ratio(v,c)

    speed, peak_size, width, ratio_res, R_minus = get_observables(c, R_m, ratio_res,
                                np.array(x_sol, order='F')[:,1], np.array(mesh))

    print('speed = {}'.format(speed))
    print('peak_size = {}'.format(peak_size))
    print('width = {}'.format(width))
    print('ratio = {}'.format(ratio_res))
    print('R_minus = {}'.format(R_minus))

    plt.subplot(2,1,1)
    plt.plot(mesh, x_sol[:,0], label='R')
    r_m_mesh = np.ones_like(mesh)*R_m
    plt.plot(mesh, r_m_mesh, label='R_m')
    r_p_mesh = np.ones_like(mesh)*R_p
    plt.plot(mesh, r_p_mesh, label='R_p')
    plt.legend()
    plt.title('R')

    plt.subplot(2,1,2)
    plt.plot(mesh, x_sol[:,1])
    plt.title(r'$\rho$')

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    test_me()