import numpy as np
import pars

def setup_2body(ecc=0.5):
    '''
    Construct the 2 body problem, initialize at apocenter.
    ecc: eccentricity
    '''
    
    # declare the parameters as zeros
    xarr = np.zeros((3, pars.Np))
    varr = np.zeros((3, pars.Np))
    marr = np.zeros(pars.Np)

    # Consider the Sun (particle 0) and Earth (particle 1)
    marr[0] = pars.mSun
    marr[1] = pars.mEarth

    # Keplerian velocity corresponding to 1 AU around Sun (particle 0)
    vKep = np.sqrt(pars.gN*marr[0] / pars.au)

    # Initialize at apocenter
    xarr[:,1] = [pars.au*(1+ecc), 0., 0.]
    varr[:,1] = [0., vKep*np.sqrt((1-ecc)/(1+ecc)), 0.]

    return xarr, varr, marr


def forces(xarr, marr):
    '''
    xarr(3,Np): positions
    marr(Np): masses
    
    Calculate the gravitational force (accelerations)
    on each particle

    returns the acclerations
    '''
    acc = np.zeros((3,pars.Np))
    for i in range(pars.Np):
        for j in range(i+1, pars.Np):
            rji = xarr[:,j] - xarr[:,i] # relative position (vector)
            r2 = (rji**2).sum(axis=0) # squared distance
            r1 = np.sqrt(r2)   # distance
            r3 = r1*r2 # cubed distance
            
            force = pars.gN*rji/r3
            acc[:,i] += force*marr[j] # add to i
            acc[:,j] -= force*marr[i] # reverse sign 
            
    return acc

def e_tot(xarr, varr, marr):
    '''
    xarr(3,Np): positions
    varr(3,Np): velocities
    marr(Np): masses
    
    Calculate the energy of the system

    returns the energy (scalar)
    '''
    kinetic = 0.5 * np.sum(np.multiply(marr, varr**2))
    
    potential = 0
    for i in range(pars.Np):
        for j in range(i+1, pars.Np):
            rji = xarr[:,j] - xarr[:,i] # relative position (vector)
            r2 = (rji**2).sum(axis=0) # squared distance
            r1 = np.sqrt(r2)   # distance

            potential -= pars.gN*marr[i]*marr[j]/r1 # r1 is already strictly positive

    return kinetic + potential # total energy

def forces_hermite(xarr, varr, marr):
    '''
    Compute acceleration and jerk terms for Hermite scheme.
    '''
    acc = np.zeros((3, pars.Np)) # accelerations
    jer = np.zeros((3, pars.Np)) # time-derivative of acc
    for i in range(pars.Np):
        for j in range(i+1, pars.Np):
            rji = xarr[:,j] - xarr[:,i] # relative position
            vji = varr[:,j] - varr[:,i] # relative velocity

            r2 = sum(rji**2)
            r1 = np.sqrt(r2)
            r3 = r1*r2

            rv = sum(rji*vji) # dot product
            rv /= r2 # divide by r**2
            
            force = pars.gN * rji / r3
            acc[:,i] += force*marr[j]
            acc[:,j] -= force*marr[i]
            
            # add the jerk terms
            jerk = pars.gN * (vji - 3*rv*rji) / r3
            jer[:,i] += jerk*marr[j]
            jer[:,j] -= jerk*marr[i]
    
    return acc, jer            

#############
## Schemes ##
#############

def update(xarr, varr, marr, dt, scheme):

    if scheme == 'Forward Euler':
        acc = forces(xarr, marr)
        xarr += varr*dt
        varr += acc*dt

    elif scheme == 'Midpoint':
        acc = forces(xarr, marr)
        xmid = xarr + varr*dt/2
        vmid = varr + acc*dt/2
        amid = forces(xmid, marr)
        xarr += vmid*dt
        varr += amid*dt

    elif scheme == 'Leapfrog':
        acc = forces(xarr, marr)
        varr += acc*dt/2
        xarr += varr*dt
        acc = forces(xarr, marr)
        varr += acc*dt/2

    elif scheme == 'RK4':

        vk1 = varr
        ak1 = forces(xarr, marr)

        vk2 = varr + ak1*dt/2.
        ak2 = forces(xarr+vk1*dt/2., marr)

        vk3 = varr + ak2*dt/2.
        ak3 = forces(xarr+vk2*dt/2., marr)

        vk4 = varr + ak3*dt
        ak4 = forces(xarr+vk3*dt, marr)

        xarr = xarr + (vk1+2*vk2+2*vk3+vk4)*dt/6.
        varr = varr + (ak1+2*ak2+2*ak3+ak4)*dt/6.


    elif scheme == 'Alternative Leapfrog' or scheme == 'Leapfrog2':
        acc = forces(xarr, marr)
        old_x, old_v, old_a = xarr.copy(), varr.copy(), acc.copy()
        
        xarr += varr*dt + acc*dt**2/2 # predicted position
        acc = forces(xarr, marr)
        varr += (acc+old_a)*dt/2
        xarr = old_x + (old_v+varr)*dt/2 + (old_a-acc)*dt**2/4 # corrected position

    elif scheme == 'Hermite':
        acc, jer = forces_hermite(xarr, varr, marr)
        old_x, old_v, old_a, old_j = xarr.copy(), varr.copy(), acc.copy(), jer.copy()

        # predict
        predict_x = old_x + old_v*dt + old_a*(dt**2)/2. + old_j*(dt**3)/6.
        predict_v = old_v + old_a*dt + old_j*(dt**2)/2.
        predict_a, predict_j = forces_hermite(predict_x, predict_v, marr)

        # correct
        varr = old_v + (old_a + predict_a)*dt/2. + (old_j-predict_j)*(dt**2)/12.
        xarr = old_x + (old_v + predict_v)*dt/2. + (old_a-predict_a)*(dt**2)/12.
        

    else:
        print '{} scheme not supported'
        raise # just stops the run

    return xarr, varr
        

        




