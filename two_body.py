# Main program file
import numpy as np
import functions as fn
import pars
import visualize
import time as time_module

reload(fn)

# List of supported schemes
schemes = [ 'Forward Euler', 'Midpoint', 'Leapfrog', 'Leapfrog2', 'Hermite', 'RK4' ]

def run(scheme='RK4', step=1e-2, pdf=False, mp4=False):

    t0 = time_module.time()
    print '\n~~~~~~~~~~~~~~ Running two_body program ~~~~~~~~~~~~~~~\n'
    print 'Using {} scheme'.format(scheme)
    # assign the initial positions, velocities and masses
    xarr, varr, marr = fn.setup_2body(ecc=0.)

    # decelerations
    time = 0 # start time
    tfinal = 1e8 #pars.yr*10
    dt = step*pars.yr
    print '\ttimestep of {:.0e} years'.format(dt/pars.yr)
    print '\tend time of {:.3f} years'.format(tfinal/pars.yr)

    # compute the total energy, used for verification
    etot0 = fn.e_tot(xarr, varr, marr)

    # start main loop
    # store positions and velocities
    xarrs, varrs = [], []
    while time < tfinal:
        # calculate forces, update positions and velocities
        # this involves calling the function that computes the accelerations

        xarr, varr = fn.update(xarr, varr, marr, dt, scheme=scheme)

        xarrs.append(xarr.copy())
        varrs.append(varr.copy())

        # increment time
        time += dt
    #end

    etot = fn.e_tot(xarr, varr, marr)
    error = (etot - etot0) / etot0


    print 'Initial total energy of {:.2e}, final total energy of {:.2e}.\nError of {:.2e}\n'.format(etot0, etot, error)
    print 'Time taken: {:.2f}s\n'.format(time_module.time()-t0)

    #visualize.orbit_plot_lines(xarrs, ['y', 'k'])

    if pdf:
        visualize.orbit_pdf(xarrs, marr, 0, dt)
    if mp4:
        fname = '{}_{:.0g}.mp4'.format(scheme.lower(), dt/pars.yr)
        visualize.make_mp4(xarrs, marr, 0, dt, mp4_file=fname)

    return error

if __name__ == '__main__':
    run()









