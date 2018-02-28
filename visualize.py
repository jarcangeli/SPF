import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.animation as manimation

import matplotlib.animation as animation
import numpy as np
import pylab as p
from matplotlib.backends.backend_pdf import PdfPages
import my_fns as f
import pars

def mass_to_size(marr):

    sizes = np.log(marr)
    sizes = sizes / np.mean(sizes)
    return sizes*10

def orbit_plot_points(xarr, marr, t, colours):

    p.figure()
    p.xlabel('Position in x')
    p.ylabel('Position in y')
    p.title('System at time {:.2e} years'.format(t/60/60/24/365))

    sizes = mass_to_size(marr)

    n_obj = pars.Np
    xarr = np.array(xarr)
    objects = [xarr[:,i] for i in range(n_obj)]
    for size, xy, c in zip(sizes, objects, colours):
        # plot planet position
        p.plot(xy[0], xy[1], marker='o', ms=size, color=c)

def orbit_plot_lines(xarrs, colours, **kwargs):
    # xarrs: time x object x xyz
    n_obj = pars.Np
 
    xarrs = np.array(xarrs)
    objects = [xarrs[:,:,i] for i in range(n_obj)]
    for obj, c in zip(objects, colours):
        xs = obj[:,0]
        ys = obj[:,1]
        p.plot(xs, ys, color=c, **kwargs)

    #p.show()

############
# PDF File #
############

def orbit_pdf(xarrs, marr, t0, dt, pdf_file='orbits.pdf'):

    print '~~~~~~~~~~~~~~ Saving orbits to pdf ~~~~~~~~~~~~~~~'
    
    f.silentremove(pdf_file)
    pdf = PdfPages(pdf_file)

    times = np.arange(len(xarrs))*dt + t0

    # set up colours for two objects
    colours = ['y', 'c']

    # Number of images to save
    n_im = 20
    skip = int(len(times)/n_im)
    indexes = range(len(times))[::skip]

    for i in indexes: # time step

        # plot the current position of the objects
        xarr, t = xarrs[i], times[i]
        orbit_plot_points(xarr, marr, t, colours)

        # plot the trajectories from time 0
        orbit_plot_lines(xarrs[:i+1], colours)

        pdf.savefig()
        p.close()

    pdf.close()
        











#################
# Movie attempt #
#################


def make_mp4(xarrs, marr, t0, dt, mp4_file='orbits.mp4'):

    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='Nbody', artist='jarcangeli',
                    comment='N-body orbit visualizer')
    writer = FFMpegWriter(fps=15, metadata=metadata)

    print '~~~~~~~~~~~~~~ Saving orbits to mp4 ~~~~~~~~~~~~~~~'
    f.silentremove(mp4_file)

    times = np.arange(len(xarrs))*dt + t0

    # set up colours for two objects
    colours = ['C1', 'C0']

    # Number of images to save
    n_im = 300
    skip = int(len(times)/n_im)
    indexes = range(len(times))[::skip]

    fig, ax = p.subplots()
    p.xlabel('Position in x')
    p.ylabel('Position in y')
    sizes = mass_to_size(marr)
    
    x_range = np.max(xarrs)*1.1
    p.xlim([-x_range,x_range])
    p.ylim([-x_range,x_range])
    ax.set_autoscale_on(False)

    bodies = []; trails = []; long_trails = [];
    for k in range(pars.Np):
        l, = p.plot([], [], ls='None', marker='o', color=colours[k], ms=sizes[k])
        bodies.append(l)

        tl, = p.plot([], [], ls='-', marker='None', color=colours[k], lw=3)
        trails.append(tl)

        ltl, = p.plot([], [], ls='-', marker='None', color=colours[k], lw=2, alpha=0.2)
        long_trails.append(ltl)

    with writer.saving(fig, mp4_file, n_im):
        for n_i in range(len(indexes)):
            i = indexes[n_i]

            # plot the current position of the objects
            xarr, t = xarrs[i], times[i]

            n_obj = pars.Np
            xarr = np.array(xarr)
            objects = [xarr[:,j] for j in range(n_obj)]
            for body, trail, xy in zip(bodies, trails, objects):
                # plot planet position
                #p.plot(xy[0], xy[1], marker='o', ms=size, color=c, ls='None')
                body.set_data(xy[0], xy[1])

            # Plot the trail
            trail_n = 5
            if n_i > trail_n:
                for j in range(n_obj):
                    # xarr is now a "steps x positions x bodies" array    
                    xarr = np.array(xarrs[indexes[n_i-trail_n]:i])
                    trail = trails[j]
                    trail_xys = xarr[:,:,j] # step x pstn for 1 body
                    trail.set_data(trail_xys[:,0], trail_xys[:,1])

            # Plot the long trail
            trail_n = 30
            if n_i > trail_n:
                for j in range(n_obj):
                    # xarr is now a "steps x positions x bodies" array    
                    xarr = np.array(xarrs[indexes[n_i-trail_n]:i])
                    trail = long_trails[j]
                    trail_xys = xarr[:,:,j] # step x pstn for 1 body
                    trail.set_data(trail_xys[:,0], trail_xys[:,1])

            p.title('System at time {:.2f} years'.format(t/pars.yr))

            # plot the trajectories from time 0
            #orbit_plot_lines(xarrs[:i+1], colours)
            writer.grab_frame()
    


           
