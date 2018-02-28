import two_body

rdict = {}
steps = [1e-2, 1e-3, 1e-4, 1e-5]
schemes = [ 'Forward Euler', 'Midpoint', 'Leapfrog', 'Leapfrog2', 'Hermite', 'RK4' ]

for scheme in schemes:
    errors = []
    for step in steps:
        error = two_body.run(scheme=scheme, step=step)
        errors.append(error)

    rdict[scheme] = errors

# Print output
col_width = max([len(scheme) for scheme in schemes])+2 #padding
print 'dt\t', '\t'.join([scheme.ljust(col_width) for scheme in schemes])
for i, step in enumerate(steps):
    errors = ['{:.1e}'.format(rdict[scheme][i]).ljust(col_width) for scheme in schemes]
    print '{:.0e}'.format(step), '\t', '\t'.join(errors)

'''
dt	    Forward Euler	Midpoint	Leapfrog	Leapfrog2	Hermite     RK4
1e-02 	-5.3e-01	    -6.1e-04	-1.0e-06	-1.0e-06	8.2e-04     5.4e-07
1e-03 	-1.7e-01	    -6.2e-07	-7.0e-11	-7.0e-11	8.2e-07     5.4e-12
1e-04 	-2.4e-02	    -6.2e-10	3.4e-13	    3.1e-13	    8.2e-10	    -0.0e+00
1e-05 	-2.5e-03	    -5.3e-13	1.7e-13	    2.6e-14	    8.2e-13	    -1.3e-13
'''
