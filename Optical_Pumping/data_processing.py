%pylab inline

from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

data = genfromtxt("data/s4_3_scope_trace.tsv", skip_header=1)
data = delete(data, data[:,2], 1)

figure(figsize = (8, 8))
plot(data[:,1], data[:,0])
plot(data[:,1], .01*data[:,2])
show()
