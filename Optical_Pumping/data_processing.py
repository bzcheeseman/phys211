%pylab inline

from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

data = genfromtxt("data/s4_4_scope_trace.tsv", skip_header=1)
data = delete(data, data[:,2], 1)

figure(figsize = (8, 8))
plot(data[:,1], data[:,0])
plot(data[:,1], .01*data[:,2])
show()


#Iv = [.514, .5, .45, .55, .6, .65, .7, .75, .8, .4, .35, .3, .25, .2, .15, .1, .05, 0, .325]

Iv = [0, .05, .1, .15, .2, .25, .3, .325, .35, .4, .45, .5, .514, .55, .6, .65, .7, .75, .8]

#f = [356, 349, 327, 372, 402, 438, 475, 512, 554, 313, 307, 307, 315, 330, 351, 380, 409, 444, 305]

f = [444, 409, 380, 351, 330, 315, 307, 305, 307, 313, 327, 349, 356, 372, 402, 438, 475, 512, 554]

plot(sort(Iv), f)
savefig("Earth_B_Field.png")
show()
