%pylab inline

ydata = np.genfromtxt("data/cross_section/countrates.tsv", usecols=(0))

xdata_1 = np.genfromtxt("data/cross_section/countrates.tsv", dtype="string", usecols=(1))

xdata = []
for i in range(0, len(xdata_1)):
    try:
        xdata.append(float(xdata_1[i].split("_")[1])+.01*float(xdata_1[i].split("_")[2]))
    except Exception:
        xdata.append(0)

print xdata

plt.plot(xdata, ydata, 'o')
plt.show()
