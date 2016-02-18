%pylab inline

ydata = np.genfromtxt("data/cross_section/countrates.tsv", usecols=(0))

xdata_1 = np.genfromtxt("data/cross_section/countrates.tsv", dtype="string", usecols=(1))

xdata_cu = [0]
ydata_cu = [23.2394366197]
xdata_al = [0]
ydata_al = [23.2394366197]
for i in range(0, len(xdata_1)):
    try:
        if xdata_1[i].split("_")[0] == "cu":
            xdata_cu.append(float(xdata_1[i].split("_")[1])+.01*float(xdata_1[i].split("_")[2]))
            ydata_cu.append(ydata[i])
        elif xdata_1[i].split("_")[0] == "al":
            xdata_al.append(float(xdata_1[i].split("_")[1])+.01*float(xdata_1[i].split("_")[2]))
            ydata_al.append(ydata[i])
    except Exception:
        pass

plt.plot(xdata_al, ydata_al, 'o')
plt.xlabel("Absorber Thickenss (cm)")
plt.ylabel("Countrate (cout/sec)")
plt.title("Countrate Plots")
#plt.savefig("plots/countrate_cu.pdf")
