import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.rcParams['lines.linewidth'] = 0.8
plt.rcParams["figure.autolayout"] = True
ratio = 0.3
width = 6
dpi = 750
form = "png"

#naca profile calc

aoa = 0
resolution = 10
param = "4412"

resolution = int(round(resolution / 2, 0))
aoa = np.radians(- aoa)
naca = list(param)
m = int(naca[0]) / 100
p = int(naca[1]) / 10
e = int(naca[2] + naca[3]) / 100

xupper = []
yupper = []
xlower = []
ylower = []
meanx = []
meany = []
normx = []
normy = []

def xy(t):
    return 5 * e * (0.2969 * t**0.5
                    - 0.1260 * t
                    - 0.3516 * t**2
                    + 0.2843 * t**3
                    - 0.1015 * t **4)

def mean(t):
    if t <= p and t >= 0:
        return (m / p**2) * (2 * p * t - t**2)
    elif t > p and t <= 1:
        return (m / (1 - p)**2) * ((1 - 2 * p) + 2 * p * t - t**2)

def theta(t):
    if t <= p and t >= 0:
        return np.arctan(((2 * m) / p**2) * (p - t))
    elif t > p and t <= 1:
        return np.arctan(((2 * m) / (1 - p)**2) * (p - t))

#generate naca coordinates

if p == 0:
    for i in range(resolution):
        t = (i / resolution)**1.5
        xupper.append(t)
        xlower.append(t)
        yupper.append(xy(t))
        ylower.append(-xy(t))
        meanx.append(t)
        meany.append(0)
else:
    for i in range(resolution):
        t = (i / resolution)**(1.5 + 0.25 * (1 / p))
        xu = t - xy(t) * np.sin(theta(t))
        yu = mean(t) + xy(t) * np.cos(theta(t))
        xl = t + xy(t) * np.sin(theta(t))
        yl = mean(t) - xy(t) * np.cos(theta(t))
        xupper.append(xu)
        yupper.append(yu)
        xlower.append(xl)
        ylower.append(yl)
        meanx.append(t)
        meany.append(mean(t))
        u = (i / resolution)**1.5
        normx.append(u)
        normy.append(xy(u))

#visual rotation of the foil (not flow relative angle of attack!)

if aoa != 0:
    for i in range(len(xupper)):
        x = xupper[i]
        y = yupper[i]
        xupper[i] = (x - 0.5) * np.cos(aoa) - y * np.sin(aoa) + 0.5
        yupper[i] = y * np.cos(aoa) + (x - 0.5) * np.sin(aoa)
        x = xlower[i]
        y = ylower[i]
        xlower[i] = (x - 0.5) * np.cos(aoa) - y * np.sin(aoa) + 0.5
        ylower[i] = y * np.cos(aoa) + (x - 0.5) * np.sin(aoa)
        x = meanx[i]
        y = meany[i]
        meanx[i] = (x - 0.5) * np.cos(aoa) - y * np.sin(aoa) + 0.5
        meany[i] = y * np.cos(aoa) + (x - 0.5) * np.sin(aoa)
    for i in range(len(normx)):
        x = normx[i]
        y = normy[i]
        normx[i] = (x - 0.5) * np.cos(aoa) - y * np.sin(aoa) + 0.5
        normy[i] = y * np.cos(aoa) + (x - 0.5) * np.sin(aoa)

#plots

fig = plt.gcf()
plt.plot(xupper, yupper, "-", color="navy", label=f"naca {param}")
plt.plot(xlower, ylower, "-", color="navy")
if len(meanx) > 0:
    plt.plot(meanx, meany, "-", color="darkred", label="ligne de cambrure")
    plt.plot(normx, normy, "--", color="darkgreen",
             label=f"normale naca 00{int(e * 100)}")
plt.hlines(0, -0.25, 1.25, linestyles="--", color="black")
plt.xlabel("")
plt.ylabel("")
#plt.yticks(np.arange(-0.2, 0.2, 0.2))
plt.grid()
#plt.text(0.44, 0.15, "extrados")
#plt.text(0.44, -0.1, "intrados")
#plt.xlabel("corde")
#plt.legend(loc='upper left', framealpha=1)
fig.set_size_inches(width, width*ratio)
plt.xlim([0, 1])
plt.ylim([-0.1, 0.2])
plt.savefig(f"plot_{param}.{form}", dpi=dpi)
plt.show()
plt.close()

xlower=xlower[::-1]
ylower=ylower[::-1]

for i in range(len(xlower)):
    print("\\filldraw[black] " + "(" + str(round(xlower[i]*10, 2)) + "," + str(round(ylower[i]*10, 2)) + ") (1.5pt)")
for i in range(len(xupper)):
    print("\\filldraw[black] " + "(" + str(round(xupper[i]*10, 2)) + "," + str(round(yupper[i]*10, 2)) + ") (1.5pt)")
