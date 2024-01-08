import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import math as m
from imports import var
from panel_integrals import ijintegral, klintegral
from streamlines import ijstreamline, klstreamline
from tqdm import tqdm

mpl.rcParams['lines.linewidth'] = 0.8
plt.rcParams["figure.autolayout"] = True
ratio = 1
width = 6
dpi = 750
form = "png"

#variables

aoa = var("angle_of_attack")
gridsize = var("grid_size")
param = str(var("naca"))

#import foil coordinates

xarr = []
yarr = []

file = open(f"foils/{param}.dat", "r", encoding="utf8")
data = file.readlines()
i = 1
for _ in range(len(data) - 1):
    linedat = data[i]
    linedat = linedat.split()
    xarr.append(float(linedat[0]))
    yarr.append(float(linedat[1]))
    i += 1
file.close()

npanels = len(xarr) - 1

#check panel direction

edge = np.zeros(npanels)
for i in range(npanels):
    edge[i] = (xarr[i + 1] - xarr[i]) * (yarr[i + 1] + yarr[i])

if np.sum(edge) < 0:
    xarr = np.flipud(xarr)
    yarr = np.flipud(yarr)

#compute panel geometry

xcpoint = np.zeros(npanels)
ycpoint = np.zeros(npanels)
lenpanel = np.zeros(npanels)
phipanel = np.zeros(npanels)

for i in range(npanels):
    xcpoint[i] = (xarr[i] + xarr[i + 1]) / 2
    ycpoint[i] = (yarr[i] + yarr[i + 1]) / 2
    dx = xarr[i + 1] - xarr[i]
    dy = yarr[i + 1] - yarr[i]
    lenpanel[i] = (dx**2 + dy**2)**0.5
    phipanel[i] = m.atan2(dy, dx)
    if phipanel[i] < 0:
        phipanel[i] = phipanel[i] + 2 * np.pi

#normal angles

delta = phipanel + np.pi / 2
beta = delta - (aoa * (np.pi / 180))

#plot

fig = plt.gcf()
naca = var("naca")
plt.fill(xarr, yarr, color="black", label=f"naca {naca}{npanels} panels")
for i in range(npanels):
    plt.plot([xcpoint[i], xcpoint[i] + 5 * lenpanel[i] * np.cos(delta[i])],
             [ycpoint[i], ycpoint[i] + 5 * lenpanel[i] * np.sin(delta[i])], color="red")
plt.ylabel("")
plt.legend(loc='upper left', framealpha=1)
plt.grid()
fig.set_size_inches(width, width * ratio)
plt.xlim([0, 1])
plt.ylim([-0.5, 0.5])
plt.savefig(f"plots/panels.{form}", dpi=dpi)
plt.show()