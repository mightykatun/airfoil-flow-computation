import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import path
import numpy as np
import math
import time as tm
import concurrent.futures as ccft
from tqdm import tqdm
from imports import var
from xfoil import xfoil
from panel_integrals import ijintegral, klintegral
from streamlines import ijstreamline, klstreamline

starta = tm.time()
print("Process started...")

mpl.rcParams['lines.linewidth'] = 0.8
plt.rcParams["figure.autolayout"] = True
ratio = 1
width = 6
dpi = 750
form = "png"

#variables

aoa = np.radians(var("angle_of_attack", "float"))
param = var("naca", "str")
vinf = var("v_infinity", "float")
resolution = var("panel_number", "int")
if resolution > 500:
    print("Panel number should be smaller than 500.")
    #exit(1)

#generate foil file with config.txt params

resolution = int(resolution // 2)
naca = list(param)
diff = 4 - len(naca)
for i in range(diff):
    naca[:0] = ["0"]
param = naca[0] + naca[1] + naca[2] + naca[3]
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

print(f"Generating naca {param}...")
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

#dump foil to .dat file

file = open(f"foils/naca_{param}.dat", "w", encoding="utf8")
i = 0
xupper.reverse()
yupper.reverse()
file.writelines(f"naca {param}" + "\n")
for _ in xupper:
    file.writelines(f"{xupper[i]:.10f} {yupper[i]:.10f}" + "\n")
    i += 1
i = 0
for _ in xlower:
    file.writelines(f"{xlower[i]:.10f} {ylower[i]:.10f}" + "\n")
    i += 1
file.close()

#visual rotation of the foil (not flow relative angle of attack!)

for i in range(len(xupper)):
    x = xupper[i]
    y = yupper[i]
    xupper[i] = (x - 0.5) * np.cos(- aoa) - y * np.sin(- aoa) + 0.5
    yupper[i] = y * np.cos(- aoa) + (x - 0.5) * np.sin(- aoa)
    x = xlower[i]
    y = ylower[i]
    xlower[i] = (x - 0.5) * np.cos(- aoa) - y * np.sin(- aoa) + 0.5
    ylower[i] = y * np.cos(- aoa) + (x - 0.5) * np.sin(- aoa)
    x = meanx[i]
    y = meany[i]
    meanx[i] = (x - 0.5) * np.cos(- aoa) - y * np.sin(- aoa) + 0.5
    meany[i] = y * np.cos(- aoa) + (x - 0.5) * np.sin(- aoa)
for i in range(len(normx)):
    x = normx[i]
    y = normy[i]
    normx[i] = (x - 0.5) * np.cos(- aoa) - y * np.sin(- aoa) + 0.5
    normy[i] = y * np.cos(- aoa) + (x - 0.5) * np.sin(- aoa)

#plots

fig = plt.figure(1)
plt.plot(xupper, yupper, "-", color="blue", label=f"naca {param} alpha {np.degrees(aoa):.1f}")
plt.plot(xlower, ylower, "-", color="blue")
if len(meanx) > 0:
    plt.plot(meanx, meany, "-", color="red", label="camber line")
    if m != 0 and p != 0:
        plt.plot(normx, normy, "--", color="green",
                 label=f"normal naca 00{int(e * 100)}")
plt.hlines(0, -0.25, 1.25, linestyles="--", color="black")
plt.legend(loc='upper left', framealpha=1)
plt.grid()
fig.set_size_inches(width, width*ratio)
plt.xlim([-0.1, 1.1])
plt.ylim([-0.6, 0.6])
plt.savefig(f"plots/plot_{param}.{form}", dpi=dpi)
if var("foilgen_plot", "str") == "True":
    plt.show()
plt.close()

#xfoil predictions

if var("xfoil_run", "str") == "True":
    print("Running XFoil...")
    try:
        xfoil(param)
    except:
        print("XFoil computation failed, try lowering alpha and/or panel number.")
        exit(1)
    else:
        x_cp_xfoil, y_cp_xfoil, cp_xfoil, lift_coeff_xfoil, drag_coeff_xfoil, moment_coeff_xfoil = xfoil(param)
else:
    print("Skipping XFoil...")
    x_cp_xfoil = []
    y_cp_xfoil = []
    cp_xfoil = []
    lift_coeff_xfoil = 0
    drag_coeff_xfoil = 0
    moment_coeff_xfoil = 0

#panel geometries

print("Panelling...")

xarr = []
yarr = []

file = open(f"foils/naca_{param}.dat", "r", encoding="utf8")
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
    phipanel[i] = math.atan2(dy, dx)
    if phipanel[i] < 0:
        phipanel[i] = phipanel[i] + 2 * np.pi

#normal angles

delta = phipanel + np.pi / 2
beta = delta - aoa
beta[beta > 2 * np.pi] = beta[beta > 2 * np.pi] - 2 * np.pi

#plot

fig = plt.figure(2)
plt.fill(xarr, yarr, color="black", label=f"naca {param}\n{npanels} panels")
for i in range(npanels):
    plt.plot([xcpoint[i], xcpoint[i] + 5 * lenpanel[i] * np.cos(delta[i])],
             [ycpoint[i], ycpoint[i] + 5 * lenpanel[i] * np.sin(delta[i])], color="red")
plt.legend(loc='upper left', framealpha=1)
plt.grid()
fig.set_size_inches(width, width * ratio)
plt.xlim([-0.1, 1.1])
plt.ylim([-0.6, 0.6])
plt.savefig(f"plots/panels_{param}.{form}", dpi=dpi)
if var("panel_plot", "str") == "True":
    plt.show()
plt.close()

#integral computation

print("Computing integrals...")
with ccft.ThreadPoolExecutor() as executor:
    ij = executor.submit(ijintegral, xarr, yarr, xcpoint, ycpoint, phipanel, lenpanel, npanels)
    kl = executor.submit(klintegral, xarr, yarr, xcpoint, ycpoint, phipanel, lenpanel, npanels)
    I, J = ij.result()
    K, L = kl.result()
executor.shutdown()

A = I + np.pi * np.eye(npanels, npanels)
A = np.hstack((A, np.zeros((npanels, 1)))) #make space for gamma
A = np.vstack((A, np.zeros((1, npanels + 1)))) #make space for kutta condition
for i in range(npanels): #populate with kutta condition
    A[npanels, i] = J[0, i] + J[npanels - 1, i]
A[npanels, npanels] = - (sum(L[0, :] + L[npanels - 1, :])) + 2 * np.pi

B = - vinf * 2 * np.pi * np.cos(beta)
B = np.append(B, - vinf * 2 * np.pi * (np.sin(beta[0]) + np.sin(beta[npanels - 1]))) #kutta condition

AB = np.linalg.solve(A, B) #solve panel source and vortex strength

l = AB[0:len(AB) - 1]
g = AB[len(AB) - 1]

#panel velocities

vtan = np.zeros(npanels)
cp = np.zeros(npanels)

for i in range(npanels):
    temp1 = vinf * np.sin(beta[i])
    temp2 = (1 / (2 * np.pi)) * sum(l * J[i, :])
    temp3 = g / 2
    temp4 = - (g / (2 * np.pi)) * sum(L[i, :])
    vtan[i] = temp1 + temp2 + temp3 + temp4
    cp[i] = 1 - (vtan[i] / vinf)**2

normal_force = - cp * lenpanel * np.sin(beta)
axial_force = - cp * lenpanel * np.cos(beta)

lift_coeff = sum(normal_force * np.cos(aoa)) - sum(axial_force * np.sin(aoa))
drag_coeff = sum(normal_force * np.sin(aoa)) + sum(axial_force * np.cos(aoa))
moment_coeff = sum(cp * (xcpoint - 0.25) * lenpanel * np.cos(phipanel))

temp_path = f"foils/results_{param}_{np.degrees(aoa):.1f}_{npanels}.txt"
results = open(temp_path, "w")
results.write(f"===== Parameters =====\n"
              f"NACA: {var('naca_foil', 'str')}\n"
              f"AoA: {var('angle_of_attack', 'str')}\n"
              f"Panel number: {var('panel_number', 'str')}\n"
              f"V infinity: {var('v_infinity', 'str')}\n\n"
              f"===== Lift Coefficient =====\n"
              f"Code: {lift_coeff:.5f}\n"
              f"XFoil: {lift_coeff_xfoil:.5f}\n\n"
              f"===== Moment Coefficient =====\n"
              f"Code: {moment_coeff:.5f}\n"
              f"XFoil: {moment_coeff_xfoil:.5f}\n\n"
              f"===== Drag Coefficient =====\n"
              f"Code: {drag_coeff:.5f}\n"
              f"XFoil: {drag_coeff_xfoil:.5f}\n\n"
              f"===== Circulation =====\n"
              f"Gamma: {g:.5f}\n")
results.close()

if var("streamline_comp", "str") == "True":
    print(f"Results exported in \"{temp_path}\", computing streamlines...", end="\n")
else:
    print(f"Results exported in \"{temp_path}\", skipping streamlines...")

#streamlines

ngridx = var("grid_size_x", "int")
ngridy = var("grid_size_y", "int")
x_interval = [-1, 2]
y_interval = [-2.5, 2.5]
y_length = (max(yarr) - min(yarr))

gridy_array = np.linspace(y_interval[0], y_interval[1], ngridy)
gridx_array = x_interval[0] * np.ones(len(gridy_array))
starting = np.vstack((gridx_array.T, gridy_array.T)).T
xgrid = np.linspace(x_interval[0], x_interval[1], ngridx)
ygrid = np.linspace(y_interval[0], y_interval[1], ngridy)
xx, yy = np.meshgrid(xgrid, ygrid)

vx = np.zeros([ngridx, ngridy])
vy = np.zeros([ngridx, ngridy])

#create an airfoil closed surface
airfoil_path = mpl.path.Path(np.vstack((xarr.T, yarr.T)).T, closed=True)

#grid point velocities
if var("streamline_comp", "str") == "True":
    for m in tqdm(range(ngridx), desc="Streamlines"):
        for n in range(ngridy):
            xtemp = xx[m, n]
            ytemp = yy[m, n]
            mx, my = ijstreamline(xtemp, ytemp, xarr, yarr, phipanel, lenpanel, npanels)
            nx, ny = klstreamline(xtemp, ytemp, xarr, yarr, phipanel, lenpanel, npanels)
            #velocity if point inside airfoil = 0
            if airfoil_path.contains_points([(xtemp, ytemp)]):
                vx[m, n] = 0
                vy[m, n] = 0
            #velocities outside airfoil
            else:
                vx[m, n] = (vinf * np.cos(aoa) + sum(l * mx / (2 * np.pi))
                            + sum(- g * nx / (2 * np.pi)))
                vy[m, n] = (vinf * np.sin(aoa) + sum(l * my / (2 * np.pi))
                            + sum(- g * ny / (2 * np.pi)))

#velocity magnitudes
v_magnitude = np.sqrt(vx**2 + vy**2)  # Compute magnitude of velocity vector []
cp_grid = 1 - (v_magnitude / vinf)**2

print("Finished calculations, plotting...")

fig = plt.figure(3)
cp_scaled = np.absolute(cp * 0.15)
x_val = np.zeros(2)
y_val = np.zeros(2)
plt.plot([-1, 0], [-0.5, 0], "-", color="red", label="positive cp")
plt.plot([-1, 0], [-0.5, 0], "-", color="blue", label="negative cp")
for i in range(len(cp_scaled)):
    x_val[0] = xcpoint[i]
    x_val[1] = xcpoint[i] + cp_scaled[i] * np.cos(delta[i])
    y_val[0] = ycpoint[i]
    y_val[1] = ycpoint[i] + cp_scaled[i] * np.sin(delta[i])
    if cp[i] < 0:
        plt.plot(x_val, y_val, color="blue")
    elif cp[i] >= 0:
        plt.plot(x_val, y_val, color="red")
plt.fill(xarr, yarr, color="black", label=f"naca {param}")
plt.legend(loc='upper left', framealpha=1)
plt.grid()
fig.set_size_inches(width, width * ratio)
plt.xlim([-0.2, 1.05])
plt.ylim([-0.625, 0.625])
plt.savefig(f"plots/cp_{param}.{form}", dpi=dpi)
if var("cp_plot", "str") == "True":
    plt.show()
plt.close()

fig = plt.figure(4)
midxfoil = int(np.floor(len(x_cp_xfoil) / 2))
mid = int(np.floor(len(xarr) / 2))
if midxfoil != 0:
    plt.plot(x_cp_xfoil[midxfoil:], cp_xfoil[midxfoil:], color="gray", label="xfoil")
    plt.plot(x_cp_xfoil[:midxfoil], cp_xfoil[:midxfoil], color="gray")
plt.plot(xcpoint[:mid], cp[:mid], "x-", color="black", label="code")
plt.plot(xcpoint[mid:], cp[mid:], "x-", color="black")
plt.legend(loc='upper left', framealpha=1)
plt.grid()
fig.set_size_inches(width, width * ratio)
plt.xlim([-0.05, 1.05])
plt.ylim([-1.5, 1.5])
plt.savefig(f"plots/cps_{param}.{form}", dpi=dpi)
if var("cp_comparaison", "str") == "True":
    plt.show()
plt.close()

if var("streamline_comp", "str") == "True":
    fig = plt.figure(5)
    vn = np.sqrt(vx**2 + vy**2)
    np.seterr(under="ignore")
    plt.fill(xarr, yarr, color="black", zorder=10)
    plt.streamplot(xx, yy, vx, vy, color=vn, linewidth=0.75, density=100, arrowstyle="-", start_points=starting, cmap="cool_r")
    plt.colorbar(label="Velocity")
    plt.gca().set_aspect("equal")
    plt.xlim([-1, 2])
    plt.ylim([-0.5, 0.5])
    fig.set_size_inches(width, width * (1 / 3))
    plt.savefig(f"plots/streamlines_{param}.{form}", dpi=dpi)
    if var("streamline_plot", "str") == "True":
        plt.show()
    plt.close()

startb = tm.time()
lenmin = (startb - starta) // 60
lensec = round((startb - starta) - lenmin * 60, 0)
print(f"Process finished in {lenmin:.0f} min {lensec:.0f} seconds.")