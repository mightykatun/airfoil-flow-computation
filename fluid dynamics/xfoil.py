import os
import sys
from imports import var
import numpy as np

def xfoil(param:int):
    param = str(param)

    airfoil_save = "xfoil/" + param + ".txt"
    airfoil_cp = "xfoil/cp_" + param + ".txt"
    airfoil_pol = "xfoil/pol_" + param + ".txt"
    for i in [airfoil_save, airfoil_cp, airfoil_pol]:
        if os.path.exists(i):
            os.remove(i)
    # generate config file
    temp = open('xfoil/temp.inp', "w")
    temp.write("PLOP\n")
    temp.write("G F\n")
    temp.write("\n")
    temp.write(f"LOAD foils/naca_{param}.dat" + "\n")
    temp.write("PPAR\n")
    temp.write("N " + var("panel_number", "str") + "\n")
    temp.write("P " + var("panel_bunching", "str") + "\n")
    temp.write("T " + var("te/le_density", "str") + "\n")
    temp.write("R " + var("panel_density", "str") + "\n")
    temp.write("XT " + var("top_y/c_lim", "str") + "\n")
    temp.write("XB " + var("bottom_y/c_lim", "str") + "\n")
    temp.write("\n")
    temp.write("\n")
    temp.write("PSAV " + airfoil_save + "\n")
    temp.write("OPER\n")
    temp.write("Pacc 1 \n")
    temp.write("\n\n")
    temp.write("Alfa " + var("angle_of_attack", "str") + "\n")
    temp.write("CPWR " + airfoil_cp + "\n")
    temp.write("PWRT\n")
    temp.write(airfoil_pol + "\n")
    if os.path.exists(airfoil_pol):
        temp.write("y \n")
    temp.write("\n")
    temp.write("QUIT\n")
    temp.close()

    ## run xfoil
    if 'linux' in sys.platform:
        os.system("wine xfoil.exe < xfoil/temp.inp > xfoil/temp.out")
    elif 'win' in sys.platform:
        os.system("xfoil.exe < xfoil/temp.inp > xfoil/temp.out")
    ##

    data_cp = np.loadtxt(airfoil_cp, skiprows=3)
    x_cp = data_cp[:,0]
    y_cp = data_cp[:,1]
    cp = data_cp[:,2]
    data_pol = np.loadtxt(airfoil_pol, skiprows=12)
    lift_coeff = data_pol[1]
    drag_coeff = data_pol[2]
    moment_coeff = data_pol[4]
    if os.path.exists("xfoil/temp.inp"):
        os.remove("xfoil/temp.inp")
    if os.path.exists("xfoil/temp.out"):
        os.remove("xfoil/temp.out")
    return x_cp, y_cp, cp, lift_coeff, drag_coeff, moment_coeff
