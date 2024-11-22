import numpy as np
import math as m

np.seterr('raise')

def ijstreamline(xgrid, ygrid, xarr, yarr, phipanel, lenpanel, npanels):
    M_x = np.zeros(npanels)
    M_y = np.zeros(npanels)
    for i in range(npanels):
        A = - (xgrid - xarr[i]) * np.cos(phipanel[i]) - (ygrid - yarr[i]) * np.sin(phipanel[i])
        B = (xgrid - xarr[i])**2 + (ygrid - yarr[i])**2
        C_x = - np.cos(phipanel[i])
        C_y = - np.sin(phipanel[i])
        D_x = xgrid - xarr[i]
        D_y = ygrid - yarr[i]
        E = np.sqrt(B - A**2)
        #manage value problems
        if E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E):
            M_x[i] = 0
            M_y[i] = 0
        #compute integrals
        else:
            temp1 = 0.5 * C_x * np.log((lenpanel[i]**2 + 2 * A * lenpanel[i] + B) / B)
            temp2 = ((D_x - A * C_x) / E) * (m.atan2((lenpanel[i] + A), E) - m.atan2(A, E))
            M_x[i] = temp1 + temp2
            temp1 = 0.5 * C_y * np.log((lenpanel[i]**2 + 2 * A * lenpanel[i] + B) / B)
            temp2 = ((D_y - A * C_y) / E) * (m.atan2((lenpanel[i] + A), E) - m.atan2(A, E))
            M_y[i] = temp1 + temp2
        #manage value problems
        if np.iscomplex(M_x[i]) or np.isnan(M_x[i]) or np.isinf(M_x[i]):
            M_x[i] = 0
        if np.iscomplex(M_y[i]) or np.isnan(M_y[i]) or np.isinf(M_y[i]):
            M_y[i] = 0

    return M_x, M_y

def klstreamline(xgrid, ygrid, xarr, yarr, phipanel, lenpanel, npanels):
    N_x = np.zeros(npanels)
    N_y = np.zeros(npanels)
    for i in range(npanels):
        A = - (xgrid - xarr[i]) * np.cos(phipanel[i]) - (ygrid - yarr[i]) * np.sin(phipanel[i])
        B = (xgrid - xarr[i])**2 + (ygrid - yarr[i])**2
        C_x = np.sin(phipanel[i])
        C_y = - np.cos(phipanel[i])
        D_x = - (ygrid - yarr[i])
        D_y = xgrid - xarr[i]
        E = np.sqrt(B - A**2)
        #manage value problems
        if E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E):
            N_x[i] = 0
            N_y[i] = 0
        #compute integrals
        else:
            temp1 = 0.5 * C_x * np.log((lenpanel[i]**2 + 2 * A * lenpanel[i] + B) / B)
            temp2 = ((D_x - A * C_x) / E) * (m.atan2((lenpanel[i] + A), E) - m.atan2(A, E))
            N_x[i] = temp1 + temp2
            temp1 = 0.5 * C_y * np.log((lenpanel[i] ** 2 + 2 * A * lenpanel[i] + B) / B)
            temp2 = ((D_y - A * C_y) / E) * (m.atan2((lenpanel[i] + A), E) - m.atan2(A, E))
            N_y[i] = temp1 + temp2
        #manage value problems
        if np.iscomplex(N_x[i]) or np.isnan(N_x[i]) or np.isinf(N_x[i]):
            N_x[i] = 0
        if np.iscomplex(N_y[i]) or np.isnan(N_y[i]) or np.isinf(N_y[i]):
            N_y[i] = 0

    return N_x, N_y
