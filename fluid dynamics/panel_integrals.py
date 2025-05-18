import numpy as np
import math as m
from tqdm import tqdm

np.seterr("raise")

def ijintegral(xarr, yarr, xcpoint, ycpoint, phipanel, lenpanel, npanels):
    I = np.zeros([npanels, npanels])
    J = np.zeros([npanels, npanels])
    for k in range(npanels):
        for l in range(npanels):
            if (l != k):
                A = - (xcpoint[k] - xarr[l]) * np.cos(phipanel[l]) - (ycpoint[k] - yarr[l]) * np.sin(phipanel[l])
                B = (xcpoint[k] - xarr[l])**2 + (ycpoint[k] - yarr[l])**2
                C_n = np.sin(phipanel[k] - phipanel[l])
                C_t = - np.cos(phipanel[k] - phipanel[l])
                D_n = - (xcpoint[k] - xarr[l]) * np.sin(phipanel[k]) + (ycpoint[k] - yarr[l]) * np.cos(phipanel[k])
                D_t = (xcpoint[k] - xarr[l]) * np.cos(phipanel[k]) + (ycpoint[k] - yarr[l]) * np.sin(phipanel[k])
                E = np.sqrt(B - A**2)
                #manage value problems
                if E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E):
                    I[k, l] = 0
                    J[k, l] = 0
                #compute integrals
                else:
                    temp1 = 0.5 * C_n * np.log((lenpanel[l]**2 + 2 * A * lenpanel[l] + B) / B)
                    temp2 = ((D_n - A * C_n) / E) * (m.atan2((lenpanel[l] + A), E) - m.atan2(A, E))
                    I[k, l] = temp1 + temp2
                    temp1 = 0.5 * C_t * np.log((lenpanel[l] ** 2 + 2 * A * lenpanel[l] + B) / B)
                    temp2 = ((D_t - A * C_t) / E) * (m.atan2((lenpanel[l] + A), E) - m.atan2(A, E))
                    J[k, l] = temp1 + temp2
            #manage value problems
            if (np.iscomplex(I[k, l]) or np.isnan(I[k, l]) or np.isinf(I[k, l])):
                I[k, l] = 0
            if (np.iscomplex(J[k, l]) or np.isnan(J[k, l]) or np.isinf(J[k, l])):
                J[k, l] = 0
    return I, J

def klintegral(xarr, yarr, xcpoint, ycpoint, phipanel, lenpanel, npanels):
    K = np.zeros([npanels, npanels])
    L = np.zeros([npanels, npanels])
    for k in tqdm(range(npanels), desc="  Integrals"):
        for l in range(npanels):
            if (l != k):
                A = - (xcpoint[k] - xarr[l]) * np.cos(phipanel[l]) - (ycpoint[k] - yarr[l]) * np.sin(phipanel[l])
                B = (xcpoint[k] - xarr[l])**2 + (ycpoint[k] - yarr[l])**2
                C_n = - np.cos(phipanel[k] - phipanel[l])
                C_t = np.sin(phipanel[k] - phipanel[l])
                D_n = (xcpoint[k] - xarr[l]) * np.cos(phipanel[k]) + (ycpoint[k] - yarr[l]) * np.sin(phipanel[k])
                D_t = (xcpoint[k] - xarr[l]) * np.sin(phipanel[k]) + (ycpoint[k] - yarr[l]) * np.cos(phipanel[k])
                E = np.sqrt(B - A**2)
                #manage value problems
                if E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E):
                    K[k, l] = 0
                    L[k, l] = 0
                #compute integrals
                else:
                    temp1 = 0.5 * C_n * np.log((lenpanel[l]**2 + 2 * A * lenpanel[l] + B) / B)
                    temp2 = ((D_n - A * C_n) / E) * (m.atan2((lenpanel[l] + A), E) - m.atan2(A, E))
                    K[k, l] = temp1 + temp2
                    temp1 = 0.5 * C_t * np.log((lenpanel[l] ** 2 + 2 * A * lenpanel[l] + B) / B)
                    temp2 = ((D_t - A * C_t) / E) * (m.atan2((lenpanel[l] + A), E) - m.atan2(A, E))
                    L[k, l] = temp1 + temp2
            #manage value problems
            if (np.iscomplex(K[k, l]) or np.isnan(K[k, l]) or np.isinf(K[k, l])):
                K[k, l] = 0
            if (np.iscomplex(L[k, l]) or np.isnan(L[k, l]) or np.isinf(L[k, l])):
                L[k, l] = 0
    return K, L

