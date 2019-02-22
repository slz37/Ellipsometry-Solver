'''
                            Table of Layer Information
                    Layer   - Description       - Variables
                    Layer 0 - Air               - x, (n = 1, k = 0)
                    Layer 1 - Film Layer        - theta1, n, k, d
                    Layer 2 - Substrate Layer   - theta2, m, l
'''

'''
We have rho values calculated from delta and psi, and
rho = r_p / r_s

r_p and r_s can be put into terms of n, k, d of all layers as seen here:
https://jameslandry.files.wordpress.com/2012/02/09-appendix-b.pdf

Then solve numerically for specific n, k, d of layers we want.
'''

import numpy as np
from lmfit import Model, Parameters, minimize
from scipy.optimize import curve_fit

def solve(psi, rho, layers):
    #Ensure proper number of layers
    if len(layers) != 2:
        raise Exception("You should have exactly two layers: substrate and film.")
    
    #Fit parameters
    xdata = psi
    ydata = rho
    guess = [row[1] for row in layers[1].get_fit_variables()]

    #Setup model and set variables constant
    single_layer_model = Model(single_layer)
    params = single_layer_model.make_params(lamb = 632.8, m = layers[0].n, \
                                            l = layers[0].k)
    params['lamb'].vary = False
    params['m'].vary = False
    params['l'].vary = False

    #Possible Variables to Fit
    if layers[1].fit_n:
        params.add('n', layers[1].n_guess)
    else:
        n = params.add('n', layers[1].n)
        params['n'].vary = False

    if layers[1].fit_k:
        params.add('k', layers[1].k_guess)
    else:
        n = params.add('k',layers[1].k)
        params['k'].vary = False
        
    if layers[1].fit_d:
        params.add('d', layers[1].d_guess)
    else:
        n = params.add('d', layers[1].d)
        params['d'].vary = False
    
    #Fit function
    results = single_layer_model.fit(ydata, params, x = xdata)
    return results

def single_layer(x, lamb, m, l, n, k, d):
    #Theta Equations
    theta1 = np.arcsin(np.sin(x) / (n - 1j * k))
    theta2 = np.arcsin(np.sin(x) / (m - 1j * l))

    #Beta Equation
    beta = -4 * np.pi * d * (n - 1j * k) * np.cos(theta1) / lamb

    #P-Polarized Equations
    r01p = ((n - 1j * k) * np.cos(x) - np.cos(theta1)) / ((n - 1j * k) \
        * np.cos(x) + np.cos(theta1))
    r12p = ((m - 1j * l) * np.cos(theta1) - (n - 1j * k) * np.cos(theta2)) \
        / ((m - 1j * l) * np.cos(theta1) + (n - 1j * k) * np.cos(theta2))
    rptot = (r01p + r12p * np.exp(1j * beta)) / (1 + r01p * r12p * \
        np.exp(1j * beta))

    #S-Polarized Equations
    r01s = (np.cos(x) - (n - 1j * k) * np.cos(theta1)) / (np.cos(x) + (n - 1j * \
        k) * np.cos(theta1))
    r12s = ((n - 1j * k) * np.cos(theta1) - (m - 1j * l) * np.cos(theta2)) / \
        ((n - 1j * k) * np.cos(theta1) + (m - 1j * l) * np.cos(theta2))
    rstot = (r01s + r12s * np.exp(1j * beta)) / (1 + r01s * r12s * \
        np.exp(1j * beta))

    #Split Function Into Real and Complex Components
    return rptot / rstot
