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
from lmfit import model
from scipy.optimize import curve_fit

def solve(psi, rho, layers):
    #Ensure proper number of layers
    if len(layers) != 2:
        raise Exception("You should have exactly two layers: substrate and film.")
    
    #Fit parameters
    xdata = psi
    ydata = rho
    guess = [row[1] for row in layers[1].get_fit_variables()]
    
    #Fit function
    popt, pcov = curve_fit(single_layer, xdata, ydata, p0 = guess)

def single_layer(x, layers):
    print(layers)
    #Constants
    lamb = 632.8        #Wavelength of laser
    m = layers[0].n     #Substrate Layer n
    l = layers[0].k     #Substrate Layer k

    #Possible Variables to Fit
    if layers[1].fit_n:
        n = layers[1].n_guess 
    else:
        n = layers[1].n

    if layers[1].fit_k:
        k = layers[1].k_guess
    else:
        k = layers[1].k
        
    if layers[1].fit_d:
        d = layers[1].d_guess
    else:
        d = layers[1].d

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
    rptot = (r01p + r12p * exp(1j * beta)) / (1 + r01p * r12p * \
        exp(1j * beta))

    #S-Polarized Equations
    r01s = (np.cos(x) - (n - 1j * k) * np.cos(theta1)) / (np.cos(x) + (n - 1j * \
        k) * np.cos(theta1))
    r12s = ((n - 1j * k) * np.cos(theta1) - (m - 1j * l) * np.cos(theta2)) / \
        ((n - 1j * k) * np.cos(theta1) + (m - 1j * l) * np.cos(theta2))
    rstot = (r01s + r12s * exp(1j * beta)) / (1 + r01s * r12s * \
        exp(1j * beta))

    #Split Function Into Real and Complex Components
    yout = [real(rptot / rstot), imag(rptot / rstot)]


'''
                            Table of Layer Information
                    Layer   - Description       - Variables
                    Layer 0 - Air               - x, (n, k = 1)
                    Layer 1 - First Film Layer  - theta1, n, k, d1
                    Layer 2 - Second Film Layer - theta2, m, l, d2
                    Layer 3 - Substrate Layer   - theta3, o, j


def double_layer():
    #Constants
    lamb = 632.8;     #Wavelength
    o = t(1);        #Substrate n
    j = t(2);        #Substrate k
    n = t(3);        #First Film Layer n
    k = t(4);        #First Film Layer k
    d1 = t(5);       #First Film Layer Thickness

    #Possible Variables to Fit
    if t(7) ~= -1
        m = t(7);            #Second Film Layer n - not fit
    else
        m = a(2);               #Second Film Layer n - fit
    

    if t(6) ~= -1
        l = t(6);            #Second Film Layer k - not fit
    else
        l = a(1);               #Second Film Layer k - fit

    if t(8) ~= -1
        d2 = t(8);           #Second Film Layer d - not fit
    else
        d2 = a(3);              #Second Film Layer d - fit

    #Theta Equations
    theta1 = arcsin(np.sin(x) / (n - 1j * k));
    theta2 = arcsin(np.sin(x) / (m - 1j * l));
    theta3 = arcsin(np.sin(x) / (o - 1j * j));

    #Beta Equations
    beta1 = -4 * np.pi * d1 * (n - 1j * k) * np.cos(theta1) / lamb;
    beta2 = -4 * np.pi * d2 * (m - 1j * l) * np.cos(theta2) / lamb;

    #P-Polarized Equations
    r01p = ((n - 1j * k) * np.cos(x) - np.cos(theta1)) / ((n - 1j * k) * \
        np.cos(x) + np.cos(theta1));
    r12p = ((m - 1j * l) * np.cos(theta1) - (n - 1j * k) * \
        np.cos(theta2)) / ((m - 1j * l) * np.cos(theta1) + (n - 1j * k) * \
        np.cos(theta2));
    r23p = ((o - 1j * j) * np.cos(theta2) - (m - 1j * l) * \
        np.cos(theta3)) / ((o - 1j * j) * np.cos(theta2) + (m - 1j * l) * \
        np.cos(theta3));
    r123p = (r12p + r23p * exp(1j * beta2)) / (1 + r12p * r23p * \
        exp(1j * beta2));
    rptot = (r01p + r123p * exp(1j * beta1)) / (1 + r01p * r123p * \
        exp(1j * beta1));

    #S-Polarized Equations
    r01s = (np.cos(x) - (n - 1j * k) * np.cos(theta1)) / (np.cos(x) + (n - 1j * \
        k) * np.cos(theta1));
    r12s = ((n - 1j * k) * np.cos(theta1) - (m - 1j * l) * np.cos(theta2)) / \
        ((n - 1j * k) * np.cos(theta1) + (m - 1j * l) * np.cos(theta2));
    r23s = ((m - 1j * l) * np.cos(theta2) - (o - 1j * j) * np.cos(theta3)) / \
        ((m - 1j * l) * np.cos(theta2) + (o - 1j * j) * np.cos(theta3));
    r123s = (r12s + r23s * exp(1j * beta2)) / (1 + r12s * r23s * \
        exp(1j * beta2));
    rstot = (r01s + r123s * exp(1j * beta1)) / (1 + r01s * r123s * \
        exp(1j * beta1));

    #Split Function Into Real and Complex Components
    yout = [real(rptot / rstot), imag(rptot / rstot)];
'''
