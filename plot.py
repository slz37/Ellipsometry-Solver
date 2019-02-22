'''
Plots two sets of data: real and complex. Each set will
have a fit line and a set of points gathered from the
ellipsometer. Also prints results of fitting function.
'''
import matplotlib.pyplot as plt
import numpy as np
from math import pi
from solver import get_rho_values
from sklearn.metrics import mean_squared_error

def print_results(fit_vars, phi, rho_act, layers):
    #Get predicted rho values at specific phi
    rho_pred = get_rho_values(phi, fit_vars, layers)

    #Need to split real and imag
    rho_pred = np.split(rho_pred, 2)
    MSE = mean_squared_error(rho_act, rho_pred)

    #Iterate over each layer and print fit values
    #if they were fit
    for layer in layers:
        #Skip if nothing was fit
        if layer.fit_n or layer.fit_k or layer.fit_d:
            print(layer.name + ':')
        else:
            continue

        #Print fit variables - always organized as [n, k, d]
        if layer.fit_n:
            print('\tn: ', fit_vars[0])
            
        if layer.fit_n and layer.fit_k:
            print('\tk: ', fit_vars[1])
        elif layer.fit_k:
            print('\tk: ', fit_vars[0])

        if layer.fit_n and layer.fit_k and layer.fit_d:
            print('\td: ', fit_vars[2])
        elif (layer.fit_n or layer.fit_k) and layer.fit_d:
            print('\td: ', fit_vars[1])
        elif layer.fit_d:
            print('\td: ', fit_vars[0])
    print('\n\tMSE: ', MSE)

def plot(phi, rho, fit_vars, layers):
    #Values for fitted functions
    phi_fit = np.linspace(min(phi * 180 / pi), max(phi * 180 / pi))
    rho_fit = get_rho_values(phi_fit * pi / 180, fit_vars, layers)

    print_results(fit_vars, phi, rho, layers)

    #Plot details
    plt.plot(phi * 180 / np.pi, rho[0], 'b.')
    plt.plot(phi * 180 / np.pi, rho[1], 'r.')
    plt.plot(phi_fit, rho_fit[:50], 'b')
    plt.plot(phi_fit, rho_fit[50:100], 'r')
    plt.xlabel('Phi');
    plt.ylabel('Rho');
    plt.legend(['Real Data', 'Complex Data', 'Real Fitted Curve', \
        'Complex Fitted Curve'], loc = 'best')
    plt.show()

