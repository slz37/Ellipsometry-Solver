'''
We have rho values calculated from delta and psi, and
rho = r_p / r_s

r_p and r_s can be put into terms of n, k, d of all layers as seen here:
https://jameslandry.files.wordpress.com/2012/02/09-appendix-b.pdf

Then solve numerically for specific n, k, d of layers we want.
'''

from math import *
import numpy as np
from lmfit import Model
from scipy.optimize import curve_fit

def solve(phi, rho, layers):
    #Fit parameters
    xdata = phi
    ydata = np.append(rho[0], rho[1])
    guess = [row[1] for row in layers[1].get_fit_variables()]
    
    #Fit function
    funct = layer_functions(layers)

    #Ensure proper number of layers
    if len(layers) == 2:
        coeffs, coeffs_cov = curve_fit(funct.single_layer,
                                       xdata,
                                       ydata,
                                       p0 = guess)
    elif len(layers) == 3:
        coeffs, coeffs_cov = curve_fit(funct.double_layer,
                                       xdata,
                                       ydata,
                                       p0 = guess)
    else:
        raise Exception('You should have either 2 or 3 layers,' + \
                        ' but there are {} layers.'.format(len(layers)))

    return coeffs, coeffs_cov

def get_rho_values(x, fit_vars, layers):
    #Function class instantiation
    funct = layer_functions(layers)
    
    #Ensure proper number of layers and grab fitted rho values
    if len(layers) == 2:
        y_vals = funct.single_layer(x, *fit_vars)
    elif len(layers) == 3:
        y_vals = funct.double_layer(x, *fit_vars)
    else:
        raise Exception('You should have either 2 or 3 layers,' + \
                        ' but there are {} layers.'.format(len(layers)))
                        
    return y_vals

class layer_functions():
    def __init__(self, layers):
        #Constants
        self.lamb = 632.8        #Wavelength of laser
        self.layers = layers
        self.sub_n = layers[0].n     #Substrate Layer n
        self.sub_k = layers[0].k     #Substrate Layer k
        
        if len(layers) == 3:
            self.top_n = layers[2].n     #Top Layer n
            self.top_k = layers[2].k     #Top Layer k
            self.top_d = layers[2].d     #Top layer d
        
    def single_layer(self, x, *guess):
        '''
                            Table of Layer Information
                    Layer   - Description       - Variables
                    Layer 0 - Air               - x, (n = 1, k = 0 both not in equs)
                    Layer 1 - Film Layer        - theta1, n, k, d
                    Layer 2 - Substrate Layer   - theta2, n, k
        '''
        #Possible Variables to Fit  - always organized as [n, k, d]
        if self.layers[1].fit_n:
            self.n = guess[0]
        else:
            self.n = self.layers[1].n

        if self.layers[1].fit_n and self.layers[1].fit_k:
            self.k = guess[1]
        elif self.layers[1].fit_k:
            self.k = guess[0]
        else:
            self.k = self.layers[1].k

        if (self.layers[1].fit_n and self.layers[1].fit_k) and self.layers[1].fit_d:
            self.d = guess[2]
        elif (self.layers[1].fit_n or self.layers[1].fit_k) and self.layers[1].fit_d:
            self.d = guess[1]
        elif self.layers[1].fit_d:
            self.d = guess[0]
        else:
            self.d = self.layers[1].d
            
        #Theta Equations
        self.theta1 = np.arcsin(np.sin(x) / (self.n - 1j * self.k))
        self.theta2 = np.arcsin(np.sin(x) / (self.sub_n - 1j * self.sub_k))

        #Beta Equation = -2gamma
        self.beta = -4 * np.pi * self.d * (self.n - 1j * self.k) * np.cos(self.theta1) / self.lamb

        #P-Polarized Equations
        self.r01p = ((self.n - 1j * self.k) * np.cos(x) - np.cos(self.theta1)) / \
                    ((self.n - 1j * self.k) * np.cos(x) + np.cos(self.theta1))
        self.r12p = ((self.sub_n - 1j * self.sub_k) * np.cos(self.theta1) - (self.n - 1j * self.k) * np.cos(self.theta2)) / \
                    ((self.sub_n - 1j * self.sub_k) * np.cos(self.theta1) + (self.n - 1j * self.k) * np.cos(self.theta2))
        self.rptot = (self.r01p + self.r12p * np.exp(1j * self.beta)) / \
                     (1 + self.r01p * self.r12p * np.exp(1j * self.beta))

        #S-Polarized Equations
        self.r01s = (np.cos(x) - (self.n - 1j * self.k) * np.cos(self.theta1)) / \
                    (np.cos(x) + (self.n - 1j * self.k) * np.cos(self.theta1))
        self.r12s = ((self.n - 1j * self.k) * np.cos(self.theta1) - (self.sub_n - 1j * self.sub_k) * np.cos(self.theta2)) / \
                    ((self.n - 1j * self.k) * np.cos(self.theta1) + (self.sub_n - 1j * self.sub_k) * np.cos(self.theta2))
        self.rstot = (self.r01s + self.r12s * np.exp(1j * self.beta)) / \
                     (1 + self.r01s * self.r12s * np.exp(1j * self.beta))

        #Split Function Into Real and Complex Components
        return np.append(np.real(self.rptot / self.rstot), np.imag(self.rptot / self.rstot))

    def double_layer(self, x, *guess):
        '''
                            Table of Layer Information
                    Layer   - Description       - Variables
                    Layer 0 - Air               - x, (n = 1, k = 0 both not in eqs)
                    Layer 1 - First Film Layer  - theta1, n, k, d1
                    Layer 2 - Second Film Layer - theta2, m, l, d2
                    Layer 3 - Substrate Layer   - theta3, n, k
        '''
        #Possible Variables to Fit  - always organized as [n, k, d]
        if self.layers[1].fit_n:
            self.n = guess[0]
        else:
            self.n = self.layers[1].n

        if self.layers[1].fit_n and self.layers[1].fit_k:
            self.k = guess[1]
        elif self.layers[1].fit_k:
            self.k = guess[0]
        else:
            self.k = self.layers[1].k

        if (self.layers[1].fit_n and self.layers[1].fit_k) and self.layers[1].fit_d:
            self.d = guess[2]
        elif (self.layers[1].fit_n or self.layers[1].fit_k) and self.layers[1].fit_d:
            self.d = guess[1]
        elif self.layers[1].fit_d:
            self.d = guess[0]
        else:
            self.d = self.layers[1].d

        #Theta Equations
        self.theta1 = np.arcsin(np.sin(x) / (self.top_n - 1j * self.top_k))
        self.theta2 = np.arcsin(np.sin(x) / (self.n - 1j * self.k))
        self.theta3 = np.arcsin(np.sin(x) / (self.sub_n - 1j * self.sub_n))

        #Beta Equations
        self.beta1 = -4 * np.pi * self.top_d * (self.top_n - 1j * self.top_k) * np.cos(self.theta1) / self.lamb
        self.beta2 = -4 * np.pi * self.d * (self.n - 1j * self.k) * np.cos(self.theta2) / self.lamb

        #P-Polarized Equations
        self.r01p = ((self.top_n - 1j * self.top_k) * np.cos(x) - np.cos(self.theta1)) / \
                    ((self.top_n - 1j * self.top_k) * np.cos(x) + np.cos(self.theta1))
        self.r12p = ((self.n - 1j * self.k) * np.cos(self.theta1) - (self.top_n - 1j * self.top_k) * np.cos(self.theta2)) / \
                    ((self.n - 1j * self.k) * np.cos(self.theta1) + (self.top_n - 1j * self.top_k) * np.cos(self.theta2))
        self.r23p = ((self.sub_n - 1j * self.sub_k) * np.cos(self.theta2) - (self.n - 1j * self.k) * np.cos(self.theta3)) / \
                    ((self.sub_n - 1j * self.sub_k) * np.cos(self.theta2) + (self.n - 1j * self.k) * np.cos(self.theta3))
        self.r123p = (self.r12p + self.r23p * np.exp(1j * self.beta2)) / \
                     (1 + self.r12p * self.r23p * np.exp(1j * self.beta2))
        self.rptot = (self.r01p + self.r123p * np.exp(1j * self.beta1)) / \
                     (1 + self.r01p * self.r123p * np.exp(1j * self.beta1))

        #S-Polarized Equations
        self.r01s = (np.cos(x) - (self.top_n - 1j * self.top_k) * np.cos(self.theta1)) / \
                    (np.cos(x) + (self.top_n - 1j * self.top_k) * np.cos(self.theta1))
        self.r12s = ((self.top_n - 1j * self.top_k) * np.cos(self.theta1) - (self.n - 1j * self.k) * np.cos(self.theta2)) / \
                    ((self.top_n - 1j * self.top_k) * np.cos(self.theta1) + (self.n - 1j * self.k) * np.cos(self.theta2))
        self.r23s = ((self.n - 1j * self.k) * np.cos(self.theta2) - (self.sub_n - 1j * self.sub_k) * np.cos(self.theta3)) / \
                    ((self.n - 1j * self.k) * np.cos(self.theta2) + (self.sub_n - 1j * self.sub_k) * np.cos(self.theta3))
        self.r123s = (self.r12s + self.r23s * np.exp(1j * self.beta2)) / \
                     (1 + self.r12s * self.r23s * np.exp(1j * self.beta2))
        self.rstot = (self.r01s + self.r123s * np.exp(1j * self.beta1)) / \
                     (1 + self.r01s * self.r123s * np.exp(1j * self.beta1))

        #Split Function Into Real and Complex Components
        return np.append(np.real(self.rptot / self.rstot), np.imag(self.rptot / self.rstot))
