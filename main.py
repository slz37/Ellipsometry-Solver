'''
Takes psi and delta data gathered from a PHE 101 ellipsometer and
fits them to fit and estimate optical properties/thickness of the layers
as desired. Note: thicknesses are in units of nm!

rho = r_p / r_s = tan(psi) * exp(i * delta)
https://en.wikipedia.org/wiki/Ellipsometry
'''

import os, sys
from math import *
import numpy as np
import pandas as pd
from plot import plot
from layers import layer
from solver import *

def filter_files(files, selection):
    file_list = []

    #Create list of files that are desired and return only them
    for i in selection:
        if str(i) + '.dat' in files:
            file_list.append(str(i) + '.dat')

    return file_list

def load_data(loc, angles = 'all'):
    phi = []
    psi = []
    delta = []
    
    #Load all .dat files from ellipsometer
    files = os.listdir(loc)

    #Filter for desired angles
    if 'all' or '' in angles:
        pass
    else:
        files = filter_files(files, angles)
    
    for i in files:
        if i.endswith('.dat'):
            data = open(loc + '/' + i)

            #Grab line with data for angle file
            for j, line in enumerate(data):
                if j == 5:
                    #Split and save data
                    data = line.split()
                    phi.append(data[1])
                    psi.append(data[2])
                    delta.append(data[3])
                
    #Angles in radians before returning
    phi = np.asarray(phi, dtype = np.float32) * pi / 180
    psi = np.asarray(psi, dtype = np.float32) * pi / 180
    delta = np.asarray(delta, dtype = np.float32) * pi / 180

    #Calculate rho
    rho = [np.real(np.tan(psi) * np.exp(1j * delta)), \
                     np.imag(np.tan(psi) * np.exp(1j * delta))]
    
    return phi, psi, delta, rho

#Testing
if __name__ == '__main__':
    #Get necessary data for fitting functions
    angles = 'all'

                        #Single layer test#
    phi, psi, delta, rho = load_data('4.19.16 fresh', angles)

    #Setup Layers
    top = layer('MgB2', 40.0781912)
    bot = layer('SiC', 390000)

    #Test different combinations of fit variables
    test_vars = ['d']
    if 'd' in test_vars:
        top.set_fit_variables([['d', 45]])
    if 'k' in test_vars:
        top.set_fit_variables([['k', 2.123]])
    if 'n' in test_vars:
        top.set_fit_variables([['n', 1.017]])
    layers = [bot, top]

    #Should always be setup substrate->layer 1->...->layer n
    fit_vars, cov = solve(phi, rho, layers)

                        #Double layer test#
    phi, psi, delta, rho = load_data('4.14.16 oxide', angles)

    #Setup layers
    top = layer('MgO', 2)
    mid = layer('MgB2', 45)
    bot = layer('SiC', 250000)
    
    #Test different combinations of fit variables
    test_vars = ['d']
    if 'd' in test_vars:
        mid.set_fit_variables([['d', 46]])
    if 'k' in test_vars:
        mid.set_fit_variables([['k', 2.123]])
    if 'n' in test_vars:
        mid.set_fit_variables([['n', 1.017]])
    layers = [bot, top]

    layers = [bot, mid, top]
    
    #Fit to function - put these in order from bot - top
    fit_vars, cov = solve(phi, rho, layers)

    #Plot fitted functions with data points
    plot(phi, rho, fit_vars, cov, layers)    
