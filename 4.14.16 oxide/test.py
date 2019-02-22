'''
Example script for MgB2 deposition on SiC and measured at a later time after an
oxide layer formed on top. This is an example of a double layer film.
MgO thickness of 2 nm was determined through estimation. MgB2 thickness of
46 nm was determined by AFM to confirm that this code works.
Read info.txt for more details.
'''

import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from main import *

#Select files to read data from
#angles = ['40', '50', '55', '60']
angles = 'all'

#Single layer test
phi, psi, delta, rho = load_data('.', angles)

#Change dir to read Material_Information.txt
#Not necessary if in same dir as main.py
os.chdir("../")

#Setup Layers - provide random thickness if unknown
top = layer('MgO', 2)
mid = layer('MgB2', 46)
bot = layer('SiC', 250000)

#Test different combinations of fit variables
#Future functionality for fitting layers other than
#main will be added
test_vars = ['d']
if 'd' in test_vars:
    mid.set_fit_variables([['d', 45]])
if 'k' in test_vars:
    mid.set_fit_variables([['k', 2.123]])
if 'n' in test_vars:
    mid.set_fit_variables([['n', 1.017]])

#Should always be setup substrate->layer 1->...->layer n
layers = [bot, mid, top]

fit_vars, cov = solve(phi, rho, layers)

#Plot fitted functions with data points
plot(phi, rho, fit_vars, cov, layers)   
