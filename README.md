# Ellipsometry-Solver
Numerically solves the ellipsometry equations to find the thickness, n, and/or k values of thin films. This program takes in data from a PHE 101 Ellipsometer containing values for the laser wavelength, angle of incident, rho, and psi values. It then uses these to obtain rho values from the equation <img src="https://latex.codecogs.com/svg.latex?\inline&space;\rho&space;=&space;\frac{r_{p}}{r_{s}}&space;=&space;tan(\psi)&space;*&space;e^{i&space;\delta}" title="\rho = \frac{r_{p}}{r_{s}} = tan(\psi) * e^{i \delta}" .

These values can be used with equations from https://jameslandry.files.wordpress.com/2012/02/09-appendix-b.pdf to numerically solve for the selected parameters. It will then output the value, MSE, and a plot of the fit along with the data gathered.
