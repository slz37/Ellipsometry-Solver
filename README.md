# Ellipsometry-Solver
Numerically solves the ellipsometry equations to find the thickness, n, and/or k values of thin films. This program takes in data from a PHE 101 Ellipsometer containing values for the laser wavelength, angle of incident, rho, and psi values. It then uses these to obtain rho values from the equation <img src="https://latex.codecogs.com/svg.latex?\inline&space;\rho&space;=&space;\frac{r_{p}}{r_{s}}&space;=&space;tan(\psi)&space;*&space;e^{i&space;\delta}" title="\rho = \frac{r_{p}}{r_{s}} = tan(\psi) * e^{i \delta}" />.

These values can be used with equations from https://jameslandry.files.wordpress.com/2012/02/09-appendix-b.pdf to numerically solve for the selected parameters. It will then output the value, MSE, and a plot of the fit along with the data gathered.

A single layer example is provided. The test.py script inside 4.19.16 fresh folder should be run without any changes necessary. A double layer example will be provided once functionality is complete for that case.

Current status: Should work for single layer films (substrate + top layer), e.g. MgB2 on SiC, for fitting any combination of n, k, d of the top layer. Currently working on implementing this for double layer films (substrate + mid + top layer), e.g. MgO on MgB2 on SiC, and possibly generalizing to any n-layer film.
