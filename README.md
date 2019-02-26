# Ellipsometry-Solver
This project was started as part of my sophomore year undergraduate research on superconducting thin films. The PHE 101 ellipsometer program was not functioning correctly, so the motivation was to develop a new program that would read the data gathered and output the desired properties. This code was originally developed in Matlab; however, I am currently porting it to Python since I would like to clean up the code and develop my Python skills further.

Numerically solves the ellipsometry equations to find the thickness, n, and/or k values of thin films. This program takes in data from a PHE 101 Ellipsometer containing values for the laser wavelength, angle of incident, rho, and psi values. It then uses these to obtain rho values from the equation <img src="https://latex.codecogs.com/svg.latex?\inline&space;\rho&space;=&space;\frac{r_{p}}{r_{s}}&space;=&space;tan(\psi)&space;*&space;e^{i&space;\delta}" title="\rho = \frac{r_{p}}{r_{s}} = tan(\psi) * e^{i \delta}" />.

These values can be used with equations from https://jameslandry.files.wordpress.com/2012/02/09-appendix-b.pdf to numerically solve for the selected parameters. It will then output the value, MSE, and a plot of the fit along with the data gathered.

Current status: Should work for single layer films (substrate + top layer), e.g. MgB2 on SiC, for fitting any combination of n, k, d of the top layer. Currently working on implementing this for double layer films (substrate + mid + top layer), e.g. MgO on MgB2 on SiC, and possibly generalizing to any n-layer film.

A single layer example is provided. The test.py script inside 4.19.16 fresh folder should be run without any changes necessary. A double layer example is also provided despite the code not being able to numerically solve this scenario yet. Output should look like below:

```
MgB2:
	k:  2.108450209193739 ± 0.061766547404060325
	d:  39.95075039115975 ± 2.355272768402803

	MSE:  0.00018533485
```
<img align="center" src="./4.19.16 fresh/4.19.16-fresh-example-plot.png" alt="Example single layer fit plot" title="Single Layer Plot" hspace="20"/>
