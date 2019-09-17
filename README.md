# Shocktube_SPH

Simulating the classic sod shocktube using SPH technique. 

The Sod shock tube problem, named after Gary A. Sod, is a common test for the accuracy of computational fluid codes, like Riemann solvers, and was heavily investigated by Sod in 1978. (Source : Wikipedia)

The initial conditions for the problem describe two states (left and right) of a quiescent gas separated by an imaginary diaphragm. The states are given as (ρl , pl , ul ) = (1.0, 1.0, 0.0) and (ρr,pr,ur) = (0.125,0.1,0) on the left and right hand sides of the diaphragm which is placed at x = 0.

Results using Summation Density for density udpation : 

![alt text](https://github.com/deeptavker/Shocktube_SPH/blob/master/figure_sd_output.png)

Results using Continuity Equation for density updation :

![alt text](https://github.com/deeptavker/Shocktube_SPH/blob/master/figure_ce_output.png)
