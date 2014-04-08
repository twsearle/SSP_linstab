# Viscoelastic self sustaining process paper plan #

## Introduction ##

* Why have I done this?
* Which other work inspires this?
    - Waleffe 
    - People talking about lift up
    - people talking about What sustains turbulence
* Summary of the point
    - get streaks in purely elastic regime
    - streaks are unstable 
    - Wavy, but also not wavy instabilities
    - Free slip boundary conditions turn out to be important

## Methods ##

* Chebyshev - Fourier decomposition and then linear stability
* Newton Rhaphson method on the solution forced with rolls
* Then do linear stability analysis on these solutions with a disturbance of a specific kx fourier mode.
    - Give the free slip boundary conditions you use
* Repeat this for all kx to build up a picture of the dispersion relation

## Results ##

### Streaky profile ###

* See the lift up mechanism when the amplitude of the rolls is sufficiently large.
* Just like Waleffe's result at low Wi 
* Rather then large inflection in velocity you see regions of high Normal stress difference.
    - Corresponds to shearing of polymers?
* Solutions give a streaky streamwise velocity profile
    - Streaks in the streamwise velocity
    - Also streaks in the 1st Normal Stress difference, with Txx being  the largest contribution to the polymer stress tensor
    - Little change on altering amplitude, etc

### dispersion relations ###

* As Reynold's number decreases, the instability noticed by Waleffe goes away.
* At intemediate Reynold's numbers and Weissenberg numbers (elasticity ~ 0.1 ?) the instability disappears.
* At low Re and high Wi (elasticity ~ 1000?), we see a purely elastic instability arise at very low kx. 
* This instability is increased by reducing Reynold's number or increasing the Wiessenberg number.
* The instability plateaus at a Wiessenberg number of ~ 18 we think.
* As the Reynold's number is reduced
    - instability moves to lower kx and reduces in strength
    - Newtonian instability disappears at about Re = 100 (more accurate?)
    - Viscoelastic instability appears at very low kx ~ 0.02, and very low Reynold's number, Re = 0.01
    - With further decreases in Reynold's number, the viscoelastic instability becomes stronger whilst maintaining its width.

### eigenmodes ###

* Although the instability is large at the walls, most of the gradients in v0,v1,w1 take place away from the walls. 
* It is the gradients that are responsible for the variations in the normal stress and so can reinforce the rolls, as in the original self sustaining process of Waleffe.
* Examine eigenmodes that correspond to these eigenvalues
    - Eigenmodes show gradients in the velocity which are concentrated in the centre of the channel.
    - These are responsible for the large stresses here.
    - This corresponds to the region where the Kelvin-Helmholtz instability is likely taking place
* Need to actually check the nonlinear forcing of rolls by the instability

## No slip case ##

* Without slip at the boundaries the instability appears infinitely amplified
* [Show var slip graph]
* A possible explanation for this is that the instability moves towards the walls and can no longer be resolved due to large gradients in the velocity.

## Discussion ##

* A purely elastic instability has been uncovered in a situation analogous to that of the self sustaining process, a fundemental component coherent structure in Newtonian turbulence in wall bounded shear flow.
* The instability is very different to that of the Newtonian case, it remains to be shown that it can reinforce the streamwise rolls and complete the self sustaining process.
* It occurs at zero kx, meaning it does not correspond to a kelvin-helmholtz like instability of the streaks in the streamwise direction?
* The potential importance of free slip boundary conditions, known for ages for polymeric fluids, has been explored. Free slip and no slip seem to behave differently. 

## Conclusions ##


## STUFF TO CHECK ##

* Am I following the template?
    - bold grad operator?
    - exponential right?
    - Weissenberg number correct?

* How do I eliminate pressure from the base profile?
    - don't understand how the boundary conditions and the pressure relate

* Have I written the decomposition of the disturbance right? 
    - signs on time?
    
* Do I want a plot of Wi = N1/Txy ? A sort of spatially dependent Weissenberg number
    - Does that even make sense? - Read Alex

* I need to actually check the non-linear forcing of the rolls by the eigenmode of the instability.

* Fix up eigenmode plot so that there is no zeroth component of N1. Also think about add txy plot as well.
    - consider not including this at all, instead have only the version with Cauchy boundary conditions.

* Varying the amplitude? Perhaps lower amplitude I won't have kx = 0 solution?
    - What can I do about kx=0 solution?

* Express var slip parameter as proper slip length of the polymeric fluid?
    - at this point I don't understand the proper way of talking about slip in fluid dynamics
    - could read about it and find out?

* Variable slip dispersion relations. Are they complete?
    - Should I extend to higher kx?
    - Lower $\alpha$? Or can I not resolve that?
    - some sort of better tracking of one eigenvalue?

* Is it the First normal stress difference or the gradients in the first normal stress difference that matter?
    - Seems from the heat maps that it is the gradients where we expect instability to happen
    - trying to convince people that it is only the magnitude that matters in the eigenmode.
    - Contradictory!

* AGAINST THE RULES: free free slip might be useful if I change the amplitude of the rolls
    - alternatively, check the most unstable eigenfunctions, see if they behave as you would expect

* Do I want to include my 3D vorticity plots?
    - vorticity? O2 or Ox2? - that is the waviness that is important?

* Can I plot streamwise independent instability on top of the base profile somehow?
    - COOL IDEA

* Should I include plots of Cxy anywhere? 
    - eigenvectors
    - base profile

