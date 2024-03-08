# ThirdSound

Solves for standing wave solutions of third sound (surface waves) on thin-film helium of liquid depth $h$, in a 1D rectangular box. Hard-wall boundary conditions are assumed on the bottom of the container and the side walls. 

Uses the numerical method developed by A. J. Roberts [IMA Journal of Applied Mathematics, Volume 31, Issue 1, July 1983, Pages 13–35.] The standing waves are found by using a method very similar to that  outlined by Smith and Roberts [Physics of Fluids 11, 1051–1064 (1999)], however, instead of arclength continuation, the nonlinear system is solved using Powell's dogeleg method [as used in e.g., Mercer and Roberts, Wave Motion Volume 19, Issue 3, May 1994, Pages 233-244 ]. Upgrading to arclength continuation would be better for dealing with resonances, as seen in Smith & Roberts. 


The equations for the surface $\eta(x,t)$ and potential $\phi(x,y,t)$ are 

$$ \partial_t \eta + \phi_x \eta_x = \phi_y$$

$$ \partial_t \phi - \alpha \left[\frac{1}{(h + \eta)^3} - \frac{1}{h^3} \right] + \frac{1}{2}( \phi_x^2 + \phi_y^2) =0  $$.

The code parameterizes $\eta(x,t) = \{X_j(t), Y_j(t)\}$, and solves for a profile which repeats itself after one period (the period is also unknown and is found iteratively). 

The code is run by calling the function `doglegSolve(N,h,A,waveType,delta_A,Amax)`. Here `N` is the number of grid points, `h` is fluid depth, `A` is the wave acceleration at $t=0$, (which parameterizes the solutions). Wavetype may be either `'VdW'` or `'Gravity'`. Gravity replaces the restoring force with the ordinary linear restoring force for water waves. `deltaA` specifices the value to increment the acceleration, and `Amax` is the largest wave acceleration sought. The acceleration is related to the height of the wave but has been shown to be more suitable for seeking numerical solutions under some circumstances. 

Note that the code assumes all inputs are strings (this is needed to run on getafix cluster). For example you can run:

`doglegSolve('16','0.6','0.001','VdW','0.001','0.01')`

which will find solutions for $h=0.6$, using $N=16$ points, for accelerations ranging from 0.001 to 0.01 in increments of 0.1.  


An example solution for the surface profile is shown below, for $h=0.2$. This partciular example is a cnoidal-like wave. 

<p align="center">
 <img width="611" alt="Screen Shot 2024-03-07 at 2 44 27 pm" src="https://github.com/Matt-Reeves/ThirdSound/assets/65841999/9a79b4a0-afb2-4a98-9194-c03124d58ca3">  
</p>

# Potential Avenues for future work

As shown in `HeliumFimStandingWaves` I calculated frequency corrections for Stokes, Cnoidal, and Solitary wave solutions. Surface tension was neglected. 

Suggestions for future work (in order of increasing difficulty)

- I mostly only looked at the profiles and frequencies. Many key optomechanical metrics could still be calculated using only the existing code (effective mass etc., see, e.g. Baker 2016).
- Upgrade to use pseudoarclength continuation as per Smith and Roberts. This would better deal with resonances (see e.g. their paper), and would allow probing of extreme amplitude waves for thin films (i.e., wave amplitudes $\gtrsim 50$% of the film depths).
- Upgrade to include surface tension. Compare with the perturbation theory results of Concus, 1960. His calculation extends on the work of T&K, which I followed to calculate the Duffing frequency for Stokes waves. (Relevant Mathematica notebook is also in this repo). Compare also with soliton predictions of Nakajimia (1980).  Solving for surface tension numerically, probably need to change RK4 timestepper to a semi-implicit method. 
- Upgrade to investiage other geometry effects, e.g., non-flat bottom of the basin. 



 
