Introduction :


Numerical methods play a pivotal role in simulating dynamic processes in astrophysics, where analytical
solutions are often elusive or computationally infeasible. These methods provide a means to simulate and
analyze the dynamic behavior of celestial bodies, gravitational interactions, and evolving cosmic structures.
One widely employed numerical technique is the Runge-Kutta method, a family of iterative algorithms
specifically designed for solving ordinary differential equations (ODEs). In astrophysics, ODEs frequently
arise when modeling celestial systems, such as the motion of planets, the behavior of stars, or the evolution
of galaxies. Runge-Kutta methods excel in providing accurate and stable solutions to ODEs by iteratively
refining approximations of the solution over discrete time steps.

The initial phase of this project, named ’Galactic Disk Heating’, involved integrating a circular orbit
within the Miyamoto-Nagai potential, leveraging a fourth-order Runge-Kutta integrator to establish crucial
physical scales and integration parameters. Building upon this foundation, we develop a class capable of
extending these orbital calculations to numerous stars, enabling a collective examination of their trajectories.
To ensure a realistic density profile for the final galaxy, initial positions are strategically selected from a
uniform distribution. The Schwarzschild Distribution Function is introduced, including random velocity
components in the radial (R) and vertical (z) directions to impart distinctive characteristics to the stellar
motion while closely aligning with circular paths.
A key addition to our simulation was the incorporation of a point-mass perturber, colloquially referred
to as a black hole in our case, with a mass ($M_{BH}$) of $10^{6}$
solar masses. The ensuing stages of the project
involve the systematic simulation of diverse perturbation types and a rigorous analysis of their impact on
the galactic disk.
Through this project, we manage to untangle the intertwined processes governing galactic disk heating,
expanding on the mechanisms responsible for the observed evolution of its dynamics.

Conclusions:


In all scenarios, we have observed that the perturbation creates a significant density wave propagating
through the galaxy. In the first case, featuring an isotropic perturbation, ring densities are formed. In the
three other cases, where the perturbation is anisotropic, we observe spiral density waves, which manifest as
spiral arms. All perturbations lead to an increase in radial velocity dispersion $\sigma_{v_{R}}$, in line with the trends observed in Fig. 20.
After the black hole’s perturbation, kinetic energy is gained on average, with a more substantial gain for
stars with low epicyclic amplitude. The consistent increase observed in both radial velocity dispersion ($\sigma_{v_{R}}$)
and kinetic energy across all scenarios involving perturbations substantiates the aptness of our project’s title,
’Galactic Disc Heating.’
