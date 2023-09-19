# PICFLIP

## A 2d FLIP/PIC fluid based on Animating Sand as a Fluid (2005), Zhu and Bridson. 
* Programmed methods for transferring velocities between a MAC grid and particles, solving the incompressibility constraint, and calculating the new velocities.  
* The simulation results in realistic water looking fluids but the ratio can be adjusted to mimic liquids with different viscosities.
### A full report made in conjunction can be seen here: [PICFLIP_Paper.pdf](PICFLIP_Paper.pdf)

Here's a video with the results:  
[![Fluids](https://img.youtube.com/vi/5IXsQ7yA04g/maxresdefault.jpg)](https://youtu.be/5IXsQ7yA04g "FLIP Fluids")  
https://youtu.be/5IXsQ7yA04g

The basic controls are available within the program window. 
To change the test particle system, navigate to GridSolver.h and change the value of the pSystem variable
You can also change the gravity applied as well as the number of particles (modify particlesX and particlesY)

Lastly, to change the cell size, use the slider then press the "r" key to reset the simulation.

Be use to have all the 3rd party tools from assignment 3 before using cmake. Follow the set up as described before.




