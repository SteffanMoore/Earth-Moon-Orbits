# Earth-Moon-Orbits
Simulation of an orbiting object around the Earth and Moon using the fourth order Runge-Kutte method to simulate orbits of a spacecraft around the Earth and Moon.
This is an numerical approach which solves the differential equations for accelaration iteratively. It uses the known x and y position (in relation to the Earth's centre) to calculate the gravitational forces acting on the spacecraft before combining these with the cartesian velocities in the so-called "K vectors" to estimate accelarations at different points within one iteration. These are then combined and weighted to find new position and velocity values.

The initial approach uses a 2D system with the Earth and Moon being treated as fixed masses; using 2D is acceptable in this scenario as all orbits in this system can be interpreted as along a plane.
