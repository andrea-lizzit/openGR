# OpenGR - Open Gravitational Renderer

OpenGR is a small software that simulates the propagation of photons from a camera to an infinitely distant point in a spacetime distorted by the presence of a wormhole. Given the parameters that specify the camera and the wormhole, OpenGR produces a map that associates the final direction of the photon in the camera frame to the asymptotic location on a set of celestial spheres from which the photon comes from.

OpenGR propagates a photon backward in time from the final location on the camera to the initial location in a curved spacetime. The specific spacetime geometry determines the path the photons take, the geodesics. Here only one type of wormhole is considered, the DNEG wormhole, which produces relatively simple and computationally fast expressions. For this type of wormhole the geodesics can be calculated by numerically solving a set of differential equations subject to initial conditions that represent the initial position and direction of the photon. This project implements ideas and equations described in [[1]](#1), which were used in the production of the movie _Interstellar_ by the VFX team [DNGR](https://www.dneg.com/)

The calculation of the map is performed in Julia. It is possible (and highly recommended) to enable multithreading by setting the environment variable `JULIA_NUM_THREADS=auto`.

An environment map to be used in any standard CG software can be produced from the coordinate map with the `compose_image.py` script. Run `compose_image.py --help` for options.


## References
<a id="1">[1]</a>
James, Oliver, Eugénie von Tunzelmann, Paul Franklin, and Kip S. Thorne. ‘Visualizing Interstellar’s Wormhole’. American Journal of Physics 83, no. 6 (1 June 2015): 486–99. https://doi.org/10.1119/1.4916949.
