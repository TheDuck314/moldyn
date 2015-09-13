# moldyn

Go play with Daniel Schroeder's excellent [molecular dynamics simulation](http://physics.weber.edu/schroeder/md/).
It simulates the behavior of bunch of atoms with an attractive force between them. Try varying the temperature to produce
a solid, a liquid, or a gas. Are you intrigued?

I was, so I wrote my own version of the simulation and played with it a lot. It turns out that this simple 2D system
shows a lot of interesting behaviors that are important in thermodynamics and statistical mechanics. You can learn about
phase transitions, latent heats, triple points, critical points, sublimation of solids into gases, heat capacities,
and probably a bunch of other stuff.

# The phase diagram
The most interesting thing to look at is the phase diagram. If you have played with Schroeder's simulation enough
you'll know that you can get the atoms to behave like a solid, a liquid, or a gas depending on the circumstances. Specifically,
the phase depends on the temperature and the pressure. The phase diagram looks very much like Wikipedia's example of a typical
phase diagram:

IMAGE

Schroeder's simulation basically lets you control two things: the total energy of the atoms and the volume
they occupy. To investigate the phase diagram it's more useful to be able to control the temperature and pressure,
which are related to the energy and volume but aren't the same thing. So my simulation code is written to allow
simulations at a given temperature and pressure. By running a bunch of these simulations we get the following phase
diagram:

IMAGE

General discussion of phase transitions. Sublimation. Triple point. Critical point.

# How to find a phase transition
Discontinuities in energy and volume. Latent heat.
Heat capacity. Thermodynamic limit. Suddenness of transition, with pictures.
Why Lvaporization > Lfusion.

# The triple point
Tvap and Tmelt come together. Discuss why liquids can only exist above a certain T and P.
Below triple point, only phase transition is sublimation.

# The critical point
Latent heat goes to zero. Supercritical fluid. No liquid-gas transition above Pcritical,
only transition is melting.

# More things
Compressibility. Thermal expansion coefficient. Speed of sound. 

# Technical details
Implementing constant-pressure. Implementing constant-volume. Interatomic potential. 

ANIMATED GIFS
show solid, liquid, gas, evaporation, freezing, sublimation
