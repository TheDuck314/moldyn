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
you'll know that you can get the atoms to behave like a solid, a liquid, or a gas depending on the circumstances.

The solid phase is a crystal with a hexagonal lattice structure:

IMAGE

If you heat the solid up enough it melts into a liquid. The liquid is almost as dense as the solid but the atoms are not locked into a crystal structure:

IMAGE

If you heat the liquid up enough it boils, producing a gas, in which the atoms mostly fly around on their own, occasionally bumping into each other:

IMAGE

Which phase you get depends on the temperature and pressure. I ran a bunch of simulations to determine the phase diagram. Compare the phase diagram for water on the right.

IMAGE, IMAGE

You can see that the phase diagram is quite similar to the phase diagram for a realistic substanc like water. One difference is that for water, the melting temperature decreases at high pressure, whereas here the melting temperature increases at high pressure. The difference is because water has the weird property of expanding when it freezes, whereas the simulated substance is more normal in shrinking when it freezes.

At the *triple point* (link), all three phases come together. The triple point is around (T=0.395, P=0.009). At temperatures or pressures below these values, the liquid phase can't survive: it either boils into gas or freezes into solid. Below the triple point, there is a line separating solid and gas phases. So if you have a solid at low enough pressure and heat it up, it *sublimates* (link) directly into a gas without passing through the liquid phase.

The line separating the liquid and gas phases corresponds to vaporization (boiling). This phase transition line has the interesting property that it ends at a certain point, the *critical point* (link). The critical point is around (T=0.50, P=0.039) (check). At temperatures and pressures above the critical point, there is actually no sharp transition between the liquid and gas phase, so in this regime we just call the substance a "supercritical fluid."

# How to find a phase transition

How can we pin down exactly where to draw the lines on the phase diagram? Is there really a sharp distinction between the liquid and solid phases, for instance? Or does the liquid just gradually solidify as it gets colder?

In fact the boundaries are sharp. At phase transitions, various properties of the substance change discontinuously. Two properties that are easy to look at in the simulation are the density and the energy per atom (really, the enthalpy (link) per atom).

IMAGE

Here's an example. I've plotted the energy per atom as a function of temperature at a pressure of ????. At low temperatures the substance is a solid. Then at T=???? the energy jumps up. This is the solid melting into a liquid. At T=???? the energy jumps again by a much larger amount. This is the liquid vaporizing into a gas. So from this graph we can read off the melting point and boiling point at P=????. 

The size of the discontinuous jump in energy at a phase transition is called the *latent heat* (link) of the phase transition. It's the amount of energy you need to pump into the system to drive it from the colder phase to the hotter phase. For the solid-to-liquid transition, the latent heat is basically the energy needed to break each atom out of its spot in the hexagonal crystal lattice.For the liquid-to-gas transition, the latent heat is the energy needed to pull each atom out of the liquid into free space. 

The latent heat is kind of large because the atoms attract one another and it takes a lot of energy to pull them apart. The latent heat of vaporization is much larger than the latent heat of fusion (melting). This is very typical of real substances. For example for water the heat of vaporization is almost 7 times the heat of fusion. 

Another property we can read off the above graph is the *heat capacity* (link) of the solid, liquid, and gas phases. The heat capacity is the amount of energy needed to increase the temperature--that is, the derivative of the energy with respect to the temperature. So the heat capacity is the slope of the graph (in the regions away from the phase transitions). In this regime the heat capacity of the solid is about ????, the heat capacity of the liquid is about ????, and the heat capacity of the gas is about ????. (There are really two heat capacities, the heat capacity at constant pressure and the heat capacity at constant volume. Here we are measuring the heat capacity at constant pressure).

Above we plotted the energy per atom. We can also plot the volume per particle. (Note that since are are in 2D instead of 3D, by "volume" we really mean "area"). 

IMAGE

Again the graph consists of smooth regions separated by discontinuous jumps at the phase transitions. We can see that the volume increases by about ????% when the solid melts into a liquid, and the volume increases by about ????% when the liquid boils into a gas. (This latter number is strongly dependend on the pressure).

In between the phase transitions the volume varies continuously with temperature. The slope of the volume versus temperature graph gives the *coefficient of thermal expansion* (link). As for most substances, the volume increases with temperature. The coefficient of thermal expansion of the solid is about ???? and the coefficient of thermal expansion of the liquid is about ????

## Supercooling and superheating

The discontinuities in the above graphs provide a nice way of locating the phase transitions. However there are some subtleties. It turns out that the phase transition occurs at a different temperature depending on if you start at a low temperature and make the substance hotter, or start at a high temperature and make the substance colder. Here's an example:

IMAGE

Basically, the substance likes to stay in its current phase and it doesn't immediately transition to the other phase when you pass the transition temperature. If you start with a liquid and heat it up, you can heat it up by a small amount past the boiling point without it actually boiling. This is called *superheating* (link). If you wait long enough, or if you increase the temperature some more, the liquid does boil. Similarly if you start from a gas, you can *supercool* (link) it to a bit below the boiling point without it condensing. 

One way to understand this is that it's kind of hard to get the phase transition started, but once it does get started it proceeds steadily. The gas condenses into a liquid by forming little initial droplets of liquid and then accreting onto the droplets. But it can be hard for the initial droplets to form. Supercooled gas condenses as soon as a large enough liquid droplet forms. Droplets form faster at lower temperatures, so the supercooled gas tends to stay gas until you cool it a bit below the boiling point. Eventually a droplet forms and the gas rapidly condenses onto it.

# The critical point

What is going on at the critical point, where the liquid-gas phase transition line ends? Basically, as you move along the liquid-gas transition line, the difference between the liquid and gas phases becomes smaller and smaller until finally at the critical point the difference vanishes. 

One way to look at this is to plot the latent heat of vaporization as a function of the temperature:

IMAGE

At the critical temperature, about T=0.50, the latent heat of vaporization becomes zero, after which there's no sharp distinction between liquid and gas. 

A similar thing happens with the volume per particle. Here's the difference in volume-per-particle of the liquid and gas phases, as a function of temperature:

IMAGE

Again the difference between the phases goes to zero at the critical temperature.

Yet another way of looking at this is to plot "isotherms," lines of constant temperature in the pressure-volume plane:

IMAGE

The discontinuities in the isotherms for T<0.50 represent the discontinuous volume change of the liquid-gas phase transition. For T>0.50 there is no discontinuity in the isotherm.

UNIVERSALITY?

# Properties of the gas phase

Equation of state. Van der Waals fit?

# Compressibility and the speed of sound

# Technical details
Implementing constant-pressure. Implementing constant-temperature. Interatomic potential. 

ANIMATED GIFS
show solid, liquid, gas, evaporation, freezing, sublimation
