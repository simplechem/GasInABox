# What Is This Program For?

I originally wrote this program to provide a molecular scale visual
illustration of the difference an ideal gas and an interacting gas. As
the code grew, I added features for:

- limited keyboard input to adjust parameters.
- live updates to display running averages and distributions.
- input file including updating conditions.
- limited output.
- and of course, optimizations
- choice between Lennard-Jones or Morse potentials.

So, does it work? I think I achieved my goal of a visual
illustration. Depending on how many particles you put in the animation
and what you want to demonstrate, it's usefulness may vary.

# Why Did I Create This.

I teach chemistry at the University of Manitoba in Canada. During the
Winter 2023 term I had the pleasure of teaching our 3000 level
thermodynamics course. After teaching online during COVID, I found it
frustrating to explain the differences between ideal gases and
interacting gases just by scribbling on the white board and waving my
hands around. I wanted something that gave a different kind of visual
illustration of the effects I was describing. I tried Manim and found
it wasn't a good tool for this job. I'd been tinkering with Rust for a
bit, so I started my journey to write this program. I could have used
C with either pthreads or openMPI and some graphics library, but
wanted to make this a learning exercise for myself as well, so I took
on the challenge of doing this in Rust with Bevy as the graphical
library (crate).

# What Can You Demonstrate with This?

It depends a bit on the situation. If you want to just look at what
happens as you change conditions, then I'd recommend just starting
the program with no input file. It defaults to 4096 particles with
infinite iterations. You can use the keyboard input to adjust
parameters live and see how the system responds to changes in things
like the temperature or the size of the box. With this number of
particles, the system should be pretty responsive on most modern
systems. You can cool it down and see things condense into a liquid or
even crystallize into a solid. You can turn the potentials on and off
and see how the ideal gas behaves differently. Just goof around with
it.

If you adjust the box size in a live setting, be careful to let the
system equilibrate from time to time, since going too fast (especially
when using Lennard Jones potentials) can cause infinities and
eventually crash the sim.

Alternatively, you can setup an input file and program a series of
conditions. Maybe trace out isotherms, for example. Runs like that
tend to take a while if you want good data, so set it running and come
back in a week. To run a specific input file, you need to run the
program from the command prompt and pass the name of the input file
with the -i or --input flag. For example:

gasinabox -i input.gbi

(I have chosen to set .gbi, Gas Box Input, as an extension, but it's
not necessary).

Some interesting aspects:
- Depending on temperature and box size, the interacting gas can show
  an equilibrium between liquid and vapor.
- Depending on temperature and box size, the interacting gas can
  crystallize into a solid.
- If either condensed phase has the potential turned off, the
  particles will scatter and behave like an ideal gas.
- Thermal aspects are based on Monte-Carlo methods which can take a
  long time to reach equilibrium. The "thermalize" key can help with
  this.
- Expanding or shrinking the box is practically adiabatic (it's fast)
  unless you use "autotherm" in an input file. This means that you can
  see some interesting heat flow after each volume change while the
  systems reaches equilibrium.
- In addition to numerical data, there are also two graphs that are
  displayed. These are the pair correlation function and the speed
  distribution function. Axis labels will be added in a future update.
  
Finally, when runs finish, they dump their state to a file. You can
use the -r or --restart flag and one of the state files to restart a
completed case from a run. For example:

gasinabox -r results.0025.state

# What Is Under the Hood?

The program runs mostly as two parts, a graphical portion based on the
Bevy engine, and a physics portion. The graphical side periodically
requests information from the physics side to update the positions of
the particles and other data displayed on the screen. The physics side
runs a loop that goes through the following main steps:
- The separation between each pair of particles is evaluated and used
  to determine contributions to the potential energy associated with
  each particle, and the acceleration of each particle.
- Currently the potential used can either be a Lennard-Jones potential
  based on literature data (see sample input file for references), or
  a Morse potential calculated from Lennard-Jones data to have the
  same minimum as the Lennard-Jones potential. The Lennard-Jones
  potential is a slightly faster calculation, but the Morse potential
  has the benefit of not having infinities at small distances so Morse
  was chosen as the default.
- Once the accelerations are all known, the velocities are updated.
- Once the velocities are updated, the positions are updated.
- Once the positions are updated, the system checks for collisions
  with the wall.
- Collisions are addressed by moving the particle back into the box,
  reversing the component of velocity normal to the wall, and choosing
  a Boltzmann weighted random magnitude for the normal velocity.
- Collisions with the wall provide the mechanism for the Monte-Carlo
  aspects of the physics.
- After all of this is done, we check if any macroscopic parameters
  need to be averaged or if the graphics side has asked for an update,
  or if there are any requests to change parameters like temperature
  or box size.
- Then we start all over again.
- Most of the effort has gone into methods to decrease the calculation
  time for the accelerations since this, in principle requires us to
  compare N(N-1)/2 pairs of particles. Several strategies have been
  used to speed this up, but it is the bottleneck. 4096 particles is
  pretty smooth on most modern systems (as of 2025), but 32768 starts
  to push most systems to the limit in the current version of the
  code. The next step of 262144 particles is beyond the capacity of
  most common computers to work at right now.

# Are There Important Differences Between L-J and Morse Results

Yes, there are. Because of the exponential decay in the Morse
potential, it is possible to set a much smaller "closescale" (see
below when the parameters are described). This speeds up the
calculations a great deal, making the Morse potential faster to work
with. However, the long range, weakly attractive tail in the L-J
potential does lead to significant differences if you try to generate
a phase diagram with this program. The Morse potential does not lead
to the correct critical point behavior in anything like the range of
conditions seen for real neon. However, the L-J potential does give
approximately the correct behavior, even if the exact temperature and
molar volume don't match that of real neon. So, if you want something
visually smooth, go with Morse. If you want something that correctly
gives critical behavior, go with L-J.

# Why Is the Default Neon-20?

Originally I started with He-4. In some ways, you would expect this to
be nice since it is the simplest noble gas from a typical chemistry
perspective. However, if you just use the Lennard-Jones inter-atomic
potential for helium, you get results for the critical point that are
way off the experimental values. The reason for this is that LJ does
not account for interactions between the nuclei that will contribute
to the interactions. On the other hand, if you go to something larger
like xenon, you will have to include three-body and higher
interactions for an accurate result. I'm not sure where the need for
higher order interactions becomes significant, so I just jumped up to
neon in the hopes that it would be good enough to be physically
realistic if I tried to generate a phase diagram. If you don't care
about accuracy, then you can do any monatomic gas.

# Why Did I Use the Rust/Bevy Combination

Initially I tried to do this in Manim. Manim is a nice tool and is
absolutely amazing for lots of things, but not this. The underlying
physics here require lots of nested loops, and this just gets slow in
python. It took my computer about 20 min. to generate a 2
min. animation of 64 atoms bouncing around. This code can run 32768
particles live on the same laptop and still give reasonable
performance. So I abandoned Manim. In the past I would have tried this
with C, OpenGL, and pthreads. But, I've been tinkering with Rust for a
little while and wanted to use this as a learning project. Overall,
I've liked the experience of coding in Rust. The built in concurrancy
made some things nice to work with, and but safety and compiler
feedback really helped me find and squash bugs a lot faster than I
typically get with C. Also, I've tried C++ and Ada for projects like
this in the past, and I found Rust easier to learn than either.

# Bevy

This code uses Bevy 15.3. It probably doesn't use it well, but it uses
it. Bevy is a great project that has been relatively easy to learn for
me, and made developing the rest of this much simpler since stupid
physics mistakes in my code became much easier to spot. I am also
amazed at how a project like this can give such effective animation of
so many objects. Great job Bevy! I am using Bevy under the MIT license.

# How Does It Work?

There are two main parts to the original plan. In the source code you
will find two modules called gasbox and eyecandy. Gasbox takes care of
the physics of the simulation while eyecandy takes care of drawing
things to the screen. Periodically (ideally 60 times a second)
eyecandy asks gasbox for updated coordinates and gasbox send them back
to eyecandy along with data that eyecandy displays on the
screen. Apart from some initialization, almost all of the Bevy code is
in the eyecandy module. Almost all of the physics code is in the
gasbox module. Crossbeam channels are used to send information back
and forth between the two. Eventually I added extra modules to keep
things like fileio and averaging straight.

# What Happens in the Gasbox Module?

The gasbox module has one main loop that it goes through over and over
again. This loop can be broken into several stages:

## Stage 1: Update Accelerations and Potentials.

This is where the main work happens. In principle we need to go
through every possible pair of particles (N*(N-1)/2 pairs total) and
work out the separation between that pair. Once the separation is
known we use that to calculate a contribution to the potential energy
of the particles in that pair and a contribution to the total force
that each particle experiences.

There are strategies to eliminate pairs that are too far away, which
really helps speed up the calculation.

## Stage 2: Update the Velocities.

Once we know all the forces experienced by each particle, we update
the velocity using v(new) = v(old) + dt*acc. Since this only scales like
O(N), there hasn't been much effort to fine tune this yet.

## Stage 3: Update the Positions.

Once we have an update on how fast everything is moving, we update the
positions based r(new) = r(old) + dt*v. Again, this only scales as
O(N), so hasn't been optimized much.

## Stage 4: Bounce Check.

Check if any particles have moved outside the box. If they have, move
them back in, flip the normal component of their velocity, and use a
random number generator with a boltzmann weight to pick a new
magnitude of the perpendicular velocity.

This gives us two things that we can do. Firstly, it gives us a way to
reach thermal equilibrium. Secondly, we can add up the impulses from
all of the collisions with the wall and build up a picture of the
pressure.

## Stage 5: Update Other Modules

Check to see if we have any requests for coordinates or data for
avergaing, etc. that has to be sent down channels to other parts of
the code. Also check for keyboard input.

And that's about it. I might write extra explanations at some point,
but for now this will give you and idea of what's under the hood.

# What's In the Log Files?

The log files are a periodic dump of the following information:

- column 1: current iteration
- column 2: number of values used in each average
- columns 3, 4, 5: current kinetic energy, mean kinetic energy and
  mean kinetic energy squared.
- columns 6, 7, 8: current, mean and mean of square for potential
  energy.
- columns 9, 10, 11: current, mean and mean of square for pressure.
- columns 12, 13, 14: current, mean and mean of square for heat flow.
- columns 15, 16, 17: current, mean and mean of square for iterations
  per second.
  
