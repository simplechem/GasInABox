# How Do You Use This?

## Just Run It

There should be binaries available for 64-bit Linux or for Windows
where you downloaded this from. Just run the executable. It defaults
to 4096 atoms of Ne-20 bouncing around in a box. Various keys can
change conditions:

- Q will quit the program.
- I will toggle the potentials on and off.
- A will "thermalize" the system (see below).
- E and D can be used to increase or decrease all velocities.
- T and G can be used to increase or decrease temperature.
- U and J can be used to increase or decrease the volume.
- R and F can be used to adjust the time-step in the physics.

"Thermalizing" looks at the current average kinetic energy and scales
all velocities to make this average be 1.5kT. If you are adjusting
things live, this can be useful to reach thermal equilibrium more
quickly.

## Run With an Input File

From the command line, you can supply an argument with the name of an
input file to use. The input file has the parameters described
below. There are three kinds of line in the input file. Comment lines
start with a semicolon (;), and anything in a line after a semicolon
will be ignored. There are "key = value" lines that specify different
parameters for a run (summarized below) and there are variable lines
that describe how to do a sequence of runs while varying one of the
parameters. The variable lines will be explained after the parameters
are all explained.

### basename

A string that describes the base to use when generating output files
to store data.

### element

A string that describes the element used. Currently this is only for
display purposes.

### sigma

A float that represents the Lennard-Jones distance parameter
(sigma). The value is assumed to be measured in metres (m).

### epsilon

A float that represents the Lennard-Jones energy parameter
(epsilon). The value is assumed to be measured in joules (J).

### mass

A float that represents the mass of a single particle, measured in
kilograms (kg).

### tcrit, pcrit, dcrit and vcrit

A group of floats that represent critical values in K, Pa, mol/m^3 or
m^3/mol.  Note that if both dcrit and vcrit are given, the last one
will override the first, so only one is needed. these are used only
for display.

### molardensity and molarvolume

Floats that represent the molar density or molar volume (mol/m^3 or
m^3/mol) to operate at. Only one or the other is needed. If both are
given, the last one will override the earlier one.

### temperature

Float for the operating temperature in kelvin.

### closescale

Float for the cutoff to stop the calculation for distant atoms as a
multiple of sigma. 3.0 is good for morse potential, between 6.0 and
10.0 is good for LJ without loss of accuracy in the calculation.

### velscale

A scaling factor for the initial randomization of the velocities. 1.6
is pretty good.

### molfact

An integer. The number of particles is set to 8^molfact. A value of 3
gives 512 atoms and a very fast moving simulation unless dtsteps is
made very large. 4 gives 4096 particles and good visuals. 5 gives
32768 particles with rapidly converging averages, but sometimes
cluttered visuals. 6 would give 262144 particles and will be uselessly
slow on most systems as off 2025.

### itermax

An integer. The maximum number of iterations to use for the run. Set
to 0 for infinite runs.

### anneal

An integer. How many iterations to run before resetting the iteration
counter and starting the main run. You might want to do this to allow
time for transient conditions to run out at the start of a run or
after a change in conditions.

### sample

An integer. How many iterations between sampling quantities that we
keep running averages for.

### ntrack

An integer. How many values are tracked for the running averages.

### ndist

An integer. How many bins do we used for distributions like the speed
distribution or the pair correleation.

### dtsteps

An integer. How many steps should a particle moving with speed that
gives kinetic energy of 0.5kT to cross a distance of sigma.

### evalmethod

An integer 0, 1, or 2. 0 and 1 are two methods for calculating the
Lennard-Jones potential based on passing in r^2 (for 0) or r^6 (for
1). Method 2 uses a Morse potential tuned to have the same minimum as
the Lennard-Jones potential.

### accmethod

An integer 0, 1, 2, 3, or 4. 0 is the naive method that calculates
everything within the limits of closescale. 1 applies a grid and only
even looks at pairs within the same grid space, or in neighboring grid
spaces. 2, 3 and 4 Use the grid and assigns different grid
combinations to different workers in concurently (see below for
differences). All five options still make use of closescale, but 1 and
2 use the grid to avoid even checking many of the cases that would
fail closescale.

The difference between 2, 3 and 4 is the strategy of how the jobs are
portioned out. The grid is a cube of "sectors". Method 2 assigns an
entire layer to each worker, 3 assigns entire lines to each worker and
4 assigns individual sectors to each worker. Method 2 is likely to be
faster when things are evenly distributed since there is less time
spent passing information back and forth between the main job and the
workers. Method 4 is likely to be faster when things are
condensed. Method 3 is a compromise. Specific cases may give varying
results. I may introduce a method 5 that tries all three during
annealing and picks the best, but for now ... user has to do the trial
and error.

### voltherm

An integer 0, 1, or 2. 0 does nothing. 1 runs the thermalizing routine
any time we change the size of the box, and any time we sample
averages during the annealing stages. 2 takes this a step further and
restarts annealing if the mean kinetic energy is more than 0.1kT away
from 1.5kT when we get half-way through the annealing. This way we
have good averages throughout the entire run after annealing
completes.

### logaves

An integer 0, or 1. 0 for off, 1 for on. Do we save the averages to a
file at regular intervals.

### periodic

An integer 0, 1, or 2. 0 The box has hard walls and exists in a
vacuum. 1 The box has hard walls and is surrounded by identical boxes
so that interactions between atoms are felt across the boundary. This
has the effect of a particle at the right wall and a particle at the
left wall being able to interact with each other. 2 is similar to 1,
but allows a particle at the right wall to teleport to the left wall
rather than bouncing off the wall.

### nworker

An integer. How many workers to use when using accmethod 2. This one
requires tinkering on most systems. To really see the effect, play
with this setting in a system near the critical point. Also, be aware
that this can really slow down your machine for other tasks if you set
it too high. My recommendation is to figure out how many cores you
have in your system, and start at half of that number.

## Defining Variables for Automatic Sequences of Runs

Variable lines start with "var:" and then have seven parameters that
are all required. Currently only floats can be adjusted. If more than
one variable is configured, the first one will vary fast.

1. The name of the parameter to be varied.
2. The minimum value of the variable.
3. The maximum value of the variable.
4. How many steps to go in.
5. Linear (lin) or geometric (geom) progression.
6. Oscillating (osc) or monotonic (mon).
7. Do we start high (high) or low (low).

For example, the following will vary the temperature and molar volume,
with the molar volume oscillating back and forth between temperature
changes.

var: molarvolume 2.070e-5 6.2115e-5 12 lin osc high
var: temperature 5.556 50.004 14 lin osc high

A default input file is supplied in the top directory of the source
code.
