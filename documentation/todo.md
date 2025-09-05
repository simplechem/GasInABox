# Todo list

Note: This is my todo list. If I recieve requests, I may consider
adding them to the list. At this point I've reached my main goals, so
anything more is "icing on the cake." This is a passion project, not
my real job so ... if you've made a request for a feature and I'm
taking too long, go ahead and fork the project ... just comply with
GPLv2.

## Performance

- Anything I can do to help performance would be nice.
- Maybe learn/add OpenCL code.
- Tweak compile parameters.

## General Features

- Selective compile of SIMD. Use cfg to detect if we are compiling
  with standard or nightly and compile differently so that without
  nightly we still fall back on non-SIMD code variants so things still
  work, even if people don't want to compile with nightly.
- Make better use of Rust. I'm not really using things like traits,
  tests, or lifetimes at the moment. Lifetimes in particular might
  help to make the threading more stable.
- Figure out how to cope with polyatomic molecules (water, CO2,
  ammonia, etc.)
- Figure out how to add multiple types of particles (e.g., mix of Ne
  and Ar)
- Build up a library of parameters for different gases.
- Build up a library of P vs. V_molar phase diagrams and an interface
  to let someone click on a point in the phase diagram to watch that
  condition.
- Ability to load in potential data (e.g., use a cubic spline to
  generate a smooth function from experimental data.)

## File IO

- Set up custom use of the fern crate for logging without colliding
  with the logging system used by Bevy.

## UX/UI

- Design a popup parameters window. Bevy has features for this, I just
  need to sit down and learn how to use them.  - When using periodic
  mode 2, have the spheres slide from one edge to the next in a smooth
  manner, rather than just teleporting. This requires drawing the
  sphere in two placing and clipping at the wall or playing with alpha
  channels. Currently (bevy 0.16) alpha channel solutions give ugly,
  glitchy results.

