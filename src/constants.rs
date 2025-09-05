/******
GasinaBox version 0.1.0: A simulation to visually illustrate the
behavior of ideal gases and interacting gases.

Copyright (c) 2025, Carl Bartels
(University of Manitoba, Department of Chemistry)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see
<https://www.gnu.org/licenses/>.
 *****/

pub const LANES: usize = 4;

pub const BOLTZMANN: f64 = 1.380649e-23;
pub const AVOGADRO: f64 = 6.02214076e23;
pub const AMU: f64 = 1.66053906660e-27;

//Ne epsilon = 3.2135 meV = 5.1486e-13 J, sigma = 0.2782 nm
//Ne Tc = 46.7 K, Pc = 27.7 bar = 2.77e6 Pa
//Ne rhoc = 13.9 kmol/m^3,
//Ne molar mass = 20.1797 amu
//Ne-20 molar mass = 19.9924401753 amu

pub const SIGMA: f64 = 2.782e-10; //m
pub const EPSILONLJ: f64 = 5.1486e-22; //J
pub const MASS: f64 = 19.9924401753 * AMU;

pub const MOLFACT: usize = 4;
pub const NMOL: usize = 1 << (3 * MOLFACT);

pub const ITERMAX: usize = 0;
pub const ANNEAL: isize = 1 << 18;
pub const SAMPLE: usize = 1 << 10;
pub const NTRACK: usize = 1 << 8;
pub const NDIST: usize = 1 << 8;

pub const TCRIT: f64 = 44.448; //Kelvin
pub const PCRIT: f64 = 2.76e6; //Pa
pub const DCRIT: f64 = 2.415e4; //mol/m^3
pub const VCRIT: f64 = 1.0 / DCRIT; // m^3/mol
pub const DCRITA: f64 = DCRIT * AVOGADRO; //atom/m^3
pub const VCRITA: f64 = 1.0 / DCRITA; //m^3/atom
pub const ELEMENT: &str = "Ne-20";

pub const BEVYSCALE: f64 = 1.0e9;
pub const DTSTEPS: usize = 1 << 8;
pub const TEMPERATURE: f64 = 1.0e0 * TCRIT;
pub const CLOSESCALE: f64 = 3.0;
pub const VELSCALE: f64 = 1.6;

//table of EVALMETHODS:
//  0: lennard-jones based on (sigma/r)**2 from r**2
//  1: lennard-jones based on (sigma/r)**6 from r**2
//  2: morse from r**2

pub const EVALMETHOD: usize = 2;

//table of ACCMETHODS:
//  0: naive
//  1: with grid
//  2: threaded

pub const ACCMETHOD: usize = 3;

pub const BASENAME: &str = "result";

pub const NWORKER: usize = 6;
