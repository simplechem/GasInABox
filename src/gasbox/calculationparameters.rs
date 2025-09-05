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

use crate::constants::*;

pub struct CalculationParameters {
    pub molfact: usize,
    pub nmol: usize,
    pub itmax: usize,
    pub anneal: isize,
    pub samplerate: usize,
    pub velscale: f64,
    pub closescale: f64,
    pub nworker: usize,
    pub evalmethod: usize,
    pub accmethod: usize,
}

impl Default for CalculationParameters {
    fn default() -> Self {
        CalculationParameters {
            molfact: MOLFACT,
            nmol: NMOL,
            itmax: ITERMAX,
            anneal: ANNEAL,
            samplerate: SAMPLE,
            velscale: VELSCALE,
            closescale: CLOSESCALE,
            nworker: NWORKER,
            evalmethod: EVALMETHOD,
            accmethod: ACCMETHOD,
        }
    }
}

impl CalculationParameters {
    pub fn _print_status(&self) {
        println!("Running with the following calculation parameters:");
        println!("  molfact: {}", self.molfact);
        println!("  nmol: {}", self.nmol);
        println!("  itmax: {}", self.itmax);
        println!("  anneal: {}", self.anneal);
        println!("  sample rate: {}", self.samplerate);
        println!("  velocity scale: {}", self.velscale);
        println!("  closescale: {}", self.closescale);
        println!("  nworker: {}\n\n", self.nworker);
    }
}
