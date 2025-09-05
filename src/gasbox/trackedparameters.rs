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

pub struct TrackedParameters {
    pub kinetic: f64,
    pub potential: f64,
    pub impulse: f64,
    pub pressure: f64,
    pub deltaq: f64,
    pub deltaqdeltat: f64,
    pub itspersec: f64,
    pub itcount: isize,
    pub ittotal: isize,
    pub elapsedtime: f64,
    pub timetotal: f64,
    pub displayactive: bool,
}

impl Default for TrackedParameters {
    fn default() -> Self {
        TrackedParameters {
            kinetic: 0.0,
            potential: 0.0,
            impulse: 0.0,
            pressure: 0.0,
            deltaq: 0.0,
            deltaqdeltat: 0.0,
            itspersec: 0.0,
            itcount: -ANNEAL,
            ittotal: 0,
            elapsedtime: 0.0f64,
            timetotal: 0.0f64,
            displayactive: true,
        }
    }
}

impl TrackedParameters {
    pub fn update_extrinsics (
        &mut self,
    ) {
        self.pressure = self.impulse ;
        self.deltaqdeltat = self.deltaq ;
        self.impulse = 0.0 ;
        self.deltaq = 0.0 ;
    }

    pub fn _print_status(&self) {
        println!("Current state of tracked parameters:");
        println!("  kinetic energy (kT/particle): {}", self.kinetic);
        println!("  potential energy (kT/particle): {}", self.potential);
        println!("  impulse (): {}", self.impulse);
        println!("  pressure (Pa): {}", self.pressure);
        println!("  delta q (J/sample): {}", self.deltaq);
        println!("  delta q / delta t (J/s): {}", self.deltaqdeltat);
        println!("  iterations per second (1/s): {}", self.itspersec);
        println!("  itcount (this case): {}", self.itcount);
        println!("  ittotal (total run): {}", self.ittotal);
        println!("  elapsed time (this case, s): {}", self.elapsedtime);
        println!("  total time (this run, s): {}\n\n", self.timetotal);
    }
}

