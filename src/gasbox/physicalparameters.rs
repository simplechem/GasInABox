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
use std::simd::Simd;


pub struct PhysicalParameters {
    pub temperature: f64,
    pub ebeta: f64,
    pub beta: f64,
    pub mass: f64,
    pub massinv: f64,
    pub epsilon: f64,
    pub sigma: f64,
    pub sigma2: f64,
    pub sigma6: f64,
    pub re: f64,
    pub alpha: f64,
    pub upre: f64,
    pub apre: f64,
    pub bound: f64,
    pub bound2: f64,
    pub boundgoal: f64,
    pub bounddir: i32,
    pub boundtweak: f64,
    pub adjust: f64,
    pub potscale: f64,
    pub dt: f64,
    pub dtsplat: Simd<f64, LANES>,
    pub dtsteps: usize,
    pub _vscale: f64,
    pub vmax: f64,
    pub rmax: f64,
    pub rmax2: f64,
}

impl Default for PhysicalParameters {
    fn default() -> Self {
        let dta = SIGMA / (((2.0 / MASS) *
                            TEMPERATURE * BOLTZMANN).sqrt() * DTSTEPS as f64);
        let dtb = SIGMA / (((2.0 / MASS) *
                            EPSILONLJ).sqrt() * DTSTEPS as f64);
        let dt = if dta > dtb { dtb } else { dta };
        let re: f64 = SIGMA * 2.0f64.powf(1.0f64 / 6.0f64);
        let alpha: f64 = 6.0/re;
        let upre: f64 = EPSILONLJ;
        let apre: f64 = 2.0 * EPSILONLJ * alpha * 1.0 / MASS;
        PhysicalParameters {
            temperature: TEMPERATURE,
            ebeta: TEMPERATURE * BOLTZMANN,
            beta: 1.0 / (TEMPERATURE * BOLTZMANN),
            mass: MASS,
            massinv: 1.0 / MASS,
            epsilon: EPSILONLJ,
            sigma: SIGMA,
            sigma2: SIGMA * SIGMA,
            sigma6: SIGMA * SIGMA * SIGMA * SIGMA * SIGMA * SIGMA,
            re,
            alpha,
            upre,
            apre,
            bound: (VCRITA * (NMOL as f64)).cbrt(),
            bound2: (VCRITA * VCRITA * ((NMOL * NMOL) as f64)).cbrt(),
            boundgoal: (VCRITA * (NMOL as f64)).cbrt(),
            bounddir: 0,
            boundtweak: (2.0f64.ln() / (32.0 * 32.0)).exp(),
            adjust: (VCRITA * (NMOL as f64)).cbrt() - 0.01 * SIGMA,
            dtsteps: DTSTEPS,
            dt,
            dtsplat: Simd::<f64, LANES>::splat(dt),
            potscale: 1.0,
            _vscale: (2.0f64.ln() / 16.0).exp(),
            vmax: (5.9*BOLTZMANN*TEMPERATURE/MASS).sqrt(),
            rmax2: CLOSESCALE * CLOSESCALE * SIGMA * SIGMA,
            rmax: CLOSESCALE * SIGMA,
        }
    }
}

impl PhysicalParameters {
    pub fn _print_status(&self) {
        println!("Running with the following physical parameters:");
        println!("  temperature: {} K", self.temperature);
        println!("  ebeta (kT): {} J", self.ebeta);
        println!("  beta (1/kT): {} 1/J", self.beta);
        println!("  mass: {} kg", self.mass);
        println!("  massinv (1/mass): {} 1/kg", self.massinv);
        println!("  epsilon (LJ): {} J", self.epsilon);
        println!("  sigma (LJ): {} m", self.sigma);
        println!("  sigma^6 (LJ): {} m^6", self.sigma6);
        println!("  bound: {} nm", self.bound);
        println!("  adjust: {} nm", self.adjust);
        println!("  potscale: {}", self.potscale);
        println!("  dt: {}", self.dt);
        println!("  vscale: {}\n\n", self._vscale);
        println!("  rmax: {}\n\n", self.rmax);
        println!("  rmax2: {}\n\n", self.rmax2);
    }

    pub fn setup_morse(
        &mut self,
    ) {
        self.re = self.sigma * 2.0f64.powf(1.0f64 / 6.0f64);
        self.alpha = 6.0/self.re;
        self.upre = self.epsilon;
        self.apre = 2.0 * self.epsilon * self.alpha * 1.0 / self.mass;
    }

    pub fn setup_rmax(
        &mut self,
        closescale: f64,
    ) {
        self.rmax2 = closescale * closescale * self.sigma * self.sigma;
        self.rmax = closescale * self.sigma;
    }
    
    pub fn reset_temperature(&mut self, temperature: f64) {
        let dta = self.sigma
            / ((2.0 * self.massinv * temperature * BOLTZMANN).sqrt() * self.dtsteps as f64);
        let dtb = self.sigma / ((2.0 * self.massinv * self.epsilon).sqrt() * self.dtsteps as f64);
        self.dt = if dta > dtb { dtb } else { dta };
        self.ebeta = BOLTZMANN * temperature;
        self.beta = 1.0 / (BOLTZMANN * temperature);
        self.vmax = (5.9*self.ebeta*self.massinv).sqrt();
        self.temperature = temperature;
    }

    pub fn reset_boundgoal(&mut self, boundgoal: f64) {
        if boundgoal > self.bound {
            self.boundgoal = boundgoal;
            self.bounddir = 1;
        } else if boundgoal < self.bound {
            self.boundgoal = boundgoal;
            self.bounddir = -1;
        }
    }

    pub fn bound_from_volume(
        &mut self,
        volume: f64,
        nmol: usize,
    ) -> f64 {
        (volume * (nmol as f64) / AVOGADRO).cbrt()*0.5
    }

    pub fn bound_from_density(
        &mut self,
        density: f64,
        nmol: usize,
    ) -> f64 {
        (nmol as f64 / density / AVOGADRO).cbrt()*0.5
    }
    
    pub fn update_bound(&mut self) {
        if self.bounddir == -1 {
            self.bound /= self.boundtweak;
            if self.bound < self.boundgoal {
                self.bounddir = 0;
                self.bound = self.boundgoal;
            }
        } else if self.bounddir == 1 {
            self.bound *= self.boundtweak;
            if self.bound > self.boundgoal {
                self.bounddir = 0;
                self.bound = self.boundgoal;
            }
        }
        self.adjust = self.bound - 0.01 * self.sigma;
        self.bound2 = self.bound * self.bound;
    }

    pub fn reset_dtsteps(&mut self, dtsteps: usize) {
        let dta = self.sigma
            / ((2.0 * self.massinv * self.temperature * BOLTZMANN).sqrt() * dtsteps as f64);
        let dtb = self.sigma / ((2.0 * self.massinv * self.epsilon).sqrt() * dtsteps as f64);
        self.dt = if dta > dtb { dtb } else { dta };
        self.dtsteps = dtsteps;
        self.dtsplat = Simd::<f64, LANES>::splat(self.dt);
    }
}
