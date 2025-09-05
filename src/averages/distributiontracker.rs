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
use std::fs::OpenOptions;
use std::io::Write;

pub struct DistributionTracker {
    dx: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
    peakx: f64,
    pub nval: usize,
    values: Vec<u128>,
    pub normalized: Vec<f64>,
    count: u128,
}

impl Default for DistributionTracker {
    fn default() -> Self {
        DistributionTracker {
            dx: 1.0,
            xmax: 0.0,
            ymin: 0.0,
            ymax: 0.0,
            peakx: 0.0,
            nval: NDIST,
            values: vec![0; NDIST],
            normalized: vec![0.0; NDIST+5],
            count: 0,
        }
    }
}

impl DistributionTracker {
    pub fn insert_value (
        &mut self,
        value: f64,
    ) {
        let i: usize = (value/self.dx).floor() as usize;
        if i < self.nval {
            self.values[i] += 1;
            self.count += 1;
        }
        if value > self.xmax {
            self.xmax = value;
        }
    }

    pub fn normalize (
        &mut self,
        norm: (usize, f64),
    ) {
        match norm.0 {
            0 => {
                self.normalize_flat(norm.1);
            }
            1 => {
                self.normalize_r2(norm.1);
            }
            _ => {
                self.normalize_direct(norm.1);
            }
        }
        self.normalized[self.nval] = self.dx;
        self.normalized[self.nval+1] = self.dx * self.nval as f64;
        self.normalized[self.nval+2] = self.ymin;
        self.normalized[self.nval+3] = self.ymax;
        self.normalized[self.nval+4] = self.peakx;
    }

    fn normalize_flat(
        &mut self,
        norm: f64,
    ) {
        let scale: f64 = norm / (self.dx * (self.count as f64));
        self.ymin = scale * self.values[0] as f64;
        self.ymax = self.ymin;
        (0..).take(self.nval).for_each(|i| {
            self.normalized[i] = scale * self.values[i] as f64;
            if self.normalized[i] < self.ymin {
                self.ymin = self.normalized[i];
            } else if self.normalized[i] > self.ymax {
                self.ymax = self.normalized[i];
                self.peakx = self.dx * i as f64;
            }
        });
    }

    fn normalize_r2(
        &mut self,
        norm: f64,
    ) {
        let scale: f64 = norm / (self.dx * self.dx * self.dx *
                                   (self.count as f64));
        self.ymin = 0.0;
        self.ymax = 0.0;
        self.normalized[0] = 0.0;
        (0..).take(self.nval - 1).for_each(|i| {
            self.normalized[i] = scale * (self.values[i] as f64) /
                ((i*i) as f64);
            if self.normalized[i] < self.ymin {
                self.ymin = self.normalized[i];
            } else if self.normalized[i] > self.ymax {
                self.ymax = self.normalized[i];
                self.peakx = self.dx * i as f64;
            }
        });
    }

    fn normalize_direct(
        &mut self,
        norm: f64,
    ) {
        self.ymin = self.values[0] as f64;
        self.ymax = self.ymin;
        (0..).take(self.nval).for_each(|i| {
            self.normalized[i] = norm * self.values[i] as f64;
            if self.ymin < self.normalized[i] {
                self.ymin = self.normalized[i];
            } else if self.ymax > self.normalized[i] {
                self.ymax = self.normalized[i];
                self.peakx = self.dx * i as f64;
            }
        });
    }

    fn clear(
        &mut self,
    ) {
        self.count = 0;
        self.xmax = 0.0;
        self.values = vec![0; self.nval];
        self.normalized = vec![0.0; self.nval+5];
    }

    pub fn reset_dx(
        &mut self,
        dx: f64,
    ) {
        self.clear();
        self.dx = dx;
    }

    pub fn set_dx_from_max(
        &mut self,
    ) {
        self.dx = self.xmax / (self.nval as f64);
        self.clear();
    }

    pub fn dump_distribution(
        &mut self,
        filename: String,
    ) {
        let mut file = OpenOptions::new()
            .append(true)
            .create(true)
            .open(filename)
            .expect("trouble openning the distribution dump file.");

        //note that we need to normalize before we dump.

        for i in (0..).take(self.nval) {
            writeln!(file, "{}, {}",
                     self.dx*(i as f64),
                     self.normalized[i],
            )
                .expect("trouble writing to distribution file");
        }
    }
}
