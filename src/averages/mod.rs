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

use crate::channelstructs::*;
use crate::constants::*;

use std::fs::OpenOptions;
use std::io::Write;

pub mod distributiontracker;

//use crate::averages::distributiontracker::*;

pub struct AverageTracker {
    pub itcount: isize,
    pub ntrack: usize,
    pub ntrackinv: f64,
    pub kinetic: f64,
    pub kinetic1mean: f64,
    pub kinetic2mean: f64,
    pub kinetic1track: Vec<f64>,
    pub kinetic2track: Vec<f64>,
    pub potential: f64,
    pub potential1mean: f64,
    pub potential2mean: f64,
    pub potential1track: Vec<f64>,
    pub potential2track: Vec<f64>,
    pub pressure: f64,
    pub pressure1mean: f64,
    pub pressure2mean: f64,
    pub pressure1track: Vec<f64>,
    pub pressure2track: Vec<f64>,
    pub deltaqdeltat: f64,
    pub deltaqdeltat1mean: f64,
    pub deltaqdeltat2mean: f64,
    pub deltaqdeltat1track: Vec<f64>,
    pub deltaqdeltat2track: Vec<f64>,
    pub itspersec: f64,
    pub itspersec1mean: f64,
    pub itspersec2mean: f64,
    pub itspersec1track: Vec<f64>,
    pub itspersec2track: Vec<f64>,
    pub count: usize,
}

impl Default for AverageTracker {
    fn default() -> Self {
        AverageTracker {
            itcount: 0,
            ntrack: NTRACK,
            ntrackinv: 1.0 / NTRACK as f64,
            kinetic: 0.0,
            kinetic1mean: 0.0,
            kinetic2mean: 0.0,
            kinetic1track: vec![0.0f64; NTRACK],
            kinetic2track: vec![0.0f64; NTRACK],
            potential: 0.0,
            potential1mean: 0.0,
            potential2mean: 0.0,
            potential1track: vec![0.0f64; NTRACK],
            potential2track: vec![0.0f64; NTRACK],
            pressure: 0.0,
            pressure1mean: 0.0,
            pressure2mean: 0.0,
            pressure1track: vec![0.0f64; NTRACK],
            pressure2track: vec![0.0f64; NTRACK],
            deltaqdeltat: 0.0,
            deltaqdeltat1mean: 0.0,
            deltaqdeltat2mean: 0.0,
            deltaqdeltat1track: vec![0.0f64; NTRACK],
            deltaqdeltat2track: vec![0.0f64; NTRACK],
            itspersec: 0.0,
            itspersec1mean: 0.0,
            itspersec2mean: 0.0,
            itspersec1track: vec![0.0f64; NTRACK],
            itspersec2track: vec![0.0f64; NTRACK],
            count: 0,
        }
    }
}

impl AverageTracker {
    pub fn update_averages(
        &mut self,
        data: AverageRawData,
    ) {
        self.itcount = data.itcount;

        self.kinetic = data.kinetic;
        self.kinetic1mean +=
            self.ntrackinv * (self.kinetic - self.kinetic1track[self.count]) ;
        self.kinetic2mean +=
            self.ntrackinv * (self.kinetic * self.kinetic -
                              self.kinetic2track[self.count]);
        self.kinetic1track[self.count] = self.kinetic;
        self.kinetic2track[self.count] = self.kinetic * self.kinetic;

        self.potential = data.potential;
        self.potential1mean +=
            self.ntrackinv * (self.potential - self.potential1track[self.count]) ;
        self.potential2mean +=
            self.ntrackinv * (self.potential * self.potential -
                              self.potential2track[self.count]);
        self.potential1track[self.count] = self.potential;
        self.potential2track[self.count] = self.potential * self.potential;

        self.pressure = data.pressure;
        self.pressure1mean +=
            self.ntrackinv * (self.pressure - self.pressure1track[self.count]) ;
        self.pressure2mean +=
            self.ntrackinv * (self.pressure * self.pressure -
                              self.pressure2track[self.count]);
        self.pressure1track[self.count] = self.pressure;
        self.pressure2track[self.count] = self.pressure * self.pressure;

        self.deltaqdeltat = data.deltaqdeltat;
        self.deltaqdeltat1mean +=
            self.ntrackinv * (self.deltaqdeltat - self.deltaqdeltat1track[self.count]) ;
        self.deltaqdeltat2mean +=
            self.ntrackinv * (self.deltaqdeltat * self.deltaqdeltat -
                              self.deltaqdeltat2track[self.count]);
        self.deltaqdeltat1track[self.count] = self.deltaqdeltat;
        self.deltaqdeltat2track[self.count] = self.deltaqdeltat * self.deltaqdeltat;

        self.itspersec = data.itspersec;
        self.itspersec1mean +=
            self.ntrackinv * (self.itspersec - self.itspersec1track[self.count]) ;
        self.itspersec2mean +=
            self.ntrackinv * (self.itspersec * self.itspersec -
                              self.itspersec2track[self.count]);
        self.itspersec1track[self.count] = self.itspersec;
        self.itspersec2track[self.count] = self.itspersec * self.itspersec;

        self.count += 1;
        if self.count >= self.ntrack {
            self.count = 0;
        }
    }

    pub fn average_result(
        &self,
    ) -> AverageResult {
        let result: AverageResult = AverageResult {
            kinetic1mean: self.kinetic1mean,
            kinetic2mean: self.kinetic2mean,
            potential1mean: self.potential1mean,
            potential2mean: self.potential2mean,
            pressure1mean: self.pressure1mean,
            pressure2mean: self.pressure2mean,
            deltaqdeltat1mean: self.deltaqdeltat1mean,
            deltaqdeltat2mean: self.deltaqdeltat2mean,
            itspersec1mean: self.itspersec1mean,
            itspersec2mean: self.itspersec2mean,
            ntrackinv: self.ntrackinv,
        };
        result
    }

    pub fn reset_tracks(
        &mut self,
        ntrack: usize,
    ) {
        self.ntrack = ntrack;
        self.ntrackinv = 1.0/ntrack as f64;
        self.kinetic1track = vec![0.0f64; ntrack];
        self.kinetic2track = vec![0.0f64; ntrack];
        self.potential1track = vec![0.0f64; ntrack];
        self.potential2track = vec![0.0f64; ntrack];
        self.pressure1track = vec![0.0f64; ntrack];
        self.pressure2track = vec![0.0f64; ntrack];
        self.deltaqdeltat1track = vec![0.0f64; ntrack];
        self.deltaqdeltat2track = vec![0.0f64; ntrack];
        self.itspersec1track = vec![0.0f64; ntrack];
        self.itspersec2track = vec![0.0f64; ntrack];
    }

    pub fn dump(
        &self,
        temperature: f64,
        molarvol: f64,
        filename: String,
    ) {
        let mut file = OpenOptions::new()
            .append(true)
            .create(true)
            .open(filename)
            .expect("trouble openning the tracking file.");

        writeln!(file, "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
                 self.itcount,
                 self.ntrack,
                 self.kinetic,
                 self.kinetic1mean,
                 self.kinetic2mean,
                 self.potential,
                 self.potential1mean,
                 self.potential2mean,
                 self.pressure,
                 self.pressure1mean,
                 self.pressure2mean,
                 self.deltaqdeltat,
                 self.deltaqdeltat1mean,
                 self.deltaqdeltat2mean,
                 self.itspersec,
                 self.itspersec1mean,
                 self.itspersec2mean,
                 temperature,
                 molarvol,
        )
            .expect("trouble writing average data to tracking file");
    }
}
