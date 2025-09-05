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

use bevy::prelude::*;

use crate::constants::*;

#[derive(Resource)]
pub struct TextParameters {
    pub itcount: isize,
    pub ittotal: isize,
    pub elapsedtime: f64,
    pub timetotal: f64,
    pub itmax: usize,
    pub samples: usize,
    pub element: String,
    pub nmol: usize,
    pub temperature: f64,
    pub volume: f64,
    pub dt: f64,
    pub kinetic: f64,
    pub kinetic1mean: f64,
    pub kinetic2mean: f64,
    pub potential: f64,
    pub potential1mean: f64,
    pub potential2mean: f64,
    pub pressure: f64,
    pub pressure1mean: f64,
    pub pressure2mean: f64,
    pub deltaqdeltat: f64,
    pub deltaqdeltat1mean: f64,
    pub deltaqdeltat2mean: f64,
    pub itspersec: f64,
    pub itspersec1mean: f64,
    pub itspersec2mean: f64,
    pub ntrackinv: f64,
    pub displayactive: bool,
}

impl Default for TextParameters {
    fn default() -> Self {
        TextParameters {
            itcount: 0,
            ittotal: 0,
            elapsedtime: 0.0f64,
            timetotal: 0.0f64,
            itmax: 0,
            samples: 0,
            element: "Ne-20\nTc = 44.448 K\nPc = 2.76 MPa\nVc = 41.7 mL/mol".to_string(),
            nmol: NMOL,
            temperature: 0.0,
            volume: 0.0,
            dt: 0.0,
            kinetic: 0.0,
            kinetic1mean: 0.0,
            kinetic2mean: 0.0,
            potential: 0.0,
            potential1mean: 0.0,
            potential2mean: 0.0,
            pressure: 0.0,
            pressure1mean: 0.0,
            pressure2mean: 0.0,
            deltaqdeltat: 0.0,
            deltaqdeltat1mean: 0.0,
            deltaqdeltat2mean: 0.0,
            itspersec: 0.0,
            itspersec1mean: 0.0,
            itspersec2mean: 0.0,
            ntrackinv: 1.0,
            displayactive: true,
        }
    }
}
