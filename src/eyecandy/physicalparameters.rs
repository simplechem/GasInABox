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
pub struct PhysicalParameters {
    pub nmol: usize,
    pub bound: f32,
    pub bevyscale: f64,
    pub ballradius: f32,
    pub ballcoords: Vec<[f32; 3]>,
}

impl Default for PhysicalParameters {
    fn default() -> Self {
        PhysicalParameters {
            nmol: NMOL,
            bound: (BEVYSCALE * (VCRITA * (NMOL as f64)).cbrt()) as f32,
            bevyscale: BEVYSCALE,
            ballradius: (1.0 * BEVYSCALE * SIGMA) as f32,
            ballcoords: vec![[0.0f32;3];NMOL],
        }
    }
}

