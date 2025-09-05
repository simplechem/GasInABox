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

pub struct GraphParameters {
    pub dx: f32,
    pub xmax: f32,
    pub ymin: f32,
    pub ymax: f32,
    pub yrng: f32,
    pub peakx: f32,
}

impl Default for GraphParameters {
    fn default() -> Self {
        GraphParameters {
            dx: 0.1,
            xmax: 1.0,
            ymin: 0.0,
            ymax: 1.0,
            yrng: 1.0,
            peakx: 0.5,
        }
    }
}

