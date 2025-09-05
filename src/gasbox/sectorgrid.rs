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

pub struct SectorGrid {
    pub n1d: usize,
    pub n2d: usize,
    pub n3d: usize,
    width: f64,
    bound: f64,
    pub periodic: usize,
    pub sectors: Vec<Vec<usize>>,
    pub crosslist: Vec<(usize, usize)>,
    pub layerbounds: Vec<(usize, usize)>,
    pub linebounds: Vec<(usize, usize)>,
    pub blockbounds: Vec<(usize, usize)>,
}

impl Default for SectorGrid {
    fn default() -> Self {
        let bound: f64 = (VCRITA * (NMOL as f64)).cbrt();
        let n1d: usize = (2.0 * bound / (CLOSESCALE * SIGMA)).floor() as usize;
        let width: f64 = 2.0 * bound / n1d as f64;
        SectorGrid {
            n1d,
            n2d: n1d * n1d,
            n3d: n1d * n1d * n1d,
            width,
            bound,
            periodic: 0,
            sectors: vec![
                Vec::<_>::with_capacity(NMOL / (n1d * n1d * n1d));
                n1d * n1d * n1d
            ],
            crosslist: vec![],
            layerbounds: vec![],
            linebounds: vec![],
            blockbounds: vec![],
        }
    }
}

impl SectorGrid {
    pub fn generate_crosslist_sequence(
        &mut self,
    ) {
        self.crosslist.clear();
        self.layerbounds.clear();
        self.linebounds.clear();
        self.blockbounds.clear();

        //Trust me or work it out for yourself. This is set up so that
        // all possible neighbor sectors are found, but not
        // duplicated. There is a secondary test that deals with
        // neighbors at the ends differently if we are periodic or not
        // periodic.
        
        let offsets: Vec<[isize;3]> = vec![
            [ 0, 0, 0],
            [ 1, 0, 0],
            [ 0, 1, 0],
            [ 1, 1, 0],
            [-1, 1, 0],
            [ 0, 0, 1],
            [ 1, 0, 1],
            [ 1, 1, 1],
            [ 0, 1, 1],
            [-1, 1, 1],
            [-1, 0, 1],
            [-1,-1, 1],
            [ 0,-1, 1],
            [ 1,-1, 1],
        ];
        
        for i in (0..).take(self.n3d) {
            let z1: isize = (i / self.n2d) as isize ;
            let y1: isize = (i as isize - z1 * (self.n2d as isize)) /
                self.n1d as isize ;
            let x1: isize = i as isize - z1 * (self.n2d as isize) -
                y1 * (self.n1d as isize) ;
            for j in (0..).take(offsets.len()) {
                let (xd, yd, zd): (isize, isize, isize) = (
                    offsets[j][0],
                    offsets[j][1],
                    offsets[j][2],
                ) ;
                let x2: isize = self.edge_check(x1, xd) ;
                let y2: isize = self.edge_check(y1, yd) ;
                let z2: isize = self.edge_check(z1, zd) ;
                let k = self.position_from_xyz((x2, y2, z2));
                if self.periodic > 0 ||
                    (
                        (x2-x1).abs() <= 1 &&
                            (y2-y1).abs() <= 1 &&
                            (z2-z1).abs() <= 1
                    )
                {
                    if i <= k as usize {
                        self.crosslist.push((i,k as usize));
                    } else {
                        self.crosslist.push((k as usize,i));
                    }
                }
            }
        }
                                
        self.prune_redundant();
        
        self.fix_layerbounds();
    }

    fn edge_check(
        &self,
        a: isize,
        d: isize,
    ) -> isize {
        let n1d1: isize = self.n1d as isize - 1;
        let b: isize ;
        if 0 == a && -1 == d {
            b = n1d1 ;
        } else if n1d1 == a && 1 == d {
            b = 0 ;
        } else {
            b = a + d ;
        }

        b
    }

    fn fix_layerbounds(
        &mut self,
    ) {
        let mut j: usize = 0;
        let mut k: usize = 0;
        let mut l: usize = 0;
        self.layerbounds.push((j,0));
        self.linebounds.push((k,0));
        self.blockbounds.push((l,0));

        for i in 0..self.crosslist.len() {
            if self.crosslist[i].0 == j + self.n2d {
                j += self.n2d;
                self.layerbounds.push((i,0));
                let length = self.layerbounds.len();
                let l0 = self.layerbounds[length-2].0;
                let l1 = self.layerbounds[length-1].0;
                self.layerbounds[length-2] = (l0, l1);
            }
            if self.crosslist[i].0 == k + self.n1d {
                k += self.n1d;
                self.linebounds.push((i,0));
                let length = self.linebounds.len();
                let l0 = self.linebounds[length-2].0;
                let l1 = self.linebounds[length-1].0;
                self.linebounds[length-2] = (l0, l1);
            }
            if self.crosslist[i].0 == l + 1 {
                l += 1;
                self.blockbounds.push((i,0));
                let length = self.blockbounds.len();
                let l0 = self.blockbounds[length-2].0;
                let l1 = self.blockbounds[length-1].0;
                self.blockbounds[length-2] = (l0, l1);
            }
        }
        
        let mut length = self.layerbounds.len();
        self.layerbounds[length-1].1 = self.crosslist.len();
        length = self.linebounds.len();
        self.linebounds[length-1].1 = self.crosslist.len();
        length = self.blockbounds.len();
        self.blockbounds[length-1].1 = self.crosslist.len();
    }

    fn fix_all_tupples(
        &mut self,
    ) {
        for i in (0..).take(self.crosslist.len()) {
            if self.crosslist[i].0 > self.crosslist[i].1 {
                (self.crosslist[i].0, self.crosslist[i].1) =
                    (self.crosslist[i].1, self.crosslist[i].0);
            }
        }
    }

    fn prune_redundant(
        &mut self,
    ) {
        self.fix_all_tupples();
        self.crosslist.sort_by(|a, b| a.1.cmp(&b.1));
        self.crosslist.sort_by(|a, b| a.0.cmp(&b.0));
        self.crosslist.dedup();

        let mut tmplist: Vec<(usize, usize)> = vec![];
        tmplist.push(self.crosslist[0]);
        let mut j: usize = 0;
        for i in (1..).take(self.crosslist.len()-1) {
            if tmplist[j].0 != self.crosslist[i].0 ||
                tmplist[j].1 != self.crosslist[i].1 {
                    tmplist.push(self.crosslist[i]);
                    j += 1;
                }
        }
        self.crosslist.clear();

        for i in (0..).take(tmplist.len()) {
            self.crosslist.push(tmplist[i]);
        }        
    }

    fn position_from_xyz(
        &self,
        xyz: (isize, isize, isize),
    ) -> isize {
        let w: isize = xyz.0 +
            xyz.1 * self.n1d as isize +
            xyz.2 * self.n2d as isize;
        w
    }

    fn _report_status(
        &self,
    ) {
        println!("we have {} sectors", self.sectors.len());
        println!("crosslist contains {} elements", self.crosslist.len());
        println!("n1d = {}, n2d = {}, n3d = {}", self.n1d, self.n2d, self.n3d);
        println!("layerbounds contains {} elements", self.layerbounds.len());
        println!("linebounds contains {} elements", self.linebounds.len());
        println!("blockbounds contanis {} elements", self.blockbounds.len());
    }

    fn _valid_check(
        &self,
        m: isize,
        n: isize,
    ) -> bool {
        if m > n || n >= self.n3d as isize || n < 0 {
            return false;
        }
        let n2d: isize = self.n2d as isize;
        let n1d: isize = self.n1d as isize;
        let n1d1: isize = n1d - 1;
        let nz: isize = n / n2d;
        let ny: isize = (n - nz * n2d) / n1d;
        let nx: isize = n - nz * n2d - ny * n1d;
        let mz: isize = m / n2d;
        let my: isize = (m - mz * n2d) / n1d;
        let mx: isize = m - mz * n2d - my * n1d;
        if (nx - mx).abs() > 1 && (nx - mx).abs() != n1d1
        {
            return false;
        }
        if (ny - my).abs() > 1 && (ny - my).abs() != n1d1
        {
            return false;
        }
        if (nz - mz).abs() > 1 && (nz - mz).abs() != n1d1
        {
            return false;
        }
        true
    }

    pub fn rescale_grid(
        &mut self,
        bound: f64,
        closescale: f64,
        sigma: f64,
    ) {
        let n1dnew = (2.0 * bound / (closescale * sigma)).floor() as usize;
        self.width = 2.0 * bound / n1dnew as f64;
        self.bound = bound;
        self.n1d = n1dnew;
        self.n2d = self.n1d * self.n1d;
        self.n3d = self.n1d * self.n2d;
        
        self.sectors = vec![
            Vec::<_>::with_capacity(NMOL / (self.n3d));
            self.n3d
        ];
        self.crosslist.clear();
        self.generate_crosslist_sequence();
    }

    pub fn add_index(
        &mut self,
        index: usize,
        coords: &Simd<f64, LANES>,
    ) {
        let i: usize = ((coords[0] + self.bound) / self.width).floor() as usize;
        let j: usize = ((coords[1] + self.bound) / self.width).floor() as usize;
        let k: usize = ((coords[2] + self.bound) / self.width).floor() as usize;

        assert!(i < self.n1d, "i = {} of {}", i, self.n1d);
        assert!(j < self.n1d, "j = {} of {}", j, self.n1d);
        assert!(k < self.n1d, "k = {} of {}", k, self.n1d);
        
        self.sectors[i + j * self.n1d + k * self.n2d].push(index);
    }

    pub fn clear_sectors(
        &mut self,
    ) {
        (0..self.n3d).for_each(|i| {
            self.sectors[i].clear();
        });  
    }
}
