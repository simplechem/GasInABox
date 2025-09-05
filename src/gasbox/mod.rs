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

use crate::averages::*;
use crate::channelstructs::*;
use crate::constants::*;
use crate::fileio::*;

mod calculationparameters;
mod physicalparameters;
mod trackedparameters;
mod sectorgrid;

use crate::gasbox::calculationparameters::*;
use crate::gasbox::physicalparameters::*;
use crate::gasbox::trackedparameters::*;
use crate::gasbox::sectorgrid::*;
use crate::averages::distributiontracker::*;
    
use crossbeam_channel::{unbounded, Receiver, Sender, TryRecvError};
use rand::{thread_rng, Rng};
use std::process;
use std::time::Instant;
use std::path::Path;

use std::fs;

use std::simd::num::SimdFloat;
use std::simd::Simd;

use std::thread;

type ThreadSender = Sender<(usize,
                            MinThreadParameters,
                            MinThreadPointers)>;
type ThreadReceiver = Receiver<(usize,
                                MinThreadParameters,
                                MinThreadPointers)>;

pub struct GasBox {
    cparms: CalculationParameters,
    pparms: PhysicalParameters,
    tparms: TrackedParameters,
    grid: SectorGrid,
    coords: Vec<Simd<f64, LANES>>,
    vels: Vec<Simd<f64, LANES>>,
    accs: Vec<Simd<f64, LANES>>,
    pots: Vec<f64>,
    averages: AverageTracker,
    speeds: DistributionTracker,
    radialdist: DistributionTracker,
    coords_requ_channel: Receiver<usize>,
    coords_send_channel: Sender<TransferCoords>,
    text_requ_channel: Receiver<usize>,
    text_send_channel: Sender<MacroStateData>,
    dist_requ_channel: Receiver<usize>,
    dist_send_channel: Sender<TransferDist>,
    work_requ_channels: Vec<Sender<(usize,
                                    MinThreadParameters,
                                    MinThreadPointers)>>,
    work_rslt_channel: Receiver<usize>,
    work_handles: Vec<thread::JoinHandle<()>>,
    config: SimConfig,
}

impl Default for GasBox {
    fn default() -> Self {
        GasBox {
            cparms: CalculationParameters::default(),
            pparms: PhysicalParameters::default(),
            tparms: TrackedParameters::default(),
            grid: SectorGrid::default(),
            coords: vec![Simd::<f64, LANES>::splat(0.0f64); NMOL],
            vels: vec![Simd::<f64, LANES>::splat(0.0f64); NMOL],
            accs: vec![Simd::<f64, LANES>::splat(0.0f64); NMOL],
            pots: vec![0.0f64; NMOL],
            averages: AverageTracker::default(),
            speeds: DistributionTracker::default(),
            radialdist: DistributionTracker::default(),
            coords_requ_channel: unbounded().1.clone(),
            coords_send_channel: unbounded().0.clone(),
            text_requ_channel: unbounded().1.clone(),
            text_send_channel: unbounded().0.clone(),
            dist_requ_channel: unbounded().1.clone(),
            dist_send_channel: unbounded().0.clone(),
            work_requ_channels: vec![],
            work_rslt_channel: unbounded().1.clone(),
            work_handles: vec![],
            config: SimConfig::default(),
        }
    }
}

pub struct GasBoxFinish {
    pub coords_requ_channel: Receiver<usize>,
    pub coords_send_channel: Sender<TransferCoords>,
    pub text_requ_channel: Receiver<usize>,
    pub text_send_channel: Sender<MacroStateData>,
    pub dist_requ_channel: Receiver<usize>,
    pub dist_send_channel: Sender<TransferDist>,
    pub config: SimConfig,
}

impl GasBox {
    pub fn finish_setup(
        &mut self,
        finish_input: GasBoxFinish,
    ) {
        self.coords_requ_channel = finish_input.coords_requ_channel.clone();
        self.coords_send_channel = finish_input.coords_send_channel.clone();
        self.text_requ_channel = finish_input.text_requ_channel.clone();
        self.text_send_channel = finish_input.text_send_channel.clone();
        self.dist_requ_channel = finish_input.dist_requ_channel.clone();
        self.dist_send_channel = finish_input.dist_send_channel.clone();

        self.config = finish_input.config.clone();

        self.update_from_simconfig();

        let bound: f64 = self.pparms.bound_from_volume(
            self.config.molarvolume,
            self.cparms.nmol,
        );

        self.pparms.boundgoal = bound ;
        self.pparms.bound = self.pparms.boundgoal;
        self.pparms.adjust = self.pparms.bound - 0.01 * self.pparms.sigma;
        
        self.initialize_coordinates();
        self.initialize_velocities();
        self.zero_acc_pot();
        self.grid.periodic = finish_input.config.periodic;

        self.grid.rescale_grid(
            2.0 * self.pparms.bound,
            self.cparms.closescale,
            self.pparms.sigma,
        );

        self.grid.rescale_grid(
            self.pparms.bound,
            self.cparms.closescale,
            self.pparms.sigma,
        );

        self.averages.reset_tracks(self.config.ntrack);

        let dspeed: f64 = ((2.0 * self.pparms.ebeta *
            self.pparms.massinv).sqrt()) * 0.01;
        self.speeds.reset_dx(dspeed);
        self.radialdist.reset_dx(self.pparms.bound * 0.01);
        self.setup_workers();
    }

    fn update_from_simconfig(
        &mut self,
    ) {
        self.cparms.molfact = self.config.molfact;
        self.cparms.nmol = 1 << (3 * self.config.molfact);
        self.coords = vec![Simd::<f64, LANES>::splat(0.0f64); self.cparms.nmol];
        self.vels = vec![Simd::<f64, LANES>::splat(0.0f64); self.cparms.nmol];
        self.accs = vec![Simd::<f64, LANES>::splat(0.0f64); self.cparms.nmol];
        self.pots = vec![0.0f64; self.cparms.nmol];
        self.cparms.itmax = self.config.itermax;
        self.cparms.anneal = self.config.anneal;
        self.tparms.itcount = -self.config.anneal;
        self.cparms.samplerate = self.config.sample;
        self.cparms.velscale = self.config.velscale;
        self.cparms.closescale = self.config.closescale;
        self.cparms.nworker = self.config.nworker;
        self.cparms.evalmethod = self.config.evalmethod;
        self.cparms.accmethod = self.config.accmethod;
        self.pparms.sigma = self.config.sigma;
        self.pparms.epsilon = self.config.epsilon;
        self.pparms.mass = self.config.mass;
        self.pparms.massinv = 1.0/self.config.mass;
        self.pparms.dtsteps = self.config.dtsteps;
        self.pparms.reset_temperature(self.config.temperature);
        self.pparms.setup_morse();
        self.pparms.setup_rmax(
            self.config.closescale,
        );

        if self.config.readin {
            for i in (0..).take(self.cparms.nmol) {
                self.coords[i] = self.config.coords[i];
                self.vels[i] = self.config.vels[i];
            }
            self.config.coords.clear();
            self.config.vels.clear();
        }

        let bound: f64 = self.pparms.bound_from_volume(
            self.config.molarvolume,
            self.cparms.nmol,
        );
        self.pparms.reset_boundgoal(bound);

        self.tparms.displayactive = self.config.displayactive;
    }
    
    fn zero_acc_pot(
        &mut self,
    ) {
        self.accs.iter_mut().for_each(|item| {
            *item = Simd::<f64, LANES>::splat(0.0f64);
        });
        self.pots.iter_mut().for_each(|item| {
            *item = 0.0f64;
        });
    }

    fn initialize_velocities(
        &mut self,
    ) {
        let vel_max: f64 = (self.pparms.ebeta * self.pparms.massinv).sqrt();
        let vel_range = -self.cparms.velscale * vel_max ..
            self.cparms.velscale * vel_max;
        let mut rng = thread_rng();
        (0..).take(self.cparms.nmol).for_each(|i| {
            for k in (0..).take(3usize) {
                self.vels[i][k] = rng.gen_range(vel_range.clone());
            }
        });
    }
    
    fn initialize_coordinates(
        &mut self,
    ) {
        let scale: f64 = 0.95;
        let perrow: usize = (self.cparms.nmol as f64).cbrt().floor() as usize;
        let dr: f64 = 2.0 * self.pparms.bound / (perrow as f64);
        let rstart: f64 = 0.5 * dr - self.pparms.bound;
        let mut r: [f64; 3] = [rstart, rstart, rstart];
        (0..).take(self.cparms.nmol).for_each(|i| {
            self.coords[i] = Simd::<f64, LANES>::splat(0.0f64);
            (0..).take(3usize).for_each(|k| {
                self.coords[i][k] = scale * r[k];
            });
            r[0] += dr;
            if r[0] > self.pparms.bound {
                r[0] = rstart;
                r[1] += dr;
                if r[1] > self.pparms.bound {
                    r[1] = rstart;
                    r[2] += dr;
                }
            }
        });
    }

    fn update_accs_pots(
        &mut self,
    ) {
        match self.cparms.accmethod {
            0 => {
                self.update_accs_pots_naive();
            }
            1 => {
                self.update_accs_pots_grid();
            }
            2 => {
                self.update_accs_pots_thread();
            }
            3 => {
                self.update_accs_pots_thread();
            }
            4 => {
                self.update_accs_pots_thread();
            }
            _ => {
                println!("bad update method for accs and pots");
                process::exit(1);
            }
        }
    }

    fn update_accs_pots_naive(
        &mut self,
    ) {
        self.zero_acc_pot();
        for i in (0..).take(self.cparms.nmol-1) {
            for j in ((i+1)..).take(self.cparms.nmol-i-1) {
                self.accs_pots_eval_pair_mut(i, j);
            }
        }
    }

    fn update_accs_pots_grid(
        &mut self,
    ) {
        self.zero_acc_pot();
        self.grid.clear_sectors();
        (0..).take(self.cparms.nmol).for_each(|i| {
            self.grid.add_index(i, &self.coords[i]);
        });

        let length: usize = self.grid.crosslist.len();
        (0..).take(length).for_each(|l| {
            if self.grid.crosslist[l].0 == self.grid.crosslist[l].1 {
                self.update_accs_pots_grid_self(l);
            } else {
                self.update_accs_pots_grid_cross(l);
            }
        });
    }

    fn update_accs_pots_grid_self(
        &mut self,
        l: usize,
    ) {
        let s: usize = self.grid.crosslist[l].0;
        let length: usize = self.grid.sectors[s].len();
        if length > 1 {
            for i in 0 .. length - 1 {
                for j in i + 1 .. length {
                    self.accs_pots_eval_pair_mut(
                        self.grid.sectors[s][i],
                        self.grid.sectors[s][j],
                    );
                }
            }
        }
    }

    fn update_accs_pots_grid_cross(
        &mut self,
        l: usize,
    ) {
        let (m,n) = self.grid.crosslist[l];
        let mlength: usize = self.grid.sectors[m].len();
        let nlength: usize = self.grid.sectors[n].len();
        for i in 0 .. mlength {
            for j in 0 .. nlength {
                self.accs_pots_eval_pair_mut(
                    self.grid.sectors[m][i],
                    self.grid.sectors[n][j],
                );
            }
        }
    }

    fn update_accs_pots_thread(
        &mut self,
    ) {
        self.zero_acc_pot();
        self.grid.clear_sectors();
        (0..).take(self.cparms.nmol).for_each(|i| {
            self.grid.add_index(i, &self.coords[i]);
        });

        let minparms: MinThreadParameters = MinThreadParameters {
            n1d: self.grid.n1d,
            n2d: self.grid.n2d,
            potscale: self.pparms.potscale,
            rmax2: self.pparms.rmax2,
            sigma2: self.pparms.sigma2,
            sigma6: self.pparms.sigma6,
            epsilon: self.pparms.epsilon,
            ebeta: self.pparms.ebeta,
            massinv: self.pparms.massinv,
            alpha: self.pparms.alpha,
            re: self.pparms.re,
            upre: self.pparms.upre,
            apre: self.pparms.apre,
            accmethod: self.cparms.accmethod,
            evalmethod: self.cparms.evalmethod,
            periodic: self.config.periodic,
            bound: self.pparms.bound,
        };

        let minpointers: MinThreadPointers = MinThreadPointers {
            coords: (&self.coords[0] as *const Simd<f64, LANES>) as usize,
            accs: (&self.accs[0] as *const Simd<f64, LANES>) as usize,
            pots: (&self.pots[0] as *const f64) as usize,
            sectors: (&self.grid.sectors[0] as *const Vec<usize>) as usize,
            cross: (&self.grid.crosslist[0] as *const (usize, usize)) as usize,
            layers: (&self.grid.layerbounds[0] as *const (usize, usize)) as usize,
            lines: (&self.grid.linebounds[0] as *const (usize, usize)) as usize,
            blocks: (&self.grid.blockbounds[0] as *const (usize, usize)) as usize,
        };

        match minparms.accmethod {
            2 => {
                self.acc_thread_work_layers(
                    &minparms,
                    &minpointers,
                );
            }
            3 => {
                self.acc_thread_work_lines(
                    &minparms,
                    &minpointers,
                );
            }
            4 => {
                self.acc_thread_work_blocks(
                    &minparms,
                    &minpointers,
                );
            }
            _ => {
                println!("somehow got a bad accmethod into update_acc_pots_thread");
                process::exit(1);
            }
        }
    }

    fn acc_thread_work_blocks(
        &mut self,
        minparms: &MinThreadParameters,
        minpointers: &MinThreadPointers,
    ) {
        let length: usize = self.grid.blockbounds.len();
        for off2 in 0usize ..= 1 {
            for off1 in 0usize ..= 1{
                for off0 in 0usize ..= 1{
                    let mut jobs_out: usize = 0;
                    let mut worker: usize = 0;
                    let mut l0: usize = off0;
                    let mut l1: usize = off1;
                    let mut l2: usize = off2;
                    while (l0 + minparms.n1d*l1 + minparms.n2d*l2) < length {
                        let block: usize = l0 + minparms.n1d*l1 + minparms.n2d*l2;
                        if jobs_out >= self.config.nworker {
                            worker = self.work_rslt_channel.recv()
                                .expect("trouble getting all clear back from worker state 1.");
                        }
                        
                        self.work_requ_channels[worker].send((
                            block,
                            minparms.clone(),
                            minpointers.clone(),
                        ))
                            .expect("trouble sending job to worker.");
                        if jobs_out < self.config.nworker {
                            worker += 1;
                            jobs_out += 1;
                        }
                        l0 += 2;
                        if l0 >= minparms.n1d {
                            l0 = off0;
                            l1 += 2;
                            if l1 >= minparms.n1d {
                                l1 = off1;
                                l2 += 2;
                            }
                        }
                    }
                    while jobs_out > 0 {
                        _ = self.work_rslt_channel.recv()
                            .expect("trouble getting all clear back from worker state 2.");
                        jobs_out -= 1;
                    }
                }
            }
        }
    }
    
    fn acc_thread_work_lines(
        &mut self,
        minparms: &MinThreadParameters,
        minpointers: &MinThreadPointers,
    ) {
        let length: usize = self.grid.linebounds.len();
        for off1 in 0usize ..= 1 {
            for off0 in 0usize ..= 1 {
                let mut jobs_out: usize = 0;
                let mut worker: usize = 0;
                let mut l0: usize = off0;
                let mut l1: usize = off1;
                while (l0 + minparms.n1d*l1) < length {
                    let block: usize = l0 + minparms.n1d*l1;
                    if jobs_out >= self.config.nworker {
                        worker = self.work_rslt_channel.recv()
                            .expect("trouble getting all clear back from worker state 1.");
                    }

                    self.work_requ_channels[worker].send((
                        block,
                        minparms.clone(),
                        minpointers.clone(),
                    ))
                        .expect("trouble sending job to worker.");
                    if jobs_out < self.config.nworker {
                        worker += 1;
                        jobs_out += 1;
                    }
                    l0 += 2;
                    if l0 >= minparms.n1d {
                        l0 = off0;
                        l1 += 2;
                    }
                }
                while jobs_out > 0 {
                    _ = self.work_rslt_channel.recv()
                        .expect("trouble getting all clear back from worker state 2.");
                    jobs_out -= 1;
                }
            }
        }
    }
    
    fn acc_thread_work_layers(
        &mut self,
        minparms: &MinThreadParameters,
        minpointers: &MinThreadPointers,
    ) {
        let length: usize = self.grid.layerbounds.len();
        for offset in 0usize ..= 1 {
            let mut jobs_out: usize = 0;
            let mut worker: usize = 0;
            let mut layer: usize = offset;
            //want to make this go to length-1
            while layer < length-1 {
                if jobs_out >= self.config.nworker {
                    worker = self.work_rslt_channel.recv()
                        .expect("trouble getting all clear back from worker.");
                }
                self.work_requ_channels[worker].send((layer,
                                                      minparms.clone(),
                                                      minpointers.clone()))
                    .expect("trouble sending job to worker.");
                if jobs_out < self.config.nworker {
                    worker += 1;
                    jobs_out += 1;
                }
                layer += 2;
            }
            // wait for worker 0 to be available. Then if lengths is
            // odd and offset is zero or if length is even and offset
            // is 1, throw the last row at worker 0.
            worker = self.work_rslt_channel.recv()
                .expect("trouble getting all clear back from worker.");
            jobs_out -= 1;
            while worker != 0 {
                worker = self.work_rslt_channel.recv()
                    .expect("trouble getting all clear back from worker.");
                jobs_out -= 1;
            }
            if (length.trailing_zeros() > 0 && offset == 1) ||
                (length.trailing_ones() > 0 && offset == 0) {
                    layer = length - 1;
                    self.work_requ_channels[worker].send((layer,
                                                          minparms.clone(),
                                                          minpointers.clone()))
                        .expect("trouble sending job to worker.");
                    jobs_out += 1;
                }
            if jobs_out > 0 {
                for _ in 0 .. jobs_out {
                    let _ = self.work_rslt_channel.recv()
                        .expect("trouble getting all clear back from worker.");
                }
            }
        }
    }
    
    fn accs_pots_eval_pair_mut(
        &mut self,
        i: usize,
        j: usize,
    ) {
        let (status, potential, acceleration) = self.accs_pots_eval_pair(i, j);
        if status {
            self.accs[i] -= acceleration ;
            self.accs[j] += acceleration ;
            self.pots[i] += 0.5 * potential ;
            self.pots[j] += 0.5 * potential ;
        }
    }

    fn accs_pots_eval_pair(
        &self,
        i: usize,
        j: usize,
    ) -> (bool, f64, Simd<f64, LANES>) {
        let mut dir: Simd<f64, LANES> = self.coords[i] - self.coords[j];
        if self.config.periodic > 0 {
            for k in (0..).take(3usize) {
                if dir[k].abs() > self.pparms.bound {
                    dir[k] = - dir[k].signum() * (2.0*self.pparms.bound - dir[k].abs());
                }
            }
        }
        let rad2 = (dir * dir).reduce_sum();
        if rad2 <= self.pparms.rmax2 {
            let (potential, acceleration) = match self.cparms.evalmethod {
                0 => {self.eval_acc_pot_lj_rad2(rad2)}
                1 => {self.eval_acc_pot_lj_rad6(rad2)}
                2 => {self.eval_acc_pot_morse_rad2(rad2)}
                _ => {
                    println!("Bad evaluation method");
                    process::exit(1);
                }
            };
            dir *= Simd::<f64, LANES>::splat(
                acceleration * self.pparms.potscale);
            (true, potential * self.pparms.potscale, dir)
        } else {
            (false, 0.0, dir)
        }
    }
    
    fn eval_acc_pot_lj_rad6(
        &self,
        rad2: f64,
    ) -> (f64, f64) {
        let radfact6: f64 = self.pparms.sigma6 / (rad2 * rad2 * rad2) ;
        let radfact12: f64 = radfact6 * radfact6;
        let mut potential: f64 =
            4.0 * self.pparms.epsilon * (radfact12 - radfact6);
        let mut acceleration: f64 =
            24.0 * self.pparms.epsilon * self.pparms.massinv *
            (radfact6 - 2.0 * radfact12) / rad2;
        if potential.abs() > 20.0*self.pparms.ebeta {
            potential = 20.0*self.pparms.ebeta * potential.signum() ;
            acceleration = 0.0;
        }
        (potential, acceleration)
    }

    fn eval_acc_pot_lj_rad2(
        &self,
        rad2: f64,
    ) -> (f64, f64) {
        let radfact2: f64 = self.pparms.sigma2 / rad2;
        let radfact6: f64 = radfact2 * radfact2 * radfact2;
        let radfact12: f64 = radfact6 * radfact6;
        let mut potential: f64 =
            4.0 * self.pparms.epsilon * (radfact12 - radfact6);
        let mut acceleration: f64 =
            24.0 * self.pparms.epsilon * self.pparms.massinv *
            (radfact6 - 2.0 * radfact12) / rad2;
        if potential.abs() > 20.0 * self.pparms.ebeta {
            potential = 20.0*self.pparms.ebeta * potential.signum() ;
            acceleration = 0.0;
        }
        (potential, acceleration)
    }

    fn eval_acc_pot_morse_rad1(
        &self,
        rad: f64,
    ) -> (f64, f64) {
        let radinv: f64 = 1.0 / rad;
        let expfact: f64
            = (self.pparms.alpha * (self.pparms.re - rad)).exp();
        let potential: f64
            = self.pparms.upre * (expfact * (expfact - 2.0) -1.0);
        let acceleration: f64
            = self.pparms.apre * (expfact * (1.0 - expfact)) * radinv;
        (potential, acceleration)
    }

    fn eval_acc_pot_morse_rad2(
        &self,
        rad2: f64,
    ) -> (f64, f64) {
        self.eval_acc_pot_morse_rad1(rad2.sqrt())
    }

    fn update_vels_coords(
        &mut self,
    ) {
        self.tparms.kinetic = 0.0;
        self.tparms.potential = 0.0;
        let mut best_centre: usize = 0;
        let mut distance: f64 = 0.0;
        (0..).take(3usize).for_each(|k| {
            distance += self.coords[0][k].abs();
        });
        for i in (0..).take(self.cparms.nmol) {
            let mut speed: f64 = 0.0;
            let mut ldist: f64 = 0.0;
            self.vels[i] += self.pparms.dtsplat * self.accs[i];
            self.coords[i] += self.pparms.dtsplat * self.vels[i];
            speed += (self.vels[i] * self.vels[i]).reduce_sum();
            self.tparms.kinetic += (self.vels[i] * self.vels[i]).reduce_sum();
            for k in (0..).take(3usize) {
                ldist += self.coords[i][k].abs();
            }
            self.tparms.potential += self.pots[i];
            self.speeds.insert_value(speed.sqrt());
            if ldist < distance {
                distance = ldist;
                best_centre = i;
            }            
        }
        for i in (0..).take(self.cparms.nmol) {
            if i != best_centre {
                let tmp: Simd<f64, LANES> = self.coords[i] - self.coords[best_centre];
                let sep: f64 = (tmp * tmp).reduce_sum();
                if sep < self.pparms.bound2 {
                    self.radialdist.insert_value(sep.sqrt());
                }                    
            }
        }
        self.tparms.kinetic *= 0.5 * self.pparms.mass /
            (BOLTZMANN * self.pparms.temperature * self.cparms.nmol as f64);
        self.tparms.potential *= 1.0 /
            (BOLTZMANN * self.pparms.temperature * self.cparms.nmol as f64);
    }

    fn bounce_check(
        &mut self,
    ) {
        match self.config.periodic {
            2 => {
                self.bounce_check_periodic();
            }
            _ => {
                self.bounce_check_nonperiodic();
            }
        }
    }

    fn bounce_check_periodic(
        &mut self
    ) {
        let mut localimpulse: f64 = 0.0;
        let mut rng = thread_rng();
        let prob: f64 = 1.0/(self.cparms.nmol * self.cparms.nmol) as f64;
        
        for i in (0..).take(self.cparms.nmol) {
            for k in (0..).take(3usize) {
                if self.coords[i][k].abs() > self.pparms.bound {
                    let posneg: f64 = self.coords[i][k].signum();
                    self.coords[i][k] = - posneg *
                        (2.0 * self.pparms.bound - self.coords[i][k].abs());
                    localimpulse += 2.0*self.vels[i][k].abs();
                    self.single_kinetic_bounce(i, k);
                }
                //move out of main check loop to avoid adding bias at
                // the bounds of the box ... like there is a blackbody
                // spectrum randomly interacting with the
                // particles. For some reason, pulling out of the wall
                // loop adds a crazy bias if it is called too often.
                
                if rng.gen_range(0.0..1.0) < prob {
                    self.single_kinetic_bounce(i, k);
                }
            }
        }
        
        self.tparms.impulse += localimpulse * self.pparms.mass /
            (24.0 * self.pparms.bound * self.pparms.bound * self.pparms.dt *
             self.cparms.samplerate as f64) ;
    }

    fn bounce_check_nonperiodic(
        &mut self
    ) {
        let mut localimpulse: f64 = 0.0;
        let mut rng = thread_rng();
        let prob: f64 = 1.0/(self.cparms.nmol * self.cparms.nmol) as f64;
        
        for i in (0..).take(self.cparms.nmol) {
            for k in (0..).take(3usize) {
                if self.coords[i][k].abs() > self.pparms.bound {
                    let posneg: f64 = self.coords[i][k].signum();
                    self.coords[i][k] = posneg * self.pparms.adjust;
                    self.vels[i][k] *= -1.0;
                    localimpulse += 2.0*self.vels[i][k].abs();
                    self.single_kinetic_bounce(i, k);
                }
                //move out of main check loop to avoid adding bias at
                // the bounds of the box ... like there is a blackbody
                // spectrum randomly interacting with the
                // particles. For some reason, pulling out of the wall
                // loop adds a crazy bias if it is called too often.

                if rng.gen_range(0.0..1.0) < prob {
                    self.single_kinetic_bounce(i, k);
                }
            }
        }

        self.tparms.impulse += localimpulse * self.pparms.mass /
            (24.0 * self.pparms.bound * self.pparms.bound * self.pparms.dt *
             self.cparms.samplerate as f64) ;
    }

    fn single_kinetic_bounce(
        &mut self,
        j: usize,
        k: usize,
    ) {
        let mut rng = thread_rng();
        let mut localdeltaq: f64 = 0.0;

        //This has been changed a lot. It is probably the physically
        // most complex aspect of the code. I tried several
        // strategies, but all of the other strategies are slower, and
        // introduce bias that causes drift away from thermal
        // equilibrium. Sometimes catastrophically so. Go ahead and
        // tinker with this if you want, but you're not likely to get
        // something better than this.

        let vref: f64 = self.vels[j][k];
        let eref: f64 = 0.5 * self.pparms.mass * vref * vref;
        let vtmp: f64 = rng.gen_range(0.0 .. self.pparms.vmax)*vref.signum();
        let etmp: f64 = 0.5 * self.pparms.mass * vref * vref;
        if etmp <= eref
            || ((eref - etmp) * self.pparms.beta).exp() >
            rng.gen_range(0.0 .. 1.0) {
                self.vels[j][k] = vtmp ;
            }
        
        localdeltaq += self.vels[j][k]*self.vels[j][k] - vref*vref;

        self.tparms.deltaq += 0.5 * localdeltaq * self.pparms.mass /
            self.pparms.dt;
    }

    fn pre_accel(
        &mut self,
    ) {
        let lbound = self.pparms.bound;
        self.pparms.update_bound();
        self.bounce_check();
        if self.pparms.bound != lbound {
            self.grid.rescale_grid(
                self.pparms.bound,
                self.cparms.closescale,
                self.pparms.sigma,
            );
            self.bounce_check();
            self.setup_workers();
            if self.config.voltherm > 0 {
                self.autotherm();
            }
        }
        if (self.config.voltherm > 0 &&
            self.config.voltherm < 3 &&
            self.tparms.itcount < 0) ||
            (3 == self.config.voltherm &&
             self.tparms.itcount > -self.cparms.anneal/2) {
                self.autotherm();
            }
    }

    fn post_accel(
        &mut self,
    ) {
        self.update_vels_coords();

        self.tparms.itcount += 1;

        if self.tparms.itcount == -self.cparms.anneal/2 {
            self.speeds.set_dx_from_max();
            self.radialdist.set_dx_from_max();
            if self.config.voltherm > 1 &&
                (1.5 - self.averages.kinetic1mean).abs() > 0.1 {
                    self.tparms.itcount = - self.cparms.anneal;
                }
        }
        
        if self.tparms.itcount == 0 {
            self.tparms.elapsedtime = 0.0;
            self.speeds.set_dx_from_max();
            self.radialdist.set_dx_from_max();
        }

        if self.tparms.itcount > 0 {
            self.tparms.elapsedtime =
                (self.tparms.itcount as f64) * self.pparms.dt;
        }
        
        if self.tparms.itcount % self.cparms.samplerate as isize == 0 {
            self.tparms.update_extrinsics() ;
            self.averages.update_averages(
                AverageRawData{
                    itcount: self.tparms.itcount,
                    kinetic: self.tparms.kinetic,
                    potential: self.tparms.potential,
                    pressure: self.tparms.pressure,
                    deltaqdeltat: self.tparms.deltaqdeltat,
                    itspersec: self.tparms.itspersec,
                });
        }
    }

    fn scale_energies(
        &mut self,
        factor: f64,
    ) {
        let mut localdeltaq: f64 = 0.0;
        let factorsplat: Simd<f64, LANES> = Simd::<f64, LANES>::splat(factor);
        for i in (0..).take(self.cparms.nmol) {
            let v2: f64 = (self.vels[i] * self.vels[i]).reduce_sum();
            self.vels[i] *= factorsplat;
            localdeltaq += (self.vels[i] * self.vels[i]).reduce_sum() - v2;
        }
        self.tparms.deltaq += 0.5 * localdeltaq * self.pparms.mass /
            self.pparms.dt;
    }
    
    fn check_channels(
        &mut self,
    ) {
        self.check_coords_channel();
        self.check_text_channel();
        self.check_dist_channel();
        self.check_dist_channel();
    }

    fn check_coords_channel(
        &mut self,
    ) {
        let scale: f64 = (2.0f64.ln() / 32.0).exp();
        let scale_1: f64 = 1.0/scale;
        match self.coords_requ_channel.try_recv() {
            Err(TryRecvError::Empty) => {}
            Ok(1) => {
                //send coordinates down the channel
                let coordsusize: usize =
                    (&self.coords[0] as *const Simd<f64, LANES>) as usize;
                let _ = self
                    .coords_send_channel
                    .send(
                        TransferCoords{
                            nmol: self.cparms.nmol,
                            coords: coordsusize,
                        });
            }
            Ok(2) => {
                //energy down
                self.scale_energies(scale_1);
                self.tparms.itcount = -self.cparms.anneal;
            }
            Ok(3) => {
                //energy up
                self.scale_energies(scale);
                self.tparms.itcount = -self.cparms.anneal;
            }
            Ok(4) => {
                //temperature up
                let temperature = self.pparms.temperature;
                self.pparms.reset_temperature(temperature * scale);
                self.setup_workers();
                self.tparms.itcount = -self.cparms.anneal;
            }
            Ok(5) => {
                //temperature down
                let temperature = self.pparms.temperature;
                self.pparms.reset_temperature(temperature * scale_1);
                self.setup_workers();
                self.tparms.itcount = -self.cparms.anneal;
            }
            Ok(6) => {
                //dt steps bigger (finer grain)
                self.pparms.reset_dtsteps(self.pparms.dtsteps * 2);
                self.tparms.itcount = -self.cparms.anneal;
            }
            Ok(7) => {
                //dt steps smaller (coarser grain)
                self.pparms.reset_dtsteps(self.pparms.dtsteps / 2);
                self.tparms.itcount = -self.cparms.anneal;
            }
            Ok(8) => {
                //box bigger
                self.pparms.reset_boundgoal(self.pparms.boundgoal * scale);
                self.setup_workers();
                self.tparms.itcount = -self.cparms.anneal;
            }
            Ok(9) => {
                //box smaller
                self.pparms.reset_boundgoal(self.pparms.boundgoal * scale_1);
                self.setup_workers();
                self.tparms.itcount = -self.cparms.anneal;
            }
            Ok(10) => {
                //toggle potential on/off
                if self.pparms.potscale > f64::EPSILON {
                    self.pparms.potscale = 0.0;
                } else {
                    self.pparms.potscale = 1.0;
                }
                self.setup_workers();
                self.tparms.itcount = -self.cparms.anneal;
            }
            Ok(11) => {
                //autotherm
                self.autotherm();
                self.tparms.itcount = -self.cparms.anneal;
            }
            Ok(12) => {
                //toggle info display
                self.tparms.displayactive = !self.tparms.displayactive;
            }
            _ => {
                println!("sent a bad signal down the coords request channel");
                process::exit(1);
            }
        }
    }
    
    fn check_text_channel(
        &mut self,
    ) {
        match self.text_requ_channel.try_recv() {
            Err(TryRecvError::Empty) => {}
            Ok(1) => {
                self.send_text_to_eyecandy();
            }
            _ => {
                println!("sent a bad signal down the text request channel");
                process::exit(1);
            }
        }
    }

    fn check_dist_channel(
        &mut self,
    ) {
        match self.dist_requ_channel.try_recv() {
            Err(TryRecvError::Empty) => {}
            Ok(1) => {
                self.speeds.normalize((0, 1.0));
                self.dist_send_channel.send(
                    TransferDist{
                        nval: self.speeds.nval,
                        vals: self.speeds.normalized.clone(),
                    }
                ).expect("Trouble sending speed data");   
            }
            Ok(2) => {
                let normal: f64 = 
                    self.pparms.bound * self.pparms.bound * self.pparms.bound
                    / 3.0;
                self.radialdist.normalize((1,normal));
                self.dist_send_channel.send(
                    TransferDist{
                        nval: self.radialdist.nval,
                        vals: self.radialdist.normalized.clone(),
                    }
                ).expect("Trouble sending radial distribution data");
            }
            _ => {
                println!("sent a bad signal down the dist request channel");
                process::exit(1);
            }
        }
    }

    fn send_text_to_eyecandy(
        &self,
    ) {
        let aves: AverageResult = self.averages.average_result();
        
        let data: MacroStateData = MacroStateData {
            nmol: self.cparms.nmol,
            itmax: self.cparms.itmax,
            itcount: self.tparms.itcount,
            ittotal: self.tparms.ittotal,
            elapsedtime: self.tparms.elapsedtime,
            timetotal: self.tparms.timetotal,
            samples: self.cparms.samplerate,
            element: "Ne-20\nTc = 44.448 K\nPc = 2.76 MPa\nVc = 41.7 mL/mol".to_string(),
            temperature: self.pparms.temperature,
            bound: self.pparms.bound,
            volume: self.pparms.bound*self.pparms.bound*self.pparms.bound*8.0,
            dt: self.pparms.dt,
            kinetic: self.tparms.kinetic,
            kinetic1mean: aves.kinetic1mean,
            kinetic2mean: aves.kinetic2mean,
            potential: self.tparms.potential,
            potential1mean: aves.potential1mean,
            potential2mean: aves.potential2mean,
            pressure: self.tparms.pressure,
            pressure1mean: aves.pressure1mean,
            pressure2mean: aves.pressure2mean,
            deltaqdeltat: self.tparms.deltaqdeltat,
            deltaqdeltat1mean: aves.deltaqdeltat1mean,
            deltaqdeltat2mean: aves.deltaqdeltat2mean,
            itspersec: self.tparms.itspersec,
            itspersec1mean: aves.itspersec1mean,
            itspersec2mean: aves.itspersec2mean,
            ntrackinv: aves.ntrackinv,
            displayactive: self.tparms.displayactive,
        };
        self.text_send_channel.send(
            data
        ).expect("Can't send text data to eyecandy");
    }

    pub fn setup_workers(
        &mut self,
    ) {
        let mut work_requ_channels: Vec<(ThreadSender, ThreadReceiver)> = vec![];
        let work_rslt_channel: (Sender<usize>,
                                Receiver<usize>) = unbounded();

        self.terminate_workers();
        
        self.work_rslt_channel = work_rslt_channel.1.clone();
        
        for i in (0..).take(self.config.nworker) {
            work_requ_channels.push(unbounded());
            self.work_requ_channels.push(work_requ_channels[i].0.clone());
        }
        
        for i in (0..).take(self.config.nworker) {
            let tmp_rslt_channel = work_rslt_channel.0.clone();
            let tmp_requ_channel = work_requ_channels[i].1.clone();
            let tmp_thread = thread::spawn(move || {
                acc_worker(
                    ACCWorkerConf {
                        index: i,
                        channelin: tmp_requ_channel,
                        channelout: tmp_rslt_channel,
                    });
            });
            self.work_handles.push(tmp_thread);
        }
        
    }

    pub fn terminate_workers(
        &mut self,
    ) {
        
        if !self.work_handles.is_empty() {
            for i in (0..).take(self.config.nworker) {
                self.work_requ_channels[i].send((self.cparms.nmol+1,
                                                MinThreadParameters::default(),
                                                MinThreadPointers::default()))
                    .expect("cannot tell worker to stop.");                
            }
            for _ in (0..).take(self.config.nworker) {
                self.work_handles.pop()
                    .expect("cannot pop a handle that should be popped.")
                    .join()
                    .expect("cannot join a handle for a worker that should have stopped.");
            }
            self.work_handles = vec![];
        }

        if !self.work_requ_channels.is_empty() {
            self.work_requ_channels = vec![];
        }

    }

    fn set_vars_to_current(
        &mut self,
    ) {
        let nvar: usize = self.config.variables.len();
        for i in (0..).take(nvar) {
            let (name, val) = self.config.variables[i]
                .get_current_value();
            self.parse_var_change(
                name,
                val,
            );
        }
    }
    
    fn parse_var_change(
        &mut self,
        name: String,
        val: f64,
    ) {
        match name.as_str() {
            "sigma" => {
                self.pparms.sigma = val;
                let temperature = self.pparms.temperature;
                self.pparms.reset_temperature(temperature);
                self.pparms.setup_morse();
                self.setup_workers();
                self.tparms.itcount = -self.cparms.anneal;
            }
            "epsilon" => {
                self.pparms.epsilon = val;
                let temperature = self.pparms.temperature;
                self.pparms.reset_temperature(temperature);
                self.pparms.setup_morse();
                self.setup_workers();
                self.tparms.itcount = -self.cparms.anneal;
            }
            "mass" => {
                self.pparms.mass = val;
                self.pparms.massinv = 1.0/val;
                let temperature = self.pparms.temperature;
                self.pparms.reset_temperature(temperature);
                self.pparms.setup_morse();
                self.setup_workers();
                self.tparms.itcount = -self.cparms.anneal;
            }
            "temperature" => {
                self.pparms.reset_temperature(val);
                self.setup_workers();
                self.tparms.itcount = -self.cparms.anneal;
            }
            "molarvolume" => {
                let tmpbound: f64 =
                    self.pparms.bound_from_volume(
                        val,
                        self.cparms.nmol,
                    );
                self.pparms.reset_boundgoal(tmpbound);
                self.setup_workers();
                self.tparms.itcount = -self.cparms.anneal;
            }
            "molardensity" => {
                let tmpbound: f64 =
                    self.pparms.bound_from_density(
                        val,
                        self.cparms.nmol,
                    );
                self.pparms.reset_boundgoal(tmpbound);
                self.setup_workers();
                self.tparms.itcount = -self.cparms.anneal;
            }
            _ => {
                println!("variable not yet supported: {name}");
                process::exit(1);
            }
        }
    }
    
    pub fn run_iterations(
        &mut self,
    ) {
        let mut timer = Instant::now();

        let nvars: usize = self.config.variables.len();

        self.set_vars_to_current();

        let mut complete: bool = false;

        if 2 == self.cparms.accmethod {
            self.setup_workers();
        }

        let mut casenumber: usize = 0;

        while !complete {

            let mut flipped: bool = true;
            let mut count: usize = 0;
            
            let filename: String = format!("{}.{:04}.state",
                                           self.config.basename,
                                           casenumber,
            );

            if !Path::new(&filename.as_str()).is_file()
            {
                
                loop {
                    self.pre_accel();
                    self.update_accs_pots();
                    self.post_accel();
                    self.check_channels();
                    
                    if self.cparms.itmax > 0 &&
                        self.tparms.itcount >= self.cparms.itmax as isize {
                            break;
                        }
                    
                    if self.tparms.itcount % self.cparms.samplerate as isize == 0 {
                        let elapsed: u128 = timer.elapsed().as_nanos();
                        
                        self.tparms.itspersec = (self.cparms.samplerate as f64) * 1.0e9 /
                            (elapsed as f64);
                        
                        if self.tparms.itcount <= 0 && self.config.voltherm > 0 {
                            self.autotherm();
                        }
                        
                        if self.tparms.itcount >= 0 &&
                            self.tparms.itcount % (self.cparms.samplerate *
                                                   self.averages.ntrack) as isize == 0 &&
                            self.config.logaves == 1 {
                                let volume: f64 = 8.0 * self.pparms.bound *
                                    self.pparms.bound * self.pparms.bound *
                                    AVOGADRO * 1.0e6 / self.cparms.nmol as f64;
                                let filename: String = format!("{}.{:04}.log",
                                                               self.config.basename,
                                                               casenumber,
                                );
                                self.averages.dump(
                                    self.pparms.temperature,
                                    volume,
                                    filename.clone(),
                                );
                            }
                        
                        timer = Instant::now();
                    }
                }
                
                let volume: f64 = 8.0 * self.pparms.bound *
                    self.pparms.bound * self.pparms.bound *
                    AVOGADRO * 1.0e6 / self.cparms.nmol as f64;
                let filename: String = format!("{}.{:04}.state",
                                               self.config.basename,
                                               casenumber,
                );
                let mut tmpconfig: SimConfig = self.config.clone();
                tmpconfig.molarvolume = volume;
                tmpconfig.temperature = self.pparms.temperature;
                tmpconfig.molardensity = 1.0 / volume;
                tmpconfig.dump_state(
                    filename,
                    self.coords.clone(),
                    self.vels.clone(),
                );
                let filename: String = format!("{}.{:04}.speeds",
                                               self.config.basename,
                                               casenumber,
                );
                self.speeds.dump_distribution(filename);
                let filename: String = format!("{}.{:04}.radialdist",
                                               self.config.basename,
                                               casenumber,
                );
                self.radialdist.dump_distribution(filename);
            } else if !(is_file_empty(filename.as_str())
                        .expect("problem checking for empty file")) {
                let mut tmpconfig: SimConfig = SimConfig::default();
                tmpconfig.restore_phase(filename);
                if self.config.molfact != tmpconfig.molfact {
                    println!("missmatch number of molecules when reading in phase from case {casenumber}");
                    println!("main config has molfact = {}", self.config.molfact);
                    println!("phase config has molfact = {}", self.config.molfact);
                    process::exit(1);
                }
                for i in (0..).take(self.cparms.nmol) {
                    self.coords[i] = tmpconfig.coords[i];
                    self.vels[i] = tmpconfig.vels[i];
                }
            }

            while flipped && count < nvars {
                let (name, val, flipper) = self.config.variables[count].get_next_value();
                self.parse_var_change(name, val);
                count += 1;
                flipped = flipper;
            }

            if flipped && count == nvars {
                complete = true;
            }
            casenumber += 1;
        }

        if 2 == self.cparms.accmethod {
            self.terminate_workers();
        }
        
        println!("mean itspersec: {}", self.averages.itspersec1mean);
        process::exit(0);
    }

    pub fn autotherm(
        &mut self
    ) {
        if self.tparms.kinetic > 0.0 {
            let ratio: f64 = (1.5f64 / self.tparms.kinetic).sqrt();
            let ratiosplat: Simd<f64, LANES> = Simd::<f64, LANES>::splat(ratio);
            let mut localdeltaq: f64 = 0.0 ;
            (0..).take(self.cparms.nmol).for_each(|i| {
                let v2old: f64 = (self.vels[i] * self.vels[i]).reduce_sum();
                self.vels[i] *= ratiosplat;
                let v2new: f64 = (self.vels[i] * self.vels[i]).reduce_sum();
                localdeltaq += v2new - v2old;
            });
            self.tparms.deltaq += 0.5 * localdeltaq * self.pparms.mass /
                self.pparms.dt;
        }
    }
}

//Helpers for the threaded acceleration

struct ACCWorkerConf{
    index: usize,
    channelin: Receiver<(usize,
                         MinThreadParameters,
                         MinThreadPointers)>,
    channelout: Sender<usize>,
}

#[derive(Clone)]
struct MinThreadParameters{
    n1d: usize,
    n2d: usize,
    potscale: f64,
    rmax2: f64,
    sigma2: f64,
    sigma6: f64,
    epsilon: f64,
    ebeta: f64,
    massinv: f64,
    alpha: f64,
    re: f64,
    upre: f64,
    apre: f64,
    accmethod: usize,
    evalmethod: usize,
    periodic: usize,
    bound: f64,
}

impl Default for MinThreadParameters {
    fn default() -> Self {
        MinThreadParameters {
            n1d: 0,
            n2d: 0,
            potscale: 0.0,
            rmax2: 0.0,
            sigma2: 0.0,
            sigma6: 0.0,
            epsilon: 0.0,
            ebeta: 0.0,
            massinv: 0.0,
            alpha: 0.0,
            re: 0.0,
            upre: 0.0,
            apre: 0.0,
            accmethod: 0,
            evalmethod: 0,
            periodic: 0,
            bound: 1.0,
        }
    }
}

#[derive(Clone,Default)]
struct MinThreadPointers{
    coords: usize,
    accs: usize,
    pots: usize,
    sectors: usize,
    cross: usize,
    layers: usize,
    lines: usize,
    blocks: usize,
}

fn acc_worker(
    conf: ACCWorkerConf,
) {
    let index: usize = conf.index;
    let channelin: Receiver<(usize,
                             MinThreadParameters,
                             MinThreadPointers)> = conf.channelin;
    let channelout: Sender<usize> = conf.channelout;

    loop {
        //channelin should get a single layer to do.
        match channelin.recv() {
            Ok(job) => {
                if job.0 > (job.1.n1d * job.1.n2d) {
                    break;
                }
                let parameters: MinThreadParameters = job.1.clone();
                let pointers: MinThreadPointers = job.2.clone();
                match parameters.accmethod {
                    2 => {
                        work_layer(
                            &parameters,
                            &pointers,
                            job.0,
                        );
                    }
                    3 => {
                        work_line(
                            &parameters,
                            &pointers,
                            job.0,
                        );
                    }
                    4 => {
                        work_block(
                            &parameters,
                            &pointers,
                            job.0,
                        );
                    }
                    _ => {
                        println!("bad accmethod for threaded cases");
                        process::exit(1);
                    }
                }
                channelout.send(index)
                    .expect("can't send back job complete signal");
            }
            _ => {
                println!("something broke the channel for worker {index}");
                process::exit(1);
            }
        }
    }
}

fn work_layer(
    parms: &MinThreadParameters,
    pointers: &MinThreadPointers,
    layer: usize,
) {
    let rawcross: *const (usize, usize) = pointers.cross as *const (usize, usize);
    let rawlayers: *const (usize, usize) = pointers.layers as *const (usize, usize);

    let start: usize;
    let stop: usize;

    unsafe {
        start = (*rawlayers.add(layer)).0;
        stop = (*rawlayers.add(layer)).1;
    }

    for i in start .. stop {
        let pair: (usize, usize);
        unsafe {
            pair = *rawcross.add(i);
        }
        if pair.0 == pair.1 {
            work_accs_pots_thread_self(
                parms,
                pointers,
                pair,
            );
        } else {
            work_accs_pots_thread_cross(
                parms,
                pointers,
                pair,
            );
        }
    }
}

fn work_line(
    parms: &MinThreadParameters,
    pointers: &MinThreadPointers,
    line: usize,
) {
    let rawcross: *const (usize, usize) = pointers.cross as *const (usize, usize);
    let rawlines: *const (usize, usize) = pointers.lines as *const (usize, usize);

    let start: usize;
    let stop: usize;

    unsafe {
        start = (*rawlines.add(line)).0;
        stop = (*rawlines.add(line)).1;
    }

    for i in start .. stop {
        let pair: (usize, usize);
        unsafe {
            pair = *rawcross.add(i);
        }
        if pair.0 == pair.1 {
            work_accs_pots_thread_self(
                parms,
                pointers,
                pair,
            ) ;
        } else {
            work_accs_pots_thread_cross(
                parms,
                pointers,
                pair,
            );
        }
    }
}

fn work_block(
    parms: &MinThreadParameters,
    pointers: &MinThreadPointers,
    block: usize,
) {
    let rawcross: *const (usize, usize) = pointers.cross as *const (usize, usize);
    let rawblocks: *const (usize, usize) = pointers.blocks as *const (usize, usize);

    let start: usize;
    let stop: usize;

    unsafe {
        start = (*rawblocks.add(block)).0;
        stop = (*rawblocks.add(block)).1;
    }

    for i in start .. stop {
        let pair: (usize, usize);
        unsafe {
            pair = *rawcross.add(i);
        }
        if pair.0 == pair.1 {
            work_accs_pots_thread_self(
                parms,
                pointers,
                pair,
            );
        } else {
            work_accs_pots_thread_cross(
                parms,
                pointers,
                pair,
            );
        }
    }
}

fn work_accs_pots_thread_self(
    parms: &MinThreadParameters,
    pointers: &MinThreadPointers,
    pair: (usize, usize),
) {
    let s: usize = pair.0;
    let rawsectors: *const Vec<usize> = pointers.sectors as *const Vec<usize>;
    let length: usize;
    unsafe {
        length = (*rawsectors.add(s)).len();
    }

    if length > 1 {
        for i in 0 .. length - 1 {
            for j in i + 1 .. length {
                let ii: usize;
                let jj: usize;
                unsafe {
                    ii = (&(*rawsectors.add(s)))[i];
                    jj = (&(*rawsectors.add(s)))[j];
                }
                work_accs_pots_eval_pair_mut(
                    parms,
                    pointers,
                    (ii, jj),
                );
            }
        }
    }
}

fn work_accs_pots_thread_cross(
    parms: &MinThreadParameters,
    pointers: &MinThreadPointers,
    pair: (usize, usize),
) {
    let m: usize = pair.0;
    let n: usize = pair.1;
    let rawsectors: *const Vec<usize> = pointers.sectors as *const Vec<usize>;
    let mlength: usize;
    let nlength: usize;
    unsafe {
        mlength = (*rawsectors.add(m)).len();
        nlength = (*rawsectors.add(n)).len();
    }

    for i in 0 .. mlength {
        for j in 0 .. nlength {
            let ii: usize;
            let jj: usize;
            unsafe {
                ii = (&(*rawsectors.add(m)))[i];
                jj = (&(*rawsectors.add(n)))[j];
            }
            work_accs_pots_eval_pair_mut(
                parms,
                pointers,
                (ii, jj),
            );
        }
    }
}

fn work_accs_pots_eval_pair_mut(
    parms: &MinThreadParameters,
    pointers: &MinThreadPointers,
    pair: (usize, usize),
) {
    let (status, potential, acceleration) =
        work_accs_pots_eval_pair(
            parms,
            pointers,
            pair,
        );
    if status {
        let rawaccs: *mut Simd<f64, LANES> = pointers.accs as *mut Simd<f64, LANES>;
        let rawpots: *mut f64 = pointers.pots as *mut f64;
        // problem is in the acceleration part of this unsafe block.
        unsafe {
            (*rawaccs.add(pair.0)) -= acceleration;
            (*rawaccs.add(pair.1)) += acceleration;
            (*rawpots.add(pair.0)) += 0.5 * potential;
            (*rawpots.add(pair.1)) += 0.5 * potential;
        }
    }
}

fn work_accs_pots_eval_pair(
    parms: &MinThreadParameters,
    pointers: &MinThreadPointers,
    pair: (usize, usize),
) -> (bool, f64, Simd<f64, LANES>) {
    let rawcoords: *const Simd<f64, LANES> = pointers.coords as *const Simd<f64, LANES>;
    let mut dir: Simd<f64, LANES>;
    unsafe {
        dir = (*rawcoords.add(pair.0)) - (*rawcoords.add(pair.1));
    }
    if parms.periodic > 0 {
        for k in (0..).take(3usize) {
            if dir[k].abs() > parms.bound {
                dir[k] = - dir[k].signum() * (2.0*parms.bound - dir[k].abs());
            }
        }
    }
    let rad2 = (dir * dir).reduce_sum();
    if rad2 <= parms.rmax2 {
        let (potential, acceleration) = match parms.evalmethod {
            0 => {work_eval_acc_pot_lj_rad2(parms, rad2)}
            1 => {work_eval_acc_pot_lj_rad6(parms, rad2)}
            2 => {work_eval_acc_pot_morse_rad2(parms, rad2)}
            _ => {
                println!("Bad evaluation method");
                process::exit(1);
            }
        };
        dir *= Simd::<f64, LANES>::splat(
            acceleration * parms.potscale
        );
        (true, potential * parms.potscale, dir)
    } else {
        (false, 0.0, dir)
    }
}

fn work_eval_acc_pot_lj_rad6(
    parms: &MinThreadParameters,
    rad2: f64,
) -> (f64, f64) {
    let radfact6: f64 = parms.sigma6 / (rad2 * rad2 * rad2);
    let radfact12: f64 = radfact6 * radfact6;
    let mut potential: f64 =
        4.0 * parms.epsilon * (radfact12 - radfact6);
    let mut acceleration: f64 =
        24.0 * parms.epsilon * parms.massinv *
        (radfact6 - 2.0 * radfact12) / rad2;
    if potential.abs() > 20.0 * parms.ebeta {
        potential = 20.0 * parms.ebeta * potential.signum();
        acceleration = 0.0;
    }
    (potential, acceleration)
}

fn work_eval_acc_pot_lj_rad2(
    parms: &MinThreadParameters,
    rad2: f64,
) -> (f64, f64) {
    let radfact2: f64 = parms.sigma2 / rad2;
    let radfact6: f64 = radfact2 * radfact2 * radfact2;
    let radfact12: f64 = radfact6 * radfact6;
    let mut potential: f64 =
        4.0 * parms.epsilon * (radfact12 - radfact6);
    let mut acceleration: f64 =
        24.0 * parms.epsilon * parms.massinv *
        (radfact6 - 2.0 * radfact12) / rad2;
    if potential.abs() > 20.0 * parms.ebeta {
        potential = 20.0 * parms.ebeta * potential.signum();
        acceleration = 0.0;
    }
    (potential, acceleration)
}

fn work_eval_acc_pot_morse_rad2(
    parms: &MinThreadParameters,
    rad2: f64,
) -> (f64, f64) {
    work_eval_acc_pot_morse_rad1(parms, rad2.sqrt())
}

fn work_eval_acc_pot_morse_rad1(
    parms: &MinThreadParameters,
    rad: f64,
) -> (f64, f64) {
    let expfact: f64
        = (parms.alpha * (parms.re - rad)).exp();
    let potential: f64
        = parms.upre * (expfact * (expfact - 2.0) - 1.0);
    let acceleration: f64
        = parms.apre * (expfact * (1.0 - expfact)) / rad;
    (potential, acceleration)
}

fn is_file_empty(
    file_path: &str
) -> Result<bool, std::io::Error> {
    let path = Path::new(file_path);

    // Get metadata for the file
    match fs::metadata(path) {
        Ok(metadata) => {
            // Check if the file's length is 0
            Ok(metadata.len() == 0)
        }
        Err(e) => {
            // Propagate any errors encountered during metadata retrieval
            Err(e)
        }
    }
}
