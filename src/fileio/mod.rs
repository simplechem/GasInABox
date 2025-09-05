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

use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, Write};
use std::path::Path;
use std::process;
use std::string::String;
use std::simd::Simd;

#[derive(Clone)]
pub struct SimConfig{
    pub basename: String,
    pub element: String,
    pub sigma: f64,
    pub epsilon: f64,
    pub mass: f64,
    pub tcrit: f64,
    pub temperature: f64,
    pub pcrit: f64,
    pub vcrit: f64,
    pub dcrit: f64,
    pub molarvolume: f64,
    pub molardensity: f64,
    pub closescale: f64,
    pub velscale: f64,
    pub molfact: usize,
    pub itermax: usize,
    pub anneal: isize,
    pub sample: usize,
    pub ntrack: usize,
    pub ndist: usize,
    pub dtsteps: usize,
    pub evalmethod: usize,
    pub accmethod: usize,
    pub voltherm: usize,
    pub logaves: usize,
    pub periodic: usize,
    pub nworker: usize,
    pub variables: Vec<VarConfig>,
    pub coords: Vec<Simd<f64, LANES>>,
    pub vels: Vec<Simd<f64, LANES>>,
    pub readin: bool,
    pub displayactive: bool,
}

impl Default for SimConfig{
    fn default() -> Self {
        SimConfig {
            basename: BASENAME.to_string(),
            element: ELEMENT.to_string(),
            sigma: SIGMA,
            epsilon: EPSILONLJ,
            mass: MASS,
            temperature: TEMPERATURE,
            tcrit: TCRIT,
            pcrit: PCRIT,
            vcrit: VCRITA,
            dcrit: DCRITA,
            molarvolume: VCRIT,
            molardensity: DCRIT,
            closescale: CLOSESCALE,
            velscale: VELSCALE,
            molfact: MOLFACT,
            itermax: ITERMAX,
            anneal: ANNEAL,
            sample: SAMPLE,
            ntrack: NTRACK,
            ndist: NDIST,
            dtsteps: DTSTEPS,
            evalmethod: EVALMETHOD,
            accmethod: ACCMETHOD,
            voltherm: 0,
            logaves: 1,
            periodic: 1,
            nworker: NWORKER,
            variables: vec![],
            coords: vec![],
            vels: vec![],
            readin: false,
            displayactive: true,
        }
    }
}

impl SimConfig {
    pub fn get_settings_from_file(
        &mut self,
        infilename: String,
    ) {
        let mut linecount: usize = 0;

        if let Ok(lines) = SimConfig::read_lines(infilename) {
            for parsable in lines.map_while(Result::ok) {
                if parsable.contains("var:") {
                    let mut tmp: VarConfig = VarConfig::default();
                    tmp.parse_input_line(&parsable);
                    self.variables.push(tmp.clone());
                } else {
                    let cleaned = SimConfig::remove_whitespace_and_comments(&parsable);
                    let parsed = SimConfig::trim_line(cleaned);
                    if parsed.1.is_empty() && !parsed.0.is_empty() {
                        println!("line {linecount}: not a key=value entry");
                        process::exit(1);
                    }
                    match &parsed.0[..] {
                        "basename" => {
                            self.basename = parsed.1.parse::<String>()
                                .expect("cannot parse basename");
                        }
                        "element" => {
                            self.element = parsed.1.parse::<String>()
                                .expect("cannot parse element");
                        }
                        "sigma" => {
                            self.sigma = parsed.1.parse::<f64>()
                                .expect("cannot parse sigma");
                        }
                        "epsilon" => {
                            self.epsilon = parsed.1.parse::<f64>()
                                .expect("cannot parse epsilon");
                        }
                        "mass" => {
                            self.mass = parsed.1.parse::<f64>()
                                .expect("cannot parse mass");
                        }
                        "temperature" => {
                            self.temperature = parsed.1.parse::<f64>()
                                .expect("cannot parse temperature");
                        }
                        "tcrit" => {
                            self.tcrit = parsed.1.parse::<f64>()
                                .expect("cannot parse tcrit");
                        }
                        "pcrit" => {
                            self.pcrit = parsed.1.parse::<f64>()
                                .expect("cannot parse pcrit");
                        }
                        "vcrit" => {
                            self.vcrit = parsed.1.parse::<f64>()
                                .expect("cannot parse vcrit");
                            self.dcrit = 1.0 / self.vcrit;
                        }
                        "dcrit" => {
                            self.dcrit = parsed.1.parse::<f64>()
                                .expect("cannot parse dcrit");
                            self.vcrit = 1.0 / self.dcrit;
                        }
                        "molarvolume" => {
                            self.molarvolume = parsed.1.parse::<f64>()
                                .expect("cannot parse molarvolume");
                            self.molardensity = 1.0 / self.molarvolume;
                        }
                        "molardensity" => {
                            self.molardensity = parsed.1.parse::<f64>()
                                .expect("cannot parse molardensity");
                            self.molarvolume = 1.0 / self.molardensity;
                        }
                        "closescale" => {
                            self.closescale = parsed.1.parse::<f64>()
                                .expect("cannot parse closescale");
                        }
                        "velscale" => {
                            self.velscale = parsed.1.parse::<f64>()
                                .expect("cannot parse velscale");
                        }
                        "molfact" => {
                            self.molfact = parsed.1.parse::<usize>()
                                .expect("cannot parse molfact");
                        }
                        "itermax" => {
                            self.itermax = parsed.1.parse::<usize>()
                                .expect("cannot parse itermax");
                        }
                        "anneal" => {
                            self.anneal = parsed.1.parse::<isize>()
                                .expect("cannot parse anneal");
                        }
                        "sample" => {
                            self.sample = parsed.1.parse::<usize>()
                                .expect("cannot parse sample");
                        }
                        "ntrack" => {
                            self.ntrack = parsed.1.parse::<usize>()
                                .expect("cannot parse ntrack");
                        }
                        "ndist" => {
                            self.ndist = parsed.1.parse::<usize>()
                                .expect("cannot parse ndist");
                        }
                        "dtsteps" => {
                            self.dtsteps = parsed.1.parse::<usize>()
                                .expect("cannot parse dtsteps");
                        }
                        "evalmethod" => {
                            self.evalmethod = parsed.1.parse::<usize>()
                                .expect("cannot parse evalmethod");
                        }
                        "accmethod" => {
                            self.accmethod = parsed.1.parse::<usize>()
                                .expect("cannot parse accmethod");
                        }
                        "voltherm" => {
                            self.voltherm = parsed.1.parse::<usize>()
                                .expect("cannot parse voltherm");
                        }
                        "logaves" => {
                            self.logaves = parsed.1.parse::<usize>()
                                .expect("cannot parse logaves");
                        }
                        "periodic" => {
                            self.periodic = parsed.1.parse::<usize>()
                                .expect("cannot parse periodic");
                        }
                        "nworker" => {
                            self.nworker = parsed.1.parse::<usize>()
                                .expect("cannot parse nworker");
                        }
                        _ => {
                            if !parsable.contains(";") && !parsable.is_empty() {
                                println!("bad input on line {linecount}");
                                process::exit(1);
                            }
                        }
                    }
                }
                linecount += 1;
            }
        }
    }

    fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where
        P: AsRef<Path>,
    {
        let file = File::open(filename)?;
        Ok(io::BufReader::new(file).lines())
    }

    fn remove_whitespace_and_comments(
        s: &str,
    ) -> String {
        let tmp1: String = s.chars().filter(|c| !c.is_whitespace()).collect();
        let tmp2: Vec<&str> = tmp1.split(';').collect();
        if tmp1.contains(";") && tmp2.len() == 1 {
            "".to_string()
        } else {
            tmp2[0].to_string()
        }
    }

    fn trim_line(
        line: String,
    ) -> (String, String) {
        let tmp: Vec<&str> = line.split('=').collect();
        if tmp.len() > 1 {
            let tmptail: Vec<&str> = tmp[1].split_whitespace().collect();
            (tmp[0].to_string(), tmptail[0].to_string())
        } else {
            (tmp[0].to_string(), "".to_string())
        }
    }

    pub fn dump_state(
        &self,
        filename: String,
        coords: Vec<Simd<f64, LANES>>,
        vels: Vec<Simd<f64, LANES>>,
    ) {
        let mut file = OpenOptions::new()
            .append(true)
            .create(true)
            .open(filename)
            .expect("trouble openning the state file.");

        writeln!(file, "basename = {}", self.basename)
            .expect("having trouble writing to state file");
        writeln!(file, "element = {}", self.element)
            .expect("having trouble writing to state file");
        writeln!(file, "sigma = {}", self.sigma)
            .expect("having trouble writing to state file");
        writeln!(file, "epsilon = {}", self.epsilon)
            .expect("having trouble writing to state file");
        writeln!(file, "mass = {}", self.mass)
            .expect("having trouble writing to state file");
        writeln!(file, "tcrit = {}", self.tcrit)
            .expect("having trouble writing to state file");
        writeln!(file, "temperature = {}", self.temperature)
            .expect("having trouble writing to state file");
        writeln!(file, "pcrit = {}", self.pcrit)
            .expect("having trouble writing to state file");
        writeln!(file, "vcrit = {}", self.vcrit)
            .expect("having trouble writing to state file");
        writeln!(file, "dcrit = {}", self.dcrit)
            .expect("having trouble writing to state file");
        writeln!(file, "molarvolume = {}", self.molarvolume)
            .expect("having trouble writing to state file");
        writeln!(file, "molardensity = {}", self.molardensity)
            .expect("having trouble writing to state file");
        writeln!(file, "closescale = {}", self.closescale)
            .expect("having trouble writing to state file");
        writeln!(file, "velscale = {}", self.velscale)
            .expect("having trouble writing to state file");
        writeln!(file, "molfact = {}", self.molfact)
            .expect("having trouble writing to state file");
        writeln!(file, "itermax = {}", self.itermax)
            .expect("having trouble writing to state file");
        writeln!(file, "anneal = {}", self.anneal)
            .expect("having trouble writing to state file");
        writeln!(file, "sample = {}", self.sample)
            .expect("having trouble writing to state file");
        writeln!(file, "ntrack = {}", self.ntrack)
            .expect("having trouble writing to state file");
        writeln!(file, "ndist = {}", self.ndist)
            .expect("having trouble writing to state file");
        writeln!(file, "dtsteps = {}", self.dtsteps)
            .expect("having trouble writing to state file");
        writeln!(file, "evalmethod = {}", self.evalmethod)
            .expect("having trouble writing to state file");
        writeln!(file, "accmethod = {}", self.accmethod)
            .expect("having trouble writing to state file");
        writeln!(file, "voltherm = {}", self.voltherm)
            .expect("having trouble writing to state file");
        writeln!(file, "logaves = {}", self.logaves)
            .expect("having trouble writing to state file");
        writeln!(file, "periodic = {}", self.periodic)
            .expect("having trouble writing to state file");
        writeln!(file, "nworker = {}", self.nworker)
            .expect("having trouble writing to state file");
        for i in (0..).take(1<<(3 * self.molfact)) {
            writeln!(file, "{}, {}, {}",
                     coords[i][0],
                     coords[i][1],
                     coords[i][2],
            ).expect("trouble writing coordinates");
        }
        for i in (0..).take(1<<(3 * self.molfact)) {
            writeln!(file, "{}, {}, {}",
                     vels[i][0],
                     vels[i][1],
                     vels[i][2],
            ).expect("trouble writing velocities");
        }
    }

    pub fn restore_state(
        &mut self,
        infilename: String,
    ) {
        let mut linecount: usize = 0;
        let nparms: usize = 27 ;
        let mut nmol: usize = 0 ;

        println!("restoring state from {infilename}");
        
        if let Ok(lines) = SimConfig::read_lines(infilename) {
            for line in lines {
                if let Ok(parsable) = line {
                    //since this is restarting from a state dump,
                    // there will be no variables.
                    if linecount < nparms {
                        let cleaned = SimConfig::remove_whitespace_and_comments(&parsable);
                        let parsed = SimConfig::trim_line(cleaned);
                        if parsed.1.is_empty() && !parsed.0.is_empty() {
                            println!("corrupted state dump. problem on line {linecount}");
                            process::exit(1);
                        }
                        match &parsed.0[..] {
                            "basename" => { // 00
                                self.basename = parsed.1.parse::<String>()
                                    .expect("cannot parse basename in restart");
                            }
                            "element" => { // 01
                                self.element = parsed.1.parse::<String>()
                                    .expect("cannot parse basename in restart");
                            }
                            "sigma" => { // 02
                                self.sigma = parsed.1.parse::<f64>()
                                    .expect("cannot parse sigma in restart");
                            }
                            "epsilon" => { // 03
                                self.epsilon = parsed.1.parse::<f64>()
                                    .expect("cannot parse epsilon in restart");
                            }
                            "mass" => { // 04
                                self.mass = parsed.1.parse::<f64>()
                                    .expect("cannot parse mass in restart");
                            }
                            "tcrit" => { // 05
                                self.tcrit = parsed.1.parse::<f64>()
                                    .expect("cannot parse tcrit in restart");
                            }
                            "temperature" => { // 06
                                self.temperature = parsed.1.parse::<f64>()
                                    .expect("cannot parse temperature in restart");
                            }
                            "pcrit" => { // 07
                                self.pcrit = parsed.1.parse::<f64>()
                                    .expect("cannot parse pcrit in restart");
                            }
                            "vcrit" => { // 08
                                self.vcrit = parsed.1.parse::<f64>()
                                    .expect("cannot parse vcrit in restart");
                            }
                            "dcrit" => { // 09
                                self.dcrit = parsed.1.parse::<f64>()
                                    .expect("cannot parse dcrit in restart");
                            }
                            "molarvolume" => { // 10
                                self.molarvolume = parsed.1.parse::<f64>()
                                    .expect("cannot parse molarvolume in restart");
                            }
                            "molardensity" => { // 11
                                self.molardensity = parsed.1.parse::<f64>()
                                    .expect("cannot parse molardensity in restart");
                            }
                            "closescale" => { // 12
                                self.closescale = parsed.1.parse::<f64>()
                                    .expect("cannot parse closescale in restart");
                            }
                            "velscale" => { // 13
                                self.velscale = parsed.1.parse::<f64>()
                                    .expect("cannot parse velscale in restart");
                            }
                            "molfact" => { // 14
                                self.molfact = parsed.1.parse::<usize>()
                                    .expect("cannot parse molfact in restart");
                                nmol = 1 << (3 * self.molfact);
                            }
                            "itermax" => { // 15
                                self.itermax = parsed.1.parse::<usize>()
                                    .expect("cannot parse itermax in restart");
                            }
                            "anneal" => { // 16
                                self.anneal = parsed.1.parse::<isize>()
                                    .expect("cannot parse anneal in restart");
                            }
                            "sample" => { // 17
                                self.sample = parsed.1.parse::<usize>()
                                    .expect("cannot parse sample in restart");
                            }
                            "ntrack" => { // 18
                                self.ntrack = parsed.1.parse::<usize>()
                                    .expect("cannot parse ntrack in restart");
                            }
                            "ndist" => { // 19
                                self.ndist = parsed.1.parse::<usize>()
                                    .expect("cannot parse ndist in restart");
                            }
                            "dtsteps" => { // 20
                                self.dtsteps = parsed.1.parse::<usize>()
                                    .expect("cannot parse dtsteps in restart");
                            }
                            "evalmethod" => { // 21
                                self.evalmethod = parsed.1.parse::<usize>()
                                    .expect("cannot parse evalmethod in restart");
                            }
                            "accmethod" => { // 22
                                self.accmethod = parsed.1.parse::<usize>()
                                    .expect("cannot parse accmethod in restart");
                            }
                            "voltherm" => { // 23
                                self.voltherm = parsed.1.parse::<usize>()
                                    .expect("cannot parse voltherm in restart");
                            }
                            "logaves" => { // 24
                                self.logaves = parsed.1.parse::<usize>()
                                    .expect("cannot parse logaves in restart");
                            }
                            "periodic" => { // 25
                                self.periodic = parsed.1.parse::<usize>()
                                    .expect("cannot parse periodic in restart");
                            }
                            "nworker" => { // 26
                                self.nworker = parsed.1.parse::<usize>()
                                    .expect("cannot parse nworker in restart");
                            }
                            _ => {
                            }
                        }
                    } else if linecount < (nparms + nmol) {
                        let cleaned = SimConfig::remove_whitespace_and_comments(&parsable);
                        let pruned: Vec<&str> = cleaned.split(',').collect();
                        let mut parsed: Vec<f64> = vec![0.0f64; 3];
                        for i in (0..).take(3) {
                            parsed[i] = pruned[i].parse::<f64>()
                                .expect("problem reading coord[{i}] on line {linecount}"); 
                        }
                        self.coords.push(Simd::<f64,LANES>::from_array(
                            [parsed[0], parsed[1], parsed[2], 0.0]));
                    } else {
                        let cleaned = SimConfig::remove_whitespace_and_comments(&parsable);
                        let pruned: Vec<&str> = cleaned.split(',').collect();
                        let mut parsed: Vec<f64> = vec![0.0f64; 3];
                        for i in (0..).take(3) {
                            parsed[i] = pruned[i].parse::<f64>()
                                .expect("problem reading vel[{i}] on line {linecount}");
                        }
                        self.vels.push(Simd::<f64,LANES>::from_array(
                            [parsed[0], parsed[1], parsed[2], 0.0]));
                    }
                }
                linecount += 1;
                self.readin = true;
            }
        }
    }

    pub fn restore_phase (
        &mut self,
        infilename: String,
    ) {
        let mut linecount: usize = 0;
        let nparms: usize = 27 ;
        let mut nmol: usize = 1<<(3 * self.molfact) ;

        println!("restoring phase from {infilename}");
        
        if let Ok(lines) = SimConfig::read_lines(infilename) {
            for line in lines {
                if let Ok(parsable) = line {
                    if linecount < nparms {
                        let cleaned = SimConfig::remove_whitespace_and_comments(&parsable);
                        let parsed = SimConfig::trim_line(cleaned);
                        if parsed.1.is_empty() && !parsed.0.is_empty() {
                            println!("corrupted state dump. problem on line {linecount}");
                            process::exit(1);
                        }
                        if &parsed.0[..] == "molfact" { //14
                            self.molfact = parsed.1.parse::<usize>()
                                .expect("cannot parse molfact in restart");
                            nmol = 1 << (3 * self.molfact);
                        }
                    } else if linecount < (nparms + nmol) {
                        let cleaned = SimConfig::remove_whitespace_and_comments(&parsable);
                        let pruned: Vec<&str> = cleaned.split(',').collect();
                        let mut parsed: Vec<f64> = vec![0.0f64; 3];
                        for i in (0..).take(3) {
                            parsed[i] = pruned[i].parse::<f64>()
                                .expect("problem reading coord[{i}] on line {linecount}");
                        }
                        self.coords.push(Simd::<f64,LANES>::from_array(
                            [parsed[0], parsed[1], parsed[2], 0.0]));
                    } else {
                        let cleaned = SimConfig::remove_whitespace_and_comments(&parsable);
                        let pruned: Vec<&str> = cleaned.split(',').collect();
                        let mut parsed: Vec<f64> = vec![0.0f64; 3];
                        for i in (0..).take(3) {
                            parsed[i] = pruned[i].parse::<f64>()
                                .expect("problem reading vel[{i}] on line {linecount}");
                        }
                        self.vels.push(Simd::<f64,LANES>::from_array(
                            [parsed[0], parsed[1], parsed[2], 0.0]));
                    }
                }
                linecount += 1;
                self.readin = true;
            }
        }
    }
}

#[derive(Clone)]
pub struct VarConfig {
    pub name: String,
    min: f64,
    max: f64,
    steps: usize,
    linear: bool,
    oscillating: bool,
    increasing: bool,
    pub curval: f64,
}

impl Default for VarConfig {
    fn default() -> Self
    {
        VarConfig {
            name: "temperature".to_string(),
            min: 5.0,
            max: 50.0,
            steps: 45,
            linear: true,
            oscillating: true,
            increasing: true,
            curval: 5.0,
        }
    }
}

//should re-impliment VarConfig with generics, but for now, do the
// conversion in the update loop.

impl VarConfig {
    fn parse_input_line(
        &mut self,
        input_string: &String,
    ) {
        let tmp: Vec<&str> = input_string
            .split_whitespace()
            .collect();
        if tmp.len() != 8 {
            println!("bad variable specification:\n{input_string}");
            process::exit(1);
        } else {
            self.name = tmp[1].to_string();
            println!("parsing variable conditions for {}:",
                     self.name,
            );
            self.min = tmp[2].parse::<f64>()
                .expect("trouble parsing min for variable.");
            self.max = tmp[3].parse::<f64>()
                .expect("trouble parsing max for variable.");
            self.steps = tmp[4].parse::<usize>()
                .expect("trouble parsing steps for variable.");
            self.linear = tmp[5].contains("lin") || tmp[5].contains("true");
            self.oscillating = tmp[6].contains("osc") || tmp[6].contains("true");
            if self.min > self.max {
                (self.min, self.max) = (self.max, self.min);
            }
            if (self.max - self.min).abs() < f64::EPSILON {
                println!("min and max are too similar");
                process::exit(1);
            }
            if tmp[7].contains("low") {
                self.curval = self.min;
                self.increasing = true;
            } else {
                self.curval = self.max;
                self.increasing = false;
            }
        }
    }

    pub fn get_current_value(
        &self,
    ) -> (String, f64) {
        (self.name.clone(), self.curval)
    }
    
    pub fn get_next_value(
        &mut self,
    ) -> (String, f64, bool) {
        let mut flipped: bool = false;
        let oldval: f64 = self.curval ;
        let dir: f64 = if self.increasing { 1.0 } else { -1.0 };
        let fact: f64;
        if self.linear {
            fact = dir * (self.max - self.min)/(self.steps as f64);
            self.curval += fact;
        } else {
            fact = (dir * ((self.max/self.min).ln())/self.steps as f64).exp();
            self.curval *= fact;
        }
        if self.increasing && self.curval > self.max {
            if self.oscillating {
                self.curval = oldval;
                self.increasing = !self.increasing;
            } else {
                self.curval = self.min;
            }
            flipped = true;
        }
        if !self.increasing && self.curval < self.min {
            if self.oscillating {
                self.curval = oldval;
                self.increasing = !self.increasing;
            } else {
                self.curval = self.max;
            }
            flipped = true;
        }

        (self.name.clone(), self.curval, flipped)
    }
}
