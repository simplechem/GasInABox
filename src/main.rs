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
/*
#![cfg_attr(feature = "nightly", feature(portable_simd))]
 */
#![feature(portable_simd)]
//note, to switch to nightly builds, we need to run
//
//  rustup override set nightly
//
//to switch back (if simd goes into the main system) run
//
//  rustup override unset nightly
//
//alternatively you can use either of the following to compile with
// nightly or standard:
//
//  cargo +nightly build
//  cargo +standard build
//
//to cross compile for Windows. Note if it ever fails strangely, "rm
// -rf Cargo.lock target" to clean things out.
//
//  cargo build --target x86_64-pc-windows-gnu --release

mod averages;
mod channelstructs;
mod constants;
mod eyecandy;
mod fileio;
mod gasbox;
//mod fakesimd;

use bevy::{
    prelude::MonitorSelection::Primary,
    prelude::*,
    window::{PresentMode, WindowMode, WindowTheme},
    winit::WinitSettings,
};
use clap::Parser;

use crate::channelstructs::*;
use crate::eyecandy::EyeCandy;
use crate::fileio::SimConfig;
use crate::gasbox::GasBox;

use std::process;
use std::thread;

use crossbeam_channel::{Receiver, Sender, unbounded};

/*
Look at

https://github.com/bevyengine/bevy/discussions/1150

to see how to link channels between bevy and nonbevy stuff.

1) Create the channels.

2) Pass some of them into bevy as part of the set process of the bevy
   app.

3) Pass some of them into the EyeCandy implimentation as part of that
   setup.

 */

fn main() {
    let args = Args::parse();

    let mut input_config: SimConfig = SimConfig::default();
    if !args.infile.is_empty() {
        input_config.get_settings_from_file(args.infile.clone());
    }
    if !args.restart.is_empty() {
        input_config.restore_state(args.restart.clone());
    }
    if args.screensaver {
        input_config.displayactive = false;
        println!("screensaver mode set to false");
    } else {
        input_config.displayactive = true;
        println!("screensaver mode set to true");
    }
    
    let coords_requ_channel: (Sender<usize>, Receiver<usize>) = unbounded();
    let coords_send_channel: (Sender<TransferCoords>, Receiver<TransferCoords>) = unbounded();
    let text_requ_channel: (Sender<usize>, Receiver<usize>) = unbounded();
    let text_send_channel: (Sender<MacroStateData>, Receiver<MacroStateData>) = unbounded();
    let dist_requ_channel: (Sender<usize>, Receiver<usize>) = unbounded();
    let dist_send_channel: (Sender<TransferDist>, Receiver<TransferDist>) = unbounded();

    let mut gasbox: GasBox = GasBox::default();

    gasbox.finish_setup(
        gasbox::GasBoxFinish{
            coords_requ_channel: coords_requ_channel.1.clone(),
            coords_send_channel: coords_send_channel.0.clone(),
            text_requ_channel: text_requ_channel.1.clone(),
            text_send_channel: text_send_channel.0.clone(),
            dist_requ_channel: dist_requ_channel.1.clone(),
            dist_send_channel: dist_send_channel.0.clone(),
            config: input_config.clone(),
        }
    );

    let gasbox_handler = thread::spawn(move || {
        gasbox.run_iterations();
    });

    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "GasBox Display Window".into(),
                resolution: (1080., 1080.).into(),
                present_mode: PresentMode::AutoVsync,
                fit_canvas_to_parent: true,
                window_theme: Some(WindowTheme::Dark),
                mode: WindowMode::BorderlessFullscreen(Primary),
                ..default()
            }),
            ..default()
        }))
        .insert_resource(WinitSettings::game())
        .insert_resource(AmbientLight {
            brightness: 200.,
            ..default()
        })
        .insert_resource(ClearColor(Color::BLACK))
        .init_resource::<EyeCandy>()
        .insert_resource(ChannelCollection{
            coords_requ_channel: coords_requ_channel.0.clone(),
            coords_send_channel: coords_send_channel.1.clone(),
            text_requ_channel: text_requ_channel.0.clone(),
            text_send_channel: text_send_channel.1.clone(),
            dist_requ_channel: dist_requ_channel.0.clone(),
            dist_send_channel: dist_send_channel.1.clone(),
        })
        .insert_resource(SimulationConfig(input_config.clone()))
        .add_systems(Startup, connect_channels)
        .add_systems(
            Startup,
            (
                EyeCandy::generate_bodies.after(connect_channels),
                EyeCandy::generate_text.after(connect_channels),
            ),
        )
        .add_systems(
            Update,
            (
                EyeCandy::update_parms,
//                EyeCandy::update_colors.after(EyeCandy::update_parms),
                EyeCandy::update_text_display.after(EyeCandy::update_parms),
                EyeCandy::update_box.after(EyeCandy::update_parms),
                EyeCandy::update_camera.after(EyeCandy::update_parms),
                EyeCandy::update_light.after(EyeCandy::update_parms),
                EyeCandy::draw_graphs.after(EyeCandy::update_text_display),
                keyboard_input_system,
            ),
        )
        .run();

    gasbox_handler.join().expect("Gasbox can't leave!");
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    //name of an input file.
    #[arg(short, long, default_value="")]
    infile: String,

    //name of a file to restart from
    #[arg(short, long, default_value="")]
    restart: String,

    //screensaver mode
    #[arg(short, long, action)]
    screensaver: bool,
}

#[derive(Resource)]
pub struct ChannelCollection{
    coords_requ_channel: Sender<usize>,
    coords_send_channel: Receiver<TransferCoords>,
    text_requ_channel: Sender<usize>,
    text_send_channel: Receiver<MacroStateData>,
    dist_requ_channel: Sender<usize>,
    dist_send_channel: Receiver<TransferDist>,
}

#[derive(Resource)]
pub struct SimulationConfig(SimConfig);

pub fn connect_channels(
    mut eye: ResMut<EyeCandy>,
    channel_collection: Res<ChannelCollection>,
    config: Res<SimulationConfig>,
) {
    eye.finish_setup(
        eyecandy::EyeCandyFinish{
            coords_requ_channel: channel_collection.coords_requ_channel.clone(),
            coords_send_channel: channel_collection.coords_send_channel.clone(),
            text_requ_channel: channel_collection.text_requ_channel.clone(),
            text_send_channel: channel_collection.text_send_channel.clone(),
            dist_requ_channel: channel_collection.dist_requ_channel.clone(),
            dist_send_channel: channel_collection.dist_send_channel.clone(),
            input_config: config.0.clone(),
        }
    );
}

fn keyboard_input_system(
    keyboard_input: Res<ButtonInput<KeyCode>>,
    channel_collection: Res<ChannelCollection>,
) {
    //figure out hot to do this with match
    //energy up
    if keyboard_input.just_pressed(KeyCode::KeyD) {
        channel_collection.coords_requ_channel
            .send(2)
            .expect("main had trouble sending 2 down coords requ channel");
    }
    //energy down
    if keyboard_input.just_pressed(KeyCode::KeyE) {
        channel_collection.coords_requ_channel
            .send(3)
            .expect("main had trouble sending 3 down coords requ channel");
    }
    //temperature up
    if keyboard_input.just_pressed(KeyCode::KeyT) {
        channel_collection.coords_requ_channel
            .send(4)
            .expect("main had trouble sending 4 down coords requ channel");
    }
    //temperature down
    if keyboard_input.just_pressed(KeyCode::KeyG) {
        channel_collection.coords_requ_channel
            .send(5)
            .expect("main had trouble sending 5 down coords requ channel");
    }
    //dtsteps bigger
    if keyboard_input.just_pressed(KeyCode::KeyR) {
        channel_collection.coords_requ_channel
            .send(6)
            .expect("main had trouble sending 6 down coords requ channel");
    }
    //dtsteps smaller
    if keyboard_input.just_pressed(KeyCode::KeyF) {
        channel_collection.coords_requ_channel
            .send(7)
            .expect("main had trouble sending 7 down coords requ channel");
    }
    //box bigger
    if keyboard_input.just_pressed(KeyCode::KeyU) {
        channel_collection.coords_requ_channel
            .send(8)
            .expect("main had trouble sending 8 down coords requ channel");
    }
    //box smaller
    if keyboard_input.just_pressed(KeyCode::KeyJ) {
        channel_collection.coords_requ_channel
            .send(9)
            .expect("main had trouble sending 9 down coords requ channel");
    }
    //toggle potential
    if keyboard_input.just_pressed(KeyCode::KeyI) {
        channel_collection.coords_requ_channel
            .send(10)
            .expect("main had trouble sending 10 down coords requ channel");
    }
    //autotherm
    if keyboard_input.just_pressed(KeyCode::KeyA) {
        channel_collection.coords_requ_channel
            .send(11)
            .expect("main had trouble sending 11 down coords requ channel");
    }
    //toggle potential
    if keyboard_input.just_pressed(KeyCode::KeyO) {
        channel_collection.coords_requ_channel
            .send(12)
            .expect("main had trouble sending 12 down coords requ channel");
    }
    //quit
    if keyboard_input.just_pressed(KeyCode::KeyQ) {
        process::exit(0);
    }
}
