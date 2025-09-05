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

mod graphparameters;
mod physicalparameters;
mod textparameters;

use crate::eyecandy::graphparameters::*;
use crate::eyecandy::physicalparameters::*;
use crate::eyecandy::textparameters::*;
use crate::fileio::*;

use bevy::{
    prelude::*,
    color::palettes::basic::{SILVER, WHITE},
    pbr::NotShadowReceiver,
};

use crossbeam_channel::{unbounded, Receiver, Sender};

use std::f32::consts::PI;
use std::simd::Simd;

use crate::channelstructs::*;
use crate::constants::*;

#[derive(Component)]
pub struct InfoText;

#[derive(Component)]
pub struct CommText;

#[derive(Component)]
pub struct GrText;

#[derive(Component)]
pub struct PvText;

#[derive(Component, Default)]
pub struct MyIndex(usize);

#[derive(Bundle, Default)]
struct BodyBundle {
    pbr: (Mesh3d, MeshMaterial3d<StandardMaterial>, Transform, NotShadowReceiver),
    myindex: MyIndex,
}

#[derive(Component, Default)]
pub struct BoxWall;

#[derive(Bundle, Default)]
struct PlaneBundle {
    pbr: (Mesh3d, MeshMaterial3d<StandardMaterial>, Transform),
    boxwall: BoxWall,
}

#[derive(Bundle, Default)]
struct CameraBundle {
    camera: (Camera3d, Transform)
}

#[derive(Bundle, Default)]
struct LightBundle {
    light: (PointLight, Transform),
}

#[derive(Resource)]
pub struct EyeCandy {
    physparm: PhysicalParameters,
    textparm: TextParameters,
    grparm: GraphParameters,
    pvparm: GraphParameters,
    coords_requ_channel: Sender<usize>,
    coords_send_channel: Receiver<TransferCoords>,
    text_requ_channel: Sender<usize>,
    text_send_channel: Receiver<MacroStateData>,
    dist_requ_channel: Sender<usize>,
    dist_send_channel: Receiver<TransferDist>,
}

impl Default for EyeCandy {
    fn default() -> Self {
        EyeCandy {
            physparm: PhysicalParameters::default(),
            textparm: TextParameters::default(),
            grparm: GraphParameters::default(),
            pvparm: GraphParameters::default(),
            coords_requ_channel: unbounded().0.clone(),
            coords_send_channel: unbounded().1.clone(),
            text_requ_channel: unbounded().0.clone(),
            text_send_channel: unbounded().1.clone(),
            dist_requ_channel: unbounded().0.clone(),
            dist_send_channel: unbounded().1.clone(),
        }
    }
}

pub struct EyeCandyFinish{
    pub coords_requ_channel: Sender<usize>,
    pub coords_send_channel: Receiver<TransferCoords>,
    pub text_requ_channel: Sender<usize>,
    pub text_send_channel: Receiver<MacroStateData>,
    pub dist_requ_channel: Sender<usize>,
    pub dist_send_channel: Receiver<TransferDist>,
    pub input_config: SimConfig,
}

impl EyeCandy {
    pub fn finish_setup(
        &mut self,
        finish_input: EyeCandyFinish,
    ) {
        self.coords_requ_channel = finish_input.coords_requ_channel.clone();
        self.coords_send_channel = finish_input.coords_send_channel.clone();
        self.text_requ_channel = finish_input.text_requ_channel.clone();
        self.text_send_channel = finish_input.text_send_channel.clone();
        self.dist_requ_channel = finish_input.dist_requ_channel.clone();
        self.dist_send_channel = finish_input.dist_send_channel.clone();

        let nmol: usize = 1 << (3 * finish_input.input_config.molfact);
        
        self.physparm.nmol = nmol;
        self.physparm.ballcoords = vec![[0.0f32;3];nmol];

        self.update_coordinates();
        self.update_textinfo();
        self.textparm.displayactive = finish_input.input_config.displayactive;
        println!("displayactive set to {}", self.textparm.displayactive);
    }

    pub fn update_coordinates(
        &mut self,
    ) {
        self.coords_requ_channel.send(1)
            .expect("trouble sending coords request from eyecandy");
        let coords_record: TransferCoords =
            self.coords_send_channel.recv().expect(
                "Problem when recieving coordinate data from channel"
            );

        let localcoords: *const Simd<f64, LANES> =
            coords_record.coords as *const Simd<f64, LANES> ;

        (0..).take(coords_record.nmol).for_each(|i| {
            (0..).take(3usize).for_each(|k| {
                unsafe {
                    self.physparm.ballcoords[i][k] =
                        ((&(*localcoords.add(i)))[k] * self.physparm.bevyscale) as f32;
                }
            });
        });
    }

    pub fn update_textinfo(
        &mut self,
    ) {
        self.text_requ_channel.send(1)
            .expect("trouble sending text request from eyecandy");
        let data: MacroStateData = self.text_send_channel.recv().expect(
            "Problem when recieving data from the text send channel"
        );

        self.textparm.itcount = data.itcount;
        self.textparm.ittotal = data.ittotal;
        self.textparm.elapsedtime = data.elapsedtime;
        self.textparm.timetotal = data.timetotal;
        self.textparm.itmax = data.itmax;
        self.textparm.nmol = data.nmol;
        self.textparm.samples = data.samples;
        self.textparm.element = data.element;
        self.textparm.temperature = data.temperature;
        self.textparm.volume = data.volume * 1.0e27;
        self.physparm.bound = (data.bound * BEVYSCALE) as f32;
        self.textparm.dt = data.dt * 1.0e9;
        self.textparm.kinetic = data.kinetic;
        self.textparm.kinetic1mean = data.kinetic1mean;
        self.textparm.kinetic2mean = data.kinetic2mean;
        self.textparm.potential = data.potential;
        self.textparm.potential1mean = data.potential1mean;
        self.textparm.potential2mean = data.potential2mean;
        self.textparm.pressure = data.pressure;
        self.textparm.pressure1mean = data.pressure1mean;
        self.textparm.pressure2mean = data.pressure2mean;
        self.textparm.deltaqdeltat = data.deltaqdeltat;
        self.textparm.deltaqdeltat1mean = data.deltaqdeltat1mean;
        self.textparm.deltaqdeltat2mean = data.deltaqdeltat2mean;
        self.textparm.itspersec = data.itspersec;
        self.textparm.itspersec1mean = data.itspersec1mean;
        self.textparm.itspersec2mean = data.itspersec2mean;
        self.textparm.ntrackinv = data.ntrackinv;
        if self.textparm.displayactive != data.displayactive {
            self.textparm.displayactive = data.displayactive;
        }
    }

    pub fn update_box(
        eye: Res<EyeCandy>,
        mut query: Query<&mut Transform, With<BoxWall>>,
    ) {
        query.iter_mut().for_each(|mut item| {
            if item.translation.x.abs() > f32::EPSILON {
                item.translation.x = item.translation.x.signum() * eye.physparm.bound;
                item.scale = Vec3::splat(eye.physparm.bound);
            }
            if item.translation.y.abs() > f32::EPSILON {
                item.translation.y = item.translation.y.signum() * eye.physparm.bound;
                item.scale = Vec3::splat(eye.physparm.bound);
            }
            if item.translation.z.abs() > f32::EPSILON {
                item.translation.z = item.translation.z.signum() * eye.physparm.bound;
                item.scale = Vec3::splat(eye.physparm.bound);
            }
        });
    }

    pub fn update_camera(
        eye: Res<EyeCandy>,
        mut query: Query<&mut Transform, With<Camera>>,
    ) {
        query.iter_mut().for_each(|mut item| {
            item.translation.x = 3.5 * eye.physparm.bound;
        });
    }

    pub fn update_light(
        eye: Res<EyeCandy>,
        mut query: Query<&mut Transform, With<PointLight>>,
    ) {
        query.iter_mut().for_each(|mut item| {
            item.translation.x = 2.5 * eye.physparm.bound;
        });
    }

    pub fn generate_bodies(
        mut commands: Commands,
        mut meshes: ResMut<Assets<Mesh>>,
        mut materials: ResMut<Assets<StandardMaterial>>,
        eye: Res<EyeCandy>,
        _asset_server: Res<AssetServer>,
    ) {
        let mesh = meshes.add(
            Sphere::default()
                .mesh()
                .ico(2).expect(
                    "Cannot create the sphere mesh"
                )
        );

        for i in (0..).take(eye.physparm.nmol) {
            let p32: [f32; 3] = [
                eye.physparm.ballcoords[i][0],
                eye.physparm.ballcoords[i][1],
                eye.physparm.ballcoords[i][2],
            ];

            let base_color: Color = if i == 0 {
                Color::srgba(0.2, 0.5, 0.8, 1.0)
            } else {
                Color::srgba(0.8, 0.5, 0.2, 1.0)
            };

            commands.spawn(BodyBundle {
                pbr: (Mesh3d(mesh.clone()),
                      MeshMaterial3d(materials.add(StandardMaterial {
                          base_color,
//                          alpha_mode: AlphaMode::Blend,
                          ..default()
                      })),
                      Transform {
                          translation: Vec3::from_array(p32),
                          scale: Vec3::splat(1.0 * eye.physparm.ballradius),
                          ..default()
                      },
                      NotShadowReceiver,
                ),
                myindex: MyIndex(i),
            }
            );
        }

        let boxlen: f32 = eye.physparm.bound;
        let wallscale: f32 = 2.0;

        commands.spawn(PlaneBundle {
            pbr: (Mesh3d(meshes.add(Plane3d::default().mesh().size(wallscale,wallscale))),
                  MeshMaterial3d(materials.add(Color::from(SILVER))),
                  Transform{
                      translation: Vec3::from_array([0., -boxlen, 0.]),
                      rotation: Quat::from_rotation_z(0.),
                      ..default()
                  },
            ),
            boxwall: BoxWall,
        });
        commands.spawn(PlaneBundle {
            pbr: (Mesh3d(meshes.add(Plane3d::default().mesh().size(wallscale,wallscale))),
                  MeshMaterial3d(materials.add(Color::from(SILVER))),
                  Transform{
                      translation: Vec3::from_array([0., boxlen, 0.]),
                      rotation: Quat::from_rotation_z(PI),
                      ..default()
                  },
            ),
            boxwall: BoxWall,
        });
        commands.spawn(PlaneBundle {
            pbr: (Mesh3d(meshes.add(Plane3d::default().mesh().size(wallscale,wallscale))),
                  MeshMaterial3d(materials.add(Color::from(SILVER))),
                  Transform{
                      translation: Vec3::from_array([0., 0., boxlen]),
                      rotation: Quat::from_rotation_x(-0.5*PI),
                      ..default()
                  },
            ),
            boxwall: BoxWall,
        });
        commands.spawn(PlaneBundle {
            pbr: (Mesh3d(meshes.add(Plane3d::default().mesh().size(wallscale,wallscale))),
                  MeshMaterial3d(materials.add(Color::from(SILVER))),
                  Transform{
                      translation: Vec3::from_array([0., 0., -boxlen]),
                      rotation: Quat::from_rotation_x(0.5*PI),
                      ..default()
                  },
            ),
            boxwall: BoxWall,
        });
        commands.spawn(PlaneBundle {
            pbr: (Mesh3d(meshes.add(Plane3d::default().mesh().size(wallscale,wallscale))),
                  MeshMaterial3d(materials.add(Color::from(SILVER))),
                  Transform{
                      translation: Vec3::from_array([-boxlen, 0., 0.]),
                      rotation: Quat::from_rotation_z(-0.5*PI),
                      ..default()
                  },
            ),
            boxwall: BoxWall,
        });

        commands.spawn(CameraBundle {
            camera: (Camera3d::default(),
                     Transform::from_xyz(3.5 * boxlen, 0., 0.)
                     .looking_at(Vec3::ZERO, Vec3::Z),
            ),
        });

        commands.spawn(LightBundle {
            light: (PointLight {
                color: Color::WHITE,
                intensity: 800000. * (MOLFACT as f32),
                range: 6.0 * boxlen,
                radius: 2.0 * eye.physparm.ballradius,
                shadows_enabled: true,
                ..default()
            },
                    Transform::from_xyz(2.5 * boxlen, 0., 0.),
            ),
        });
    }

    pub fn generate_text(
        mut commands: Commands,
        _asset_server: Res<AssetServer>,
    ) {
        let root_uinode = commands
            .spawn(
                Node {
                    width: Val::Percent(100.),
                    height: Val::Percent(100.),
                    justify_content: JustifyContent::SpaceBetween,
                    ..default()
                },
            ).id();

        let left_column = commands.spawn((
            Text::default(),
            TextLayout::new_with_justify(JustifyText::Left),
            TextColor::from(Color::WHITE),
            Node {
                align_self: AlignSelf::FlexStart,
                position_type: PositionType::Absolute,
                top: Val::Px(5.0),
                left: Val::Px(5.0),
                ..default()
            },
            InfoText,
        )).with_children(|parent| {
            parent.spawn(
                TextSpan::new("Max iterations: 0"));
            parent.spawn(
                TextSpan::new("\nCurrent Iteration and time:  \n0\n0.0"));
            parent.spawn(
                TextSpan::new("\nSampling Rate: 0"));
            parent.spawn(
                TextSpan::new("\nElement:  \nname\nTc\nPc\nVc"));
            parent.spawn(
                TextSpan::new("\nNumber of Molecules: 0"));
            parent.spawn(
                TextSpan::new("\nTemperature: 0.0"));
            parent.spawn(
                TextSpan::new("\nBox Volume (nm^3): 0.0"));
            parent.spawn(
                TextSpan::new("\nMolar Volume (mL/mol): 0.0"));
            parent.spawn(
                TextSpan::new("\nTime Increment (ns): 0.0"));
            parent.spawn(
                TextSpan::new("\nKinetic Energy (kT/particle):  \n0.0\n0.0 +/- 0.0"));
            parent.spawn(
                TextSpan::new("\nPotential Energy (kt/particle):  \n0.0\n0.0 +/- 0.0"));
            parent.spawn(
                TextSpan::new("\nPressure (Pa):  \n0.0\n0.0 +/- 0.0"));
            parent.spawn(
                TextSpan::new("\nHeat Flow Into Box (J/s):  \n0.0\n0.0 +/- 0.0"));
            parent.spawn(
                TextSpan::new("\nIterations Per Second:  \n0.0\n0.0 +/- 0.0"));
        }).id();

        let right_column = commands.spawn((
            Text::default(),
            TextLayout::new_with_justify(JustifyText::Left),
            TextColor::from(Color::WHITE),
            Node {
                align_self: AlignSelf::FlexStart,
                position_type: PositionType::Absolute,
                top: Val::Px(5.0),
                right: Val::Px(5.0),
                ..default()
            },
            CommText,
        )).with_children(|parent| {
            parent.spawn(TextSpan::new("Keys for Actions:\n  Q: Quit"));
            parent.spawn(TextSpan::new("\n\n  E: Scale velocities by 1.022"));
            parent.spawn(TextSpan::new("\n  D: Scale velocities by 1.00/1.022"));
            parent.spawn(TextSpan::new("\n  T: Scale temperature by 1.022"));
            parent.spawn(TextSpan::new("\n  G: Scale temperature 1.00/1.022"));
            parent.spawn(TextSpan::new("\n  U: Scale bound by 1.022"));
            parent.spawn(TextSpan::new("\n  J: Scale bound by 1.00/1.022"));
            parent.spawn(TextSpan::new("\n  R: Scale the step size by 1/2"));
            parent.spawn(TextSpan::new("\n  F: Scale the step size by 2"));
            parent.spawn(TextSpan::new("\n\n  I: Toggle potentials on/off"));
            parent.spawn(TextSpan::new("\n  O: Toggle info display on/off"));
            parent.spawn(TextSpan::new("\n  A: Thermalize"));
        }).id();

        let _grplotlabel = commands.spawn((
            Text::default(),
            TextLayout::new_with_justify(JustifyText::Center),
            TextColor::from(Color::WHITE),
            Node {
                align_self: AlignSelf::FlexStart,
                position_type: PositionType::Absolute,
                top: Val::Px(420.0),
                right: Val::Px(30.0),
                ..default()
            },
            GrText,
        )).with_children(|parent| {
            parent.spawn(TextSpan::new("g(r)"));
            parent.spawn(TextSpan::new("\n\nPeak: r=0.0 nm, g=0.0"));
        }).id();

        let _pvplotlabel = commands.spawn((
            Text::default(),
            TextLayout::new_with_justify(JustifyText::Center),
            TextColor::from(Color::WHITE),
            Node {
                align_self: AlignSelf::FlexStart,
                position_type: PositionType::Absolute,
                top: Val::Px(740.0),
                right: Val::Px(30.0),
                ..default()
            },
            PvText,
        )).with_children(|parent| {
            parent.spawn(TextSpan::new("P(v)"));
            parent.spawn(TextSpan::new("\n\nPeak: v=0.0 m/s"));
        }).id();
        
        commands
            .entity(root_uinode)
            .add_children(&[left_column, right_column]);
    }
    
    pub fn update_parms(
        mut query: Query<(&mut Transform, &MyIndex)>,
        mut eye: ResMut<EyeCandy>,
    ) {
        eye.update_coordinates();

        query.iter_mut().for_each(|(mut transform, index)| {
            transform.translation = Vec3::from(
                <[f32; 3]>::try_from(&(eye.physparm.ballcoords[index.0])[0..3]).expect(
                    "cannot update coordinate in update_parms"
                ),
            );
        });
    }

    //The following does what it is supposed to do, but since shadows
    // are not impacted by the alpha channel of the object casting the
    // shadow, it leads to wonky side effects that actually make
    // things look worse. So, until Bevy addresses this issue, I think
    // I'll leave it out.

    pub fn _update_colors(
        query: Query<(&MeshMaterial3d<StandardMaterial>, &MyIndex)>,
        eye: Res<EyeCandy>,
        mut materials: ResMut<Assets<StandardMaterial>>,
    ) {
        let bound: f32 = eye.physparm.bound;
        let radius: f32 = eye.physparm.ballradius;
        let radius_1: f32 = 1.0/radius;
        
        query.iter().for_each(|(handle, index)| {
            let mut alpha: f32 = 1.0;
            for i in (0..).take(3usize) {
                let sep: f32 = bound - eye.physparm.ballcoords[index.0][i].abs();
                if sep < radius {
                    alpha *= sep * radius_1;
                }
            }
            if (1.0-alpha).abs() > f32::EPSILON {
                let material = materials.get_mut(handle)
                    .expect("can't get handle to change color");
                material.base_color.set_alpha(alpha);
            }
        });
    }
    
    pub fn draw_graphs(
        mut eye: ResMut<EyeCandy>,
        mut gizmos: Gizmos,
    ) {
        if eye.textparm.displayactive {
            let boxlen: f32 = eye.physparm.bound;
            let width: f32 = 0.78;
            let height: f32 = 0.44;
            let offset: f32 = 0.02;
            let y0: f32 = 1.0+offset;
            let y1: f32 = 1.0+offset+width;
            let mut z0: f32 = -0.6*height;
            let mut z1: f32 = z0+height;
            
            eye.dist_requ_channel.send(1)
                .expect("trouble asking for speeds");
            let speeddata: TransferDist = eye.dist_send_channel.recv()
                .expect("trouble recieving speeds");
            eye.dist_requ_channel.send(2)
                .expect("trouble asking for radial distribution");
            let radialdata: TransferDist = eye.dist_send_channel.recv()
                .expect("trouble recieving radial dist");
            
            gizmos.linestrip([
                Vec3::from_array([boxlen, y0*boxlen, z0*boxlen]),
                Vec3::from_array([boxlen, y1*boxlen, z0*boxlen]),
                Vec3::from_array([boxlen, y1*boxlen, z1*boxlen]),
                Vec3::from_array([boxlen, y0*boxlen, z1*boxlen]),
                Vec3::from_array([boxlen, y0*boxlen, z0*boxlen]),
            ],
                             WHITE,
            );
            
            let mut dx: f32 = radialdata.vals[radialdata.nval] as f32;
            let mut xmax: f32 = radialdata.vals[radialdata.nval+1] as f32;
            let mut ymin: f32 = radialdata.vals[radialdata.nval+2] as f32;
            let mut ymax: f32 = radialdata.vals[radialdata.nval+3] as f32;
            let mut yrng: f32 = ymax - ymin;
            let mut peakx: f32 = radialdata.vals[radialdata.nval+4] as f32;
            
            eye.grparm.dx = dx;
            eye.grparm.xmax = xmax;
            eye.grparm.ymin = ymin;
            eye.grparm.ymax = ymax;
            eye.grparm.yrng = yrng;
            eye.grparm.peakx = peakx;
            
            let mut verts: Vec<Vec3> = vec![Vec3::from_array([0.0f32; 3]);radialdata.nval];
            
            let zlevel: f32 = z0 + height * (1.0 - ymin) / yrng;
            
            (0..).take(radialdata.nval).for_each(|i| {
                verts[i] = Vec3::from_array([
                    boxlen,
                    boxlen * (y0 + width * (dx * i as f32) / xmax),
                    boxlen * (z0 + height * (radialdata.vals[i] as f32 - ymin)/yrng),
                ])
            });
            
            gizmos.linestrip(
                verts,
                WHITE,
            );
            
            gizmos.line(
                Vec3::from_array([boxlen, y0*boxlen, zlevel*boxlen]),
                Vec3::from_array([boxlen, y1*boxlen, zlevel*boxlen]),
                WHITE,
            );
            
            dx = speeddata.vals[speeddata.nval] as f32;
            xmax = speeddata.vals[speeddata.nval+1] as f32;
            ymin = speeddata.vals[speeddata.nval+2] as f32;
            ymax = speeddata.vals[speeddata.nval+3] as f32;
            yrng = ymax - ymin;
            peakx = speeddata.vals[speeddata.nval+4] as f32;
            
            eye.pvparm.dx = dx;
            eye.pvparm.xmax = xmax;
            eye.pvparm.ymin = ymin;
            eye.pvparm.ymax = ymax;
            eye.pvparm.yrng = yrng;
            eye.pvparm.peakx = peakx;
            
            let mut verts: Vec<Vec3> = vec![Vec3::from_array([0.0f32; 3]);speeddata.nval];
            
            z0 = -2.0 * height;
            z1 = z0 + height;
            
            gizmos.linestrip([
                Vec3::from_array([boxlen, y0*boxlen, z0*boxlen]),
                Vec3::from_array([boxlen, y1*boxlen, z0*boxlen]),
                Vec3::from_array([boxlen, y1*boxlen, z1*boxlen]),
                Vec3::from_array([boxlen, y0*boxlen, z1*boxlen]),
                Vec3::from_array([boxlen, y0*boxlen, z0*boxlen]),
            ],
                             WHITE,
            );
            
            (0..).take(speeddata.nval).for_each(|i| {
                verts[i] = Vec3::from_array([
                    boxlen,
                    boxlen * (y0 + width * (dx * i as f32) / xmax),
                    boxlen * (z0 + height * (speeddata.vals[i] as f32 - ymin)/yrng),
                ])
            });
            gizmos.linestrip(
                verts,
                WHITE,
            );
        }
    }
    
    pub fn update_text_display(
        mut eye: ResMut<EyeCandy>,
        infotext: Query<Entity, With<InfoText>>,
        grtext: Query<Entity, With<GrText>>,
        pvtext: Query<Entity, With<PvText>>,
        commtext: Query<Entity, With<CommText>>,
        mut writer: TextUiWriter,
    ) {
        eye.update_textinfo();

        let vmolar: f64 = (eye.textparm.volume / eye.textparm.nmol as f64)
            * AVOGADRO
            * 1.0e-27
            * 1.0e6;
        
        if eye.textparm.displayactive {
            for entity in &grtext {
                writer.for_each_color(entity, |mut text_color| {
                    text_color.0 = Color::WHITE;
                });
            }
            for entity in &pvtext {
                writer.for_each_color(entity, |mut text_color| {
                    text_color.0 = Color::WHITE;
                });
            }
            for entity in &infotext {
                writer.for_each_color(entity, |mut text_color| {
                    text_color.0 = Color::WHITE;
                });
            }
            for entity in &commtext {
                writer.for_each_color(entity, |mut text_color| {
                    text_color.0 = Color::WHITE;
                });
            }
        } else {
            for entity in &grtext {
                writer.for_each_color(entity, |mut text_color| {
                    text_color.0 = Color::BLACK;
                });
            }
            for entity in &pvtext {
                writer.for_each_color(entity, |mut text_color| {
                    text_color.0 = Color::BLACK;
                });
            }
            for entity in &infotext {
                writer.for_each_color(entity, |mut text_color| {
                    text_color.0 = Color::BLACK;
                });
            }
            for entity in &commtext {
                writer.for_each_color(entity, |mut text_color| {
                    text_color.0 = Color::BLACK;
                });
            }
        }

        for entity in &grtext {
            *writer.text(entity, 2) = format!(
                "\n\nPeak: r={:.3e} nm, g={:.3e}",
                eye.grparm.peakx*1.0e9,
                eye.grparm.ymax,
            );
        }

        for entity in &pvtext {
            *writer.text(entity, 2) = format!(
                "\n\nPeak: v={:.3e} m/s",
                eye.pvparm.peakx,
            );
        }
        
        for entity in &infotext {
            *writer.text(entity, 1) = format!(
                "Max iterations: {:07}",
                eye.textparm.itmax,
            );
            *writer.text(entity, 2) = format!(
                "\n\nCurrent Iteration and time:\n{:07} ({:07})\n{:.05e} ({:.05e}) s",
                eye.textparm.itcount,
                eye.textparm.ittotal,
                eye.textparm.elapsedtime,
                eye.textparm.timetotal,
            );
            *writer.text(entity, 3) = format!(
                "\n\nSampling Rate: {}",
                eye.textparm.samples,
            );
            *writer.text(entity, 4) = format!(
                "\n\n{}",
                eye.textparm.element,
            );
            *writer.text(entity, 5) = format!(
                "\n\nNumber of Molecules: {}",
                eye.textparm.nmol,
            );
            *writer.text(entity, 6) = format!(
                "\n\nTemperature (K): {:.05e}",
                eye.textparm.temperature,
            );
            *writer.text(entity, 7) = format!(
                "\n\nBox Volume (nm^3): {:.05e}",
                eye.textparm.volume,
            );
            *writer.text(entity, 8) = format!(
                "\n\nMolar Volume (mL/mol): {vmolar:.05e}",
            );
            *writer.text(entity, 9) = format!(
                "\n\nTime Increment (ns): {:.05e}",
                eye.textparm.dt,
            );
            *writer.text(entity, 10) = format!(
                "\n\nKinetic Energy (kT/particle): \n{:.05e}\n{:.05e} +/- {:.05e}",
                eye.textparm.kinetic,
                eye.textparm.kinetic1mean,
                ((eye.textparm.kinetic2mean
                  - eye.textparm.kinetic1mean
                  * eye.textparm.kinetic1mean)
                 * eye.textparm.ntrackinv).sqrt(),
            );
            *writer.text(entity, 11) = format!(
                "\n\nPotential Energy (kT/particle): \n{:.05e}\n{:.05e} +/- {:.05e}",
                eye.textparm.potential,
                eye.textparm.potential1mean,
                ((eye.textparm.potential2mean
                  - eye.textparm.potential1mean
                  * eye.textparm.potential1mean)
                 * eye.textparm.ntrackinv).sqrt(),
            );
            *writer.text(entity, 12) = format!(
                "\n\nPressure (Pa): \n{:.05e}\n{:.05e} +/- {:.05e}",
                eye.textparm.pressure,
                eye.textparm.pressure1mean,
                ((eye.textparm.pressure2mean
                  - eye.textparm.pressure1mean
                  * eye.textparm.pressure1mean)
                 * eye.textparm.ntrackinv).sqrt(),
            );
            *writer.text(entity, 13) = format!(
                "\n\nHeat Flow in (J/s): \n{:.05e}\n{:.05e} +/- {:.05e}",
                eye.textparm.deltaqdeltat,
                eye.textparm.deltaqdeltat1mean,
                ((eye.textparm.deltaqdeltat2mean
                  - eye.textparm.deltaqdeltat1mean
                  * eye.textparm.deltaqdeltat1mean)
                 * eye.textparm.ntrackinv).sqrt(),
            );
            *writer.text(entity, 14) = format!(
                "\n\nIterations per second (s^-1): \n{:.05e}\n{:.05e} +/- {:.05e}\n\n\n",
                eye.textparm.itspersec,
                eye.textparm.itspersec1mean,
                ((eye.textparm.itspersec2mean
                  - eye.textparm.itspersec1mean
                  * eye.textparm.itspersec1mean)
                 * eye.textparm.ntrackinv).sqrt(),
            );
        }
    }
}

           
