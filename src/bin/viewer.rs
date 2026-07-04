use eframe::egui;
use onedpipes::{OrganPipeConfig, OrganPipeRun, ScalarField, Snapshot};
use plotters::prelude::*;

const PLOT_W: u32 = 820;
const PLOT_H: u32 = 300;

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([1180.0, 900.0]),
        ..Default::default()
    };
    eframe::run_native(
        "onedpipes validation viewer",
        options,
        Box::new(|cc| Ok(Box::new(ViewerApp::new(cc)))),
    )
}

struct ViewerApp {
    run: OrganPipeRun,
    running: bool,
    steps_per_frame: usize,
    selected_field: ScalarField,
    config: OrganPipeConfig,
}

impl ViewerApp {
    fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        let config = OrganPipeConfig::default();
        Self {
            run: OrganPipeRun::new(config),
            running: false,
            steps_per_frame: 25,
            selected_field: ScalarField::Pressure,
            config,
        }
    }

    fn reset(&mut self) {
        self.run = OrganPipeRun::new(self.config);
    }

    fn advance(&mut self, steps: usize) {
        for _ in 0..steps {
            self.run.step();
        }
    }
}

impl eframe::App for ViewerApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        if self.running {
            self.advance(self.steps_per_frame);
            ctx.request_repaint();
        }

        egui::TopBottomPanel::top("top").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.heading("onedpipes validation viewer");
                ui.separator();
                if ui
                    .button(if self.running { "Pause" } else { "Run" })
                    .clicked()
                {
                    self.running = !self.running;
                }
                if ui.button("Step").clicked() {
                    self.advance(1);
                }
                if ui.button("Reset").clicked() {
                    self.reset();
                }
            });
        });

        egui::SidePanel::left("controls")
            .resizable(false)
            .default_width(260.0)
            .show(ctx, |ui| {
                ui.heading("Case");
                ui.label("Closed-closed organ pipe");
                ui.add_space(8.0);

                ui.heading("Run Controls");
                ui.add(egui::Slider::new(&mut self.steps_per_frame, 1..=500).text("steps/frame"));
                ui.add(egui::Slider::new(&mut self.config.cfl, 0.05..=0.95).text("CFL"));
                ui.add(
                    egui::Slider::new(&mut self.config.artificial_viscosity, 0.0..=0.08)
                        .text("artificial viscosity"),
                );
                ui.add(
                    egui::Slider::new(&mut self.config.perturbation_amplitude, 1.0e-5..=1.0e-2)
                        .logarithmic(true)
                        .text("pressure perturbation"),
                );
                ui.add(egui::Slider::new(&mut self.config.cells, 16..=256).text("cells"));
                if ui.button("Apply settings").clicked() {
                    self.reset();
                }

                ui.add_space(12.0);
                ui.heading("Field");
                ui.radio_value(&mut self.selected_field, ScalarField::Pressure, "Pressure");
                ui.radio_value(&mut self.selected_field, ScalarField::Density, "Density");
                ui.radio_value(
                    &mut self.selected_field,
                    ScalarField::Temperature,
                    "Temperature",
                );

                ui.add_space(12.0);
                ui.heading("Status");
                let report = self.run.report();
                ui.label(format!("time: {:.6} s", self.run.time()));
                ui.label(format!(
                    "expected f1: {:.2} Hz",
                    self.run.expected_frequency()
                ));
                if let Some(measured) = self.run.measured_frequency() {
                    let err = ((measured - self.run.expected_frequency())
                        / self.run.expected_frequency())
                    .abs();
                    ui.label(format!("measured f1: {:.2} Hz", measured));
                    ui.label(format!("frequency error: {:.2}%", 100.0 * err));
                } else {
                    ui.label("measured f1: waiting for peaks");
                }
                ui.label(format!("clipped cells: {}", report.clipped_cells));
                ui.label(format!("fallback faces: {}", report.fallback_faces));
                ui.label(format!("snapshots: {}", self.run.history().len()));
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.label(format!(
                    "{} profile at t = {:.6} s ({})",
                    self.selected_field.label(),
                    self.run.latest_snapshot().time,
                    self.selected_field.units()
                ));
                show_plot(ui, ctx, "profile", |buf, size| {
                    draw_profile(buf, size, self.run.latest_snapshot(), self.selected_field)
                });
                ui.separator();

                ui.label(format!(
                    "x-t history: {} ({})",
                    self.selected_field.label(),
                    self.selected_field.units()
                ));
                show_plot(ui, ctx, "history", |buf, size| {
                    draw_xt_history(buf, size, self.run.history(), self.selected_field)
                });
                ui.separator();

                ui.label("Characteristic invariants J+ and J-");
                show_plot(ui, ctx, "characteristics", |buf, size| {
                    draw_characteristics(buf, size, self.run.latest_snapshot())
                });
                ui.separator();

                ui.label("Left-end pressure probe (Pa over time)");
                show_plot(ui, ctx, "probe", |buf, size| {
                    draw_probe_pressure(buf, size, self.run.probe_pressure())
                });
            });
        });
    }
}

fn show_plot<F>(ui: &mut egui::Ui, ctx: &egui::Context, id: &str, draw: F)
where
    F: FnOnce(&mut [u8], (u32, u32)),
{
    let mut rgb = vec![255; (PLOT_W * PLOT_H * 3) as usize];
    draw(&mut rgb, (PLOT_W, PLOT_H));
    let image = egui::ColorImage::from_rgb([PLOT_W as usize, PLOT_H as usize], &rgb);
    let texture = ctx.load_texture(id, image, egui::TextureOptions::LINEAR);
    ui.image((texture.id(), egui::vec2(PLOT_W as f32, PLOT_H as f32)));
}

fn draw_profile(buf: &mut [u8], size: (u32, u32), snapshot: &Snapshot, field: ScalarField) {
    let root = BitMapBackend::with_buffer(buf, size).into_drawing_area();
    let _ = root.fill(&WHITE);
    let values = snapshot.values(field);
    let (min_y, max_y) = padded_range(values);
    let Ok(mut chart) = ChartBuilder::on(&root)
        .margin(12)
        .build_cartesian_2d(x_range(snapshot), min_y..max_y)
    else {
        return;
    };
    let _ = chart.draw_series(LineSeries::new(
        snapshot.x.iter().copied().zip(values.iter().copied()),
        &BLUE,
    ));
}

fn draw_characteristics(buf: &mut [u8], size: (u32, u32), snapshot: &Snapshot) {
    let root = BitMapBackend::with_buffer(buf, size).into_drawing_area();
    let _ = root.fill(&WHITE);
    let mut all = snapshot.c_plus.clone();
    all.extend_from_slice(&snapshot.c_minus);
    let (min_y, max_y) = padded_range(&all);
    let Ok(mut chart) = ChartBuilder::on(&root)
        .margin(12)
        .build_cartesian_2d(x_range(snapshot), min_y..max_y)
    else {
        return;
    };
    let _ = chart.draw_series(LineSeries::new(
        snapshot
            .x
            .iter()
            .copied()
            .zip(snapshot.c_plus.iter().copied()),
        &RED,
    ));
    let _ = chart.draw_series(LineSeries::new(
        snapshot
            .x
            .iter()
            .copied()
            .zip(snapshot.c_minus.iter().copied()),
        &BLUE,
    ));
}

fn draw_probe_pressure(buf: &mut [u8], size: (u32, u32), samples: &[(f64, f64)]) {
    let root = BitMapBackend::with_buffer(buf, size).into_drawing_area();
    let _ = root.fill(&WHITE);
    if samples.len() < 2 {
        return;
    }
    let values: Vec<f64> = samples.iter().map(|(_, p)| *p).collect();
    let (min_y, max_y) = padded_range(&values);
    let min_x = samples.first().map(|(t, _)| *t).unwrap_or(0.0);
    let max_x = samples.last().map(|(t, _)| *t).unwrap_or(1.0);
    let Ok(mut chart) = ChartBuilder::on(&root)
        .margin(12)
        .build_cartesian_2d(min_x..max_x, min_y..max_y)
    else {
        return;
    };
    let _ = chart.draw_series(LineSeries::new(samples.iter().copied(), &GREEN));
}

fn draw_xt_history(buf: &mut [u8], size: (u32, u32), history: &[Snapshot], field: ScalarField) {
    let root = BitMapBackend::with_buffer(buf, size).into_drawing_area();
    let _ = root.fill(&WHITE);
    if history.len() < 2 {
        return;
    }
    let first = &history[0];
    let last = history.last().unwrap();
    let all_values: Vec<f64> = history
        .iter()
        .flat_map(|snapshot| snapshot.values(field).iter().copied())
        .collect();
    let (min_value, max_value) = padded_range(&all_values);
    let Ok(mut chart) = ChartBuilder::on(&root).margin(12).build_cartesian_2d(
        first.x[0]..*first.x.last().unwrap(),
        first.time..last.time.max(first.time + 1.0e-12),
    ) else {
        return;
    };

    let dx = if first.x.len() > 1 {
        first.x[1] - first.x[0]
    } else {
        1.0
    };
    let dt = (last.time - first.time).max(1.0e-12) / history.len() as f64;
    for snapshot in history {
        let values = snapshot.values(field);
        let cells = first
            .x
            .iter()
            .copied()
            .zip(values.iter().copied())
            .map(|(x, value)| {
                let color = heat_color(value, min_value, max_value);
                Rectangle::new(
                    [
                        (x - 0.5 * dx, snapshot.time - 0.5 * dt),
                        (x + 0.5 * dx, snapshot.time + 0.5 * dt),
                    ],
                    color.filled(),
                )
            });
        let _ = chart.draw_series(cells);
    }
}

fn x_range(snapshot: &Snapshot) -> std::ops::Range<f64> {
    let dx = if snapshot.x.len() > 1 {
        snapshot.x[1] - snapshot.x[0]
    } else {
        1.0
    };
    (snapshot.x[0] - 0.5 * dx)..(*snapshot.x.last().unwrap() + 0.5 * dx)
}

fn padded_range(values: &[f64]) -> (f64, f64) {
    let min = values.iter().copied().fold(f64::INFINITY, f64::min);
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !min.is_finite() || !max.is_finite() {
        return (0.0, 1.0);
    }
    let span = (max - min).abs();
    if span <= f64::EPSILON {
        (min - 1.0, max + 1.0)
    } else {
        (min - 0.08 * span, max + 0.08 * span)
    }
}

fn heat_color(value: f64, min: f64, max: f64) -> RGBColor {
    let t = if (max - min).abs() <= f64::EPSILON {
        0.5
    } else {
        ((value - min) / (max - min)).clamp(0.0, 1.0)
    };
    let r = (255.0 * t) as u8;
    let b = (255.0 * (1.0 - t)) as u8;
    let g = (180.0 * (1.0 - (2.0 * t - 1.0).abs())) as u8;
    RGBColor(r, g, b)
}
