use onedpipes::{
    DuctConfig, ExternalBoundaryControl, ExternalBoundaryId, Model, ModelBoundary, State,
    TemperatureDependentAir,
};

#[test]
fn public_api_supports_engine_style_external_coupling_loop() {
    let gas = TemperatureDependentAir::new();
    let initial = State::from_primitive(1.2, 0.0, 101_325.0, gas);
    let mut model = Model::new(0.5);
    let pipe = model.add_uniform_duct(
        gas,
        DuctConfig::new(0.6, 16, 3.0e-4),
        initial,
        ModelBoundary::external(0),
        ModelBoundary::open(101_325.0),
    );

    let ports = model.external_ports();
    assert_eq!(ports.len(), 1);
    assert_eq!(ports[0].external_id, ExternalBoundaryId(0));
    assert_eq!(ports[0].pipe_id, pipe);

    model.set_external_boundary_control(
        ports[0].external_id,
        ExternalBoundaryControl::Flow {
            mass_flow_out: 0.0,
            energy_flow_out: 0.0,
        },
    );
    let report = model.step_with_dt(1.0e-7);

    assert_eq!(report.clipped_cells, 0);
    assert_eq!(report.fallback_faces, 0);
    assert_eq!(model.pipe_primitive_cells(pipe).len(), 16);
    assert!(model.pipe_total_mass(pipe) > 0.0);
}
