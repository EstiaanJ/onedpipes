mod closed_end;
mod junction;
mod open_end;

pub use closed_end::ClosedEnd;
pub use junction::{JunctionPort, JunctionSolution, MultiPipeJunction, PortFlow};
pub use open_end::OpenEnd;

use crate::{gas_properties::GasProperties, state::State};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DuctEnd {
    Left,
    Right,
}

pub trait BoundaryCondition<G: GasProperties> {
    fn ghost_state(&self, interior: State, end: DuctEnd, gas: G) -> State;
}
