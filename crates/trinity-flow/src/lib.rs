//! M2 流分解输出层: 量化 de Bruijn 图上的最重路径剥离。

pub mod engine;
pub mod graph;
pub mod paths;
pub mod peel;

pub use engine::{run_flow_component, FlowParams};
