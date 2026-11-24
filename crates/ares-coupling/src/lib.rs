#![doc = include_str!("../README.md")]
#![forbid(unsafe_code)]
#![deny(missing_docs)]

/// The coupling interface: exchange ordering, conformity, and interface work.
mod interface;
/// The structural partition Harmonia drives.
mod partition;

pub use interface::{InvalidInterface, StructuralInterface};
pub use partition::{PartitionError, StructuralPartition};
