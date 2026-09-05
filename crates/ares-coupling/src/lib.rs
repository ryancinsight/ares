//! The Harmonia coupling partition for the Ares solid momentum balance.
//!
//! Atlas ADR 0059 splits Phase 0 fluid-structure interaction so that `CFDrs`
//! produces the interface traction from its own flow state and a structural
//! solver consumes it. Neither depends on the other; this is the structural
//! half, joined to the fluid half by a Harmonia driver.
//!
//! # What this crate is and is not
//!
//! It owns no physics and no coupling mechanics. Ares balances, Athena solves,
//! and Harmonia drives the partitions; this crate is the adapter that lets
//! Harmonia see a structural solve as a `Partition`.
//!
//! It carries no dependency on `CFDrs`, `kwavers`, or any other balance domain,
//! and cannot: atlas ADR 0055 R7 forbids a direct balance-to-balance edge and
//! routes coupling through Harmonia, which is exactly the shape here.
//!
//! # Phase 0 is one-way
//!
//! Traction flows fluid to solid; displacement is exported but does not move
//! the fluid mesh, because `CFDrs` has no ALE. The partition is therefore
//! quasi-static: `advance` solves an equilibrium problem for the traction it
//! was handed and does not integrate in time.

#![forbid(unsafe_code)]
#![deny(missing_docs)]

/// The coupling interface: exchange ordering, conformity, and interface work.
mod interface;
/// The structural partition Harmonia drives.
mod partition;

pub use interface::{InvalidInterface, StructuralInterface};
pub use partition::{PartitionError, StructuralPartition};
