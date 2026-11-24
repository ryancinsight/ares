//! Compare 3D symmetric projection against the full matrix midpoint oracle.
use ares::SymmetricTensor;
use std::{
    hint::black_box,
    io::{self, Write},
    time::{Duration, Instant},
};

type Gradient = [[f64; 3]; 3];
const BATCH: u32 = 65_536;
const PASSES: u32 = 16;

fn elapsed<R>(gradients: &[Gradient], project: impl Fn(&Gradient) -> R) -> Duration {
    let start = Instant::now();
    for _ in 0..PASSES {
        for gradient in black_box(gradients) {
            black_box(project(gradient));
        }
    }
    start.elapsed()
}

fn main() -> io::Result<()> {
    let batch: Vec<Gradient> = (0..BATCH)
        .map(|index| {
            let offsets = [[1_u32, 2, 3], [4, 5, 6], [7, 8, 9]];
            offsets.map(|row| row.map(|offset| f64::from(index + offset).sin() * 1e-5))
        })
        .collect();
    let mut samples = [Duration::ZERO; 15];
    let mut reference = samples;
    for (production, oracle) in samples.iter_mut().zip(&mut reference) {
        *production = elapsed(&batch, SymmetricTensor::from_symmetrised);
        *oracle = elapsed(&batch, |gradient| -> Gradient {
            std::array::from_fn(|i| {
                std::array::from_fn(|j| gradient[i][j].midpoint(gradient[j][i]))
            })
        });
    }
    for (label, values) in [
        ("projection", &mut samples),
        ("full midpoint oracle", &mut reference),
    ] {
        values.sort_unstable();
        let ns = values[values.len() / 2].as_secs_f64() * 1e9 / f64::from(BATCH * PASSES);
        writeln!(io::stdout().lock(), "{label} {ns:.3} ns/projection")?;
    }
    Ok(())
}
