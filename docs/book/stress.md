# 2. Stress, and who owns the material

Strain says how the body deformed. **Stress** says what internal forces that
deformation produced. The relationship between them is the material's
contract with the world, and it is the one thing in this chapter that Ares
does not own.

## What a stress tensor is

Cut the body along an imaginary plane. The material on one side pulls on the
material on the other. Divide that force by the area of the cut and you have a
**traction** — a force per unit area, a vector.

The traction depends on how you oriented the cut. The **Cauchy stress tensor**
`sigma` is the object that answers, for any orientation, what traction you get:

```text
t = sigma . n
```

where `n` is the unit normal of the cut. So stress is not a number and not a
vector: it is the linear map from a surface orientation to the force per area
that crosses it. Angular momentum balance forces it to be symmetric.

Two readings of it recur, mirroring the strain split of chapter 1.

**Mean stress** is `tr(sigma) / D`, the pressure-like part — the average pull
in all directions. **Deviatoric stress** is what remains, and it is what drives
yielding in metals. The **von Mises** equivalent stress condenses the
deviatoric part into a single number that can be compared against a yield
strength.

```rust
# extern crate ares;
# extern crate aequitas;
use ares::{CauchyStress, SymmetricTensor};
use aequitas::systems::si::quantities::Pressure;

// Uniaxial tension: pulling along x with nothing else.
let tensor = SymmetricTensor::<f64, 2>::from_symmetrised(&[[200.0e6, 0.0], [0.0, 0.0]]);
let stress = CauchyStress::from_tensor(tensor);
assert_eq!(*stress.mean_stress().as_base(), 100.0e6);
```

## Hooke's law, and where it comes from

For a linear elastic isotropic material the relationship is

```text
sigma = lambda tr(eps) I + 2 mu eps
```

Two constants. `mu`, the **shear modulus**, is the resistance to shape change.
`lambda`, the first **Lamé parameter**, has no direct physical reading on its
own but together with `mu` it fixes the resistance to volume change.

Engineers usually quote a different pair: **Young's modulus** `E` (stiffness in
simple tension) and **Poisson's ratio** `nu` (how much a stretched bar thins).
The two pairs describe the same material and convert exactly:

```text
lambda = E nu / ((1 + nu)(1 - 2 nu))
mu     = E / (2 (1 + nu))
```

Ares does not implement that conversion. Proteus does, and Ares consumes the
result:

```rust
# extern crate ares;
# extern crate aequitas;
# extern crate proteus;
use aequitas::systems::si::quantities::{Dimensionless, Pressure};
use ares::{SmallStrain, isotropic_hooke};
use proteus::IsotropicModuli;

let moduli = IsotropicModuli::from_young_poisson(
    Pressure::from_base(200.0e9),      // 200 GPa, a steel-like stiffness
    Dimensionless::from_base(0.3),
).expect("a physically admissible material");

let strain = SmallStrain::<f64, 2>::from_displacement_gradient(&[[1.0e-3, 0.0], [0.0, 0.0]]);
let stress = isotropic_hooke(&moduli, &strain);
// Uniaxial strain, not uniaxial stress: the lateral direction is held, so it
// carries stress too. That is the plane-strain condition of chapter 3.
assert!(*stress.mean_stress().as_base() > 0.0);
```

Notice that the value `200.0e9` appears in *your* code, not in Ares. That is
the ownership boundary doing its job.

## Why the boundary is drawn there

It would be convenient for Ares to carry a small table of common materials.
The reason it does not is that a material property is a different kind of
claim from a balance law.

`sigma = lambda tr(eps) I + 2 mu eps` is a mathematical statement. It is as
true for rubber as for steel; only the constants change. It needs no citation
and has no validity range.

"Steel has `E = 200 GPa`" is a claim about the world. It needs a source, it
depends on the alloy — 304 and 316L differ, and calling either of them "steel"
conflates them — and it is only valid over some temperature range. Those are
Proteus's concerns, and Proteus has the type machinery for them.

There is also a merge-shaped argument. If Ares carried a material table and so
did a fluid solver, the two tables would drift, and they would drift silently
because each is internally consistent. That is not hypothetical: before the
consolidation this stack found steel agreeing across two such tables and
aluminium disagreeing, 70 GPa against 69 GPa, with nothing failing.

## Typed quantities

Stress values are `aequitas` quantities, not bare floats. `Pressure<T>` carries
its dimension in the type, so a length cannot be assigned to a stress and a
stress cannot be added to a strain.

This matters more than it sounds. Traction and pressure have the same physical
dimension and mean different things, and the mistake of treating one as the
other is invisible to a compiler that sees only `f64`. Aequitas carries
semantics markers precisely so that the distinction survives into the type.

## Validity

A material is not just two numbers, it is two numbers in an admissible range.
`IsotropicModuli` refuses to construct outside the positive-definite domain:
`mu > 0` and the bulk modulus `K = lambda + 2 mu / 3 > 0`. Outside it, the
material releases energy when deformed, and the assembled system is no longer
positive definite, and conjugate gradients has no reason to converge on it.

That domain is wider than a naive `lambda >= 0` check would allow, and
deliberately: **auxetic** materials — which get fatter when stretched, with
negative Poisson's ratio — are real, manufacturable, and have negative
`lambda`. Excluding them would be an implementation restriction wearing the
costume of a physical law.
