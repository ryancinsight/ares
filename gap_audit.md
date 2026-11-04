# Ares correctness audit

## 2026-11-04: traction coordinate identity

The public boundary constructor validates node indices and facet measure, but
formerly retained only the facets. On unchanged main, a later assembly call
with one node for facet `[0, 1]` passes its replacement-based shape check and
panics during coordinate gathering. Equal-length coincident replacement nodes
instead produce `Ok` and no added traction. Both witnesses execute in f32/f64.

The repair retains the exact immutable coordinate slice and removes the second
coordinate authority from assembly. Continued-use compile-fail examples reject
replacement input, mutation and lifetime escape. Positive-D const assertions
cover both public paths; D=1, empty boundaries, exact 2-D/3-D loads, shared-node
addition and untouched mismatched output have generic native/wrapper controls.
Existing force, moment, energy and coupling oracles retain their assertions.

The registry-normalized baseline exposes a separate E0277 material-provider
identity mismatch (Eunomia 0.12 versus published Proteus's 0.10). Successful
library packaging alone does not qualify those tests or release closure.
The [board](backlog.md) retains final-source gates and provider/first-publication
prerequisites; no source, hosted, consumer or publication acceptance is inferred.
