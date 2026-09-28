# Stage 19: excitation production

`ExcitationModel` converts excitation collision events into a Cartesian field
`P_j` with units `excitations m^-3 s^-1`. Each excited state is identified by
the pair `(gas species, state label)` from `CrossSectionDatabase`. Events from
the primary and all transported secondary electrons are included.

The input `markerRate_s1` is the **number of physical primary electrons per
second represented by one simulated marker**. `IsotropicElectronSource`
assigns this rate to `InitialCondition::weight`, and
`PrecipitationCollector` preserves it as `PrecipitatingElectron::weight`.
Every excitation event in that marker's cascade contributes this rate to its
cell. If the same physical source is sampled with twice as many markers, each
marker must carry half the rate. The model does not multiply by `sigma * v`:
the Monte Carlo collision sampler has already used the cross section.

```cpp
#include "src/emission/ExcitationModel.hpp"

ExcitationModel excitation(crossSections,
    {gridOrigin_m, cellSize_m, cellsX, cellsY, cellsZ});
excitation.accumulate(cascadeResult, precipitatingElectron.weight);
double production_m3_s1 = excitation.productionRate_m3_s1(
    {"O", "O -> O*((2p4)1D)"}, x, y, z);
```

Grid cells are half-open in each dimension: `[lower, upper)`. By default an
excitation outside the grid raises an error without partially changing the
grid. `ExcitationOutOfBounds::Ignore` explicitly drops it and records its
rate in `ignoredEventRate_s1()`. Cascades with `limitReached` are rejected so
that an incomplete collision history is not silently treated as a complete
excitation field.

The model stores weighted event rates in `s^-1` and divides by the cell volume
in `m^3` only when `productionRate_m3_s1()` is requested. The state list is
available through `states()` for the state-kinetics stage. This field is an
excitation-production rate, **not** an excited-state density or photon
emissivity; quenching and radiative decay follow in stages 20 and 21.

Stage 20 is implemented in `ExcitedStateKinetics.hpp`; its coefficient files,
aliases, sources and limitations are documented in `data/kinetics/README.md`.
Stages 21–22 are documented in `RENDERING.md`.
