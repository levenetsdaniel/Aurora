# Stages 21–22: spectral emission and camera

`VolumeEmissionGrid::fromKinetics` turns each stage-20 population into
`j_line = A_line * N*` in photons m^-3 s^-1. It sums production from different
event-state aliases that resolve to the same canonical state. A state with
positive production but missing kinetic or line data is an error. The two
provided outputs are O_5577 and O_6300. Their transition rates are loaded
from `data/kinetics/lines.csv`, not compiled into the model.

`OpticalRenderer` treats each voxel as constant, clips each pixel-center ray
to the grid box, splits the ray at voxel boundaries, and integrates exactly
within that piecewise-constant field. Spectral-image pixels are
`I = integral(j ds)/(4*pi)` in photons m^-2 s^-1 sr^-1. The model assumes
isotropic, optically thin emission with no absorption, scattering, or
occlusion by Earth. The first image row is the top of the camera view.

```cpp
#include "src/emission/OpticalRenderer.hpp"
#include "src/io/KineticsCsvLoader.hpp"
#include "src/io/SpectralLineCsvLoader.hpp"

auto states = KineticsCsvLoader::load("data/kinetics");
auto lines = SpectralLineCsvLoader::load("data/kinetics/lines.csv", states);
ExcitedStateKinetics kinetics(states, atmosphere);
auto volume = VolumeEmissionGrid::fromKinetics(excitation, kinetics, states, lines);
OpticalRenderer camera({cameraPosition_m, forward, up, 60.0, 640, 480});
auto images = camera.render(volume);
OpticalRenderer::saveSpectralCsv(images, "output");
auto preview = OpticalRenderer::makeRgbPreview(images[1], images[0], exposure);
preview.savePpm("output/aurora_preview.ppm");
```

The camera orthonormalizes `forward` and `up`; parallel or invalid vectors
are rejected. `exposure` is a user-selected radiance scale in the same units
as the spectral image. The RGB preview uses a nonlinear display transform;
its 8-bit color values are **not physical radiance**. The CSV files retain
the separate physical spectral data. Stage 23 must still connect these
modules to a full simulation run.
