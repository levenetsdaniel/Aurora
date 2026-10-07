# VTK electron-flow scene

`aurora_vtk_visualizer` is a standalone 3D visualizer. It builds a reproducible
illustrative source of 5–25 keV electrons at 450 km altitude, propagates them
with the project's dipole field and Boris integrator, and displays their
calculated trajectories and event locations. The default source covers the
entire downward pitch-angle hemisphere (0–90°), uniformly in cos(pitch), so
some electrons enter the 150 km atmosphere and others mirror above it.
With `--aurora`, the source is concentrated in a polar arc and the
precipitating electrons also drive the atmospheric collision/secondary
cascade, excitation grid, stationary-state kinetics, spectral emission grid,
camera images, and VTK glow. The source is a chosen scenario, not an observed
event.

The scene contains Earth, a translucent atmosphere boundary, magnetic field
lines numerically traced along the dipole, color-coded electron trajectories,
cyan injection points, orange precipitation points, and magenta mirror points.
Gold trajectories precipitate, magenta trajectories mirror, pale-blue
trajectories escape, and gray trajectories remain unresolved. Electron energy
is saved in the particle CSV rather than encoded by trajectory color. In aurora
mode the scene also shows separate green 557.7 nm and red 630.0 nm volumes
from the calculated photon-emission grid. VTK coordinates are Earth
radii; input positions stay in metres. Line widths and marker sizes are
enlarged for readability and are not physical particle sizes.

Build and save a PNG (requires VTK 9 and a working OpenGL context):

```sh
cmake -S . -B build
cmake --build build --target aurora_vtk_visualizer
./build/aurora_vtk_visualizer --markers 48 --seed 20260928 --output build/electron_flow.png
```

Add `--interactive` to open the VTK window after saving the PNG. Rotate and
zoom with the mouse. If VTK is unavailable, CMake disables only this optional
target; the physics library and its tests still build. On macOS, offscreen
PNG rendering may still require a GUI/OpenGL-capable user session.

To see the flight as a GIF, install Pillow (`python3 -m pip install Pillow`)
and run the wrapper. The example below uses the same verified IAA archive
required for a glowing scene:

```sh
python3 src/visualization/render_gif.py \
  --gif build/aurora_flight.gif --frames 80 --fps 12 \
  --markers 200 --cascade-markers 64 --threads 8 \
  --seed 20260928 --aurora \
  --iaa-root /absolute/path/to/extracted/archive
```

The command also saves `build/aurora_flight.png` as the final-frame preview
and `<stem>_particles.csv`; a full `--aurora` run also writes the usual emission
and manifest files. Without
`--aurora`/`--iaa-root`, the GIF shows trajectories and outcomes but **no
glow**, unless you explicitly replay a previously calculated emission grid.
For example, if `build/aurora_3d_emission.csv` and its companion manifest are
still present from an earlier run, the following avoids rerunning atmospheric
cascades:

```sh
python3 src/visualization/render_gif.py \
  --gif build/aurora_replay.gif --frames 80 --fps 12 \
  --markers 1000 --seed 20260928 --source-rate 1e20 --threads 8 \
  --emission-csv build/aurora_3d_emission.csv
```

Replay checks that the seed, marker count, source rate and pitch range match
the companion manifest. It reuses the **final** emission grid and does not
recompute the atmosphere. `--frames`, `--fps`, `--width` and `--height` set GIF quality and size;
the default frame size is 960×600. The C++ visualizer can also write numbered
PNG frames directly with `--frames-dir EMPTY_DIR --frames N`.

For a **paired** baseline/anomaly GIF, use `--compare-anomaly`. This samples
the electron source once in the baseline dipole and reuses exactly the same
initial positions, momenta, IDs and weights in both fields. The two panels
use the same simulation clock and camera. With `--aurora`, atmospheric cascades
and emission are computed separately for both panels, with the same source-ID
subsample and shared display brightness reference:

```sh
python3 src/visualization/render_gif.py \
  --gif build/paired_anomaly.gif --frames 80 --fps 12 \
  --markers 200 --cascade-markers 64 --threads 8 \
  --seed 20260928 --aurora --iaa-root /absolute/path/to/extracted/archive \
  --compare-anomaly --anomaly-surface-fraction 0.01 --anomaly-depth-km 500
```

The anomaly is a radially oriented dipole buried below colatitude 0.32 rad;
`--anomaly-surface-fraction` is its axial field at the surface divided by the
baseline field there. `0.01` is a weak model scenario; `0.30` is an intentionally
amplified teaching scenario. `--cascade-markers` picks the **same source IDs**
for atmospheric Monte Carlo in both fields and rescales their weights. The
output includes `<stem>_comparison.csv`, separate particle and emission CSVs,
separate spectral outputs, and a comparison manifest. `--emission-csv` is
deliberately rejected in comparison mode: baseline glow cannot represent an
anomalous-field atmospheric result. Without `--aurora`, the paired GIF still
shows trajectories and outcomes, but neither panel has glow.

Each frame uses recorded simulation times: moving points and growing trails
show electrons from their initial state to their first mirror or 150 km entry.
The final 15% of frames hold the completed scene. Electron sources currently
all start at `t=0`; the simulation stops at the **first** mirror, so return
flight is not shown. The color identifies the eventual outcome. The glow
appears only after precipitation events, scaled by the cumulative arriving
marker weight. This is a **display-only reveal of the final steady-state
emission grid**, not a time-dependent atmospheric excitation or photon-lifetime
calculation. The clock is simulation time, while GIF playback speed is chosen
for viewing.

For an aurora, first obtain the verified IAA archive described in
`data/cross_sections/README.md`. Then run, for example:

```sh
./build/aurora_vtk_visualizer --aurora \
  --iaa-root /absolute/path/to/extracted/archive \
  --markers 1000 --cascade-markers 64 --threads 8 \
  --source-rate 1e20 --seed 20260928 \
  --output build/aurora_3d.png --interactive
```

`--markers` controls the number of rendered Boris trajectories (maximum
5000). `--cascade-markers` is an optional Monte Carlo subsample for the more
expensive atmospheric calculation. If omitted or not smaller than the number
of precipitating markers, *every* precipitating marker gets a cascade. When
subsampling, the selected markers receive equal weights summing to the
precipitating source rate. This speeds up rendering a dense flow but increases
the statistical noise of the emission estimate. `--threads` parallelizes
independent orbit and atmospheric runs. Marker-specific random seeds and
ordered merging make outputs independent of thread scheduling.
Use `--pitch-min-deg` and `--pitch-max-deg` to limit the injected downward
pitch angles; their defaults are 0 and 90 degrees. The source is an illustrative
direction distribution, not a measured flux. A wide-angle run keeps the
original total source rate but only the precipitating fraction reaches the
atmospheric cascade. The orbit calculation stops at the first mirror point:
magenta paths do not yet show the return flight.

During a run, progress is printed to stderr for electron trajectories and,
with `--aurora`, atmospheric cascades. Each stage starts at 0%, reports
approximately every 5%, and ends at 100%. Cascade progress counts the
selected primaries, not all rendered trajectories when `--cascade-markers`
subsamples them.

Aurora mode writes:

- the 3D PNG at `--output`;
- `<stem>_emission.csv`: **unsmoothed** nonzero voxels for both spectral lines,
  with grid geometry and `photons m^-3 s^-1` units;
- `<stem>_spectral/O_5577.csv` and `O_6300.csv`: pixel radiance in
  `photons m^-2 s^-1 sr^-1`, plus a display-only RGB `aurora_preview.ppm`;
- `<stem>_manifest.txt`: seed, marker counts, input locations, source rate,
  neglected-excitation diagnostic, and energy-balance residual.

Every run, including one without `--aurora`, writes `<stem>_particles.csv` with
each marker's initial position, momentum, energy, pitch and gyro phase,
terminal outcome, event time and position, and physical rate weight. These
initial conditions can be reused for paired field-comparison experiments.
If no marker precipitates, aurora mode skips the emission calculation and
saves only the trajectory scene and particle CSV.

The VTK volume uses a photon-conserving spatial blur **only for display** and
compresses the red/green peak contrast for visibility. Thus its color/opacity is not
an absolute photometric prediction. The CSV emission and camera radiances
remain separate numerical outputs. The current collision database extends
most IAA tables to 25 keV using documented *model tails*; those are not
measurements. Only direct atomic-O 557.7/630.0 nm excitation is converted
to visible light. Other excitation channels consume energy but have no
kinetics/line data yet; their rate is reported. The atmosphere is a fixed
illustrative NRLMSIS profile, and the renderer omits atmospheric absorption
and Earth occlusion in the spectral camera. This is an end-to-end synthetic
scene, **not a quantitatively validated aurora forecast**.

The rendering backend accepts `ElectronSceneData`, so another simulation can
provide its own recorded tracks and field lines instead of the demo source.
`SimulationEngine::run(initial, observer)` reports the initial state, each
accepted step and the *located* terminal state. `ElectronSceneSampler`
subsamples these points for display while always retaining the terminal
point. The scene is interactive but not time-animated; path colors encode
the terminal outcome, while the glow comes from the calculated emission grid.
