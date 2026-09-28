# VTK electron-flow scene

`aurora_vtk_visualizer` is a standalone 3D visualizer. It builds a reproducible
illustrative source of 5–25 keV electrons at 450 km altitude, propagates them
with the project's dipole field and Boris integrator, and displays their
calculated trajectories and actual event locations at the 150 km atmospheric
boundary. With `--aurora`, the source is concentrated in a polar arc and the
precipitating electrons also drive the atmospheric collision/secondary
cascade, excitation grid, stationary-state kinetics, spectral emission grid,
camera images, and VTK glow. The source is a chosen scenario, not an observed
event.

The scene contains Earth, a translucent atmosphere boundary, magnetic field
lines numerically traced along the dipole, color-coded electron trajectories,
cyan injection points, orange precipitation points, and an energy scale in
keV. In aurora mode it also shows separate green 557.7 nm and red 630.0 nm
volumes from the calculated photon-emission grid. VTK coordinates are Earth
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

Aurora mode writes:

- the 3D PNG at `--output`;
- `<stem>_emission.csv`: **unsmoothed** nonzero voxels for both spectral lines,
  with grid geometry and `photons m^-3 s^-1` units;
- `<stem>_spectral/O_5577.csv` and `O_6300.csv`: pixel radiance in
  `photons m^-2 s^-1 sr^-1`, plus a display-only RGB `aurora_preview.ppm`;
- `<stem>_manifest.txt`: seed, marker counts, input locations, source rate,
  neglected-excitation diagnostic, and energy-balance residual.

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
point. The scene is interactive but not time-animated; the path colors encode
electron energy, while the glow comes from the calculated emission grid.
