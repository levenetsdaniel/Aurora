# Atmosphere reference profile

`nrlmsis20_2020-03-20_70N_0E_demo.csv` is a fixed, **illustrative** vertical profile from the original [NRLMSIS 2.0 release](https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSIS2.0/). It is model output, not an observation of an aurora. The solar and geomagnetic indices were chosen as inputs, not fetched as historical values.

The unmodified NRL Fortran `MSISCALC` interface was evaluated with double-precision internals (`-DDBLE`), default model switches, day 80, UT 12:00, geodetic latitude 70°, longitude 0°, previous-day F10.7 = 70, centered 81-day F10.7 = 70, and daily Ap = 4. Each altitude from 70 to 300 km was evaluated in 5 km increments. The call returned O, O₂, N₂ number densities in m⁻³ and neutral temperature in K. The CSV header names are also the units expected by `AtmosphereCsvLoader`. The same compiler build reproduced the package's `msis2.0_test_ref_dp.txt` byte for byte.

The table uses **geodetic** altitude. The current simulation uses the spherical approximation `|r| − Rₑ`; these coordinates are not exactly identical. This fixed 1D profile also does not change with latitude, longitude, or time when queried by `AtmosphereModel`. Do not use it as a universal or event-specific atmosphere.

The NRL software package is not included in this repository. Consult its academic/non-commercial license before using the model software or this derived profile in another context.

In C++:

```cpp
const AtmosphereModel atmosphere = AtmosphereCsvLoader::load(
    "data/atmosphere/nrlmsis20_2020-03-20_70N_0E_demo.csv");
const AtmosphereSample at150km = atmosphere.sample(150000.0, 0.0, 0.0, 0.0);
```
