# Excited-state kinetics data (stage 20)

The three CSV files in this directory are independent of the C++ code:

- `states.csv`: canonical state, species, total radiative decay rate `A_tot`
  in `s^-1`, and a source URL. The radiative lifetime is `1/A_tot`.
- `quenching.csv`: one or more terms per state and quencher. Values are in
  `m^3 s^-1` and evaluated as
  `k(T) = k0 * (T/T_ref)^beta * exp(-activation_K/T)`.
- `aliases.csv`: exact event-state labels emitted by the cross-section
  database and the canonical state they represent.
- `lines.csv`: line ID, air wavelength in nm, canonical upper state,
  absolute transition rate `A_line` in `s^-1`, and provenance. These are
  **not** branching ratios.

The supplied records cover **O(¹D)** and **O(¹S)** only. These are the upper
states of the red 630.0/636.4 nm and green 557.7 nm oxygen lines. Their
radiative rates are drawn from [Sakai et al. (2014), Table 1](https://doi.org/10.1002/2013JA019197)
for O(¹D), and from the two O(¹S) branches (1.26 and 0.0754 s⁻¹) listed by
[Gronoff et al. (2012)](https://doi.org/10.1029/2011JA017308). Thus the
O(¹S) total of 1.3354 s⁻¹ is a **sum of both radiative branches**, not the
557.7 nm line coefficient alone.

The 557.7339 nm and 630.0304 nm wavelengths and transition rates in
`lines.csv` are from the [NIST Atomic Spectra Database](https://physics.nist.gov/cgi-bin/ASD/lines1.pl?A_out=0&I_scale_type=1&J_out=on&conf_out=on&de=0&en_unit=1&enrg_out=on&forbid_out=1&format=0&g_out=on&intens_out=on&line_out=0&low_w=500&max_low_enrg=&max_str=&max_upp_enrg=&min_accur=&min_intens=&min_str=&order_out=0&output=0&output_type=0&page_size=15&plot_out=0&show_av=2&show_calc_wl=1&show_obs_wl=1&spectra=O+i&submit=Retrieve+Data&term_out=on&tsb_value=0&unc_out=1&unit=1&upp_w=700).
The line catalog rejects a sum of transition rates exceeding the state's
`A_tot`; 630.0 nm uses 5.63e-3 s⁻¹, **not** the O(¹D) total 7.45e-3 s⁻¹.

Neutral quenching terms are transcribed from [Gordiets et al. (1995),
Table V](https://doi.org/10.1109/27.467998), except O(¹D)+O, for which the
9.2 × 10⁻¹³ cm³ s⁻¹ estimate from [Sobral et al. (1993)](https://doi.org/10.1029/92JA01839)
is used. Its published uncertainty is large; another assumption about the
O(¹D) production yield gives a different quenching estimate. The Gordiets
rates were compiled for a low-pressure discharge, so their transfer to all
auroral temperatures and altitudes is a **model choice**, not a validation of
the complete auroral chemistry. Every rate in cm³ s⁻¹ was multiplied by
10⁻⁶ to obtain m³ s⁻¹. Multiple O(¹D)+O₂ channels are intentionally summed.

The current atmosphere supplies only O, O₂ and N₂. All three quenchers must
have explicit entries for each canonical state; a missing entry is an error,
not an implicit zero. The catalog includes aliases for the two built-in
atomic-O cross sections and the corresponding local IAA labels. Other
excitation states (including molecular states) need their own documented
radiative and quenching records before they can be passed to stage 20.

```cpp
#include "src/io/KineticsCsvLoader.hpp"

auto catalog = KineticsCsvLoader::load("data/kinetics");
ExcitedStateKinetics kinetics(catalog, atmosphere);
auto population = kinetics.solve(excitationModel, {"O", "O(1D)"}, x, y, z);
// population.density_m3 is N*; population.effectiveLifetime_s is 1/(A+Q).
```

The calculation is a local **steady-state** solution. It does not yet include
transport of excited neutrals, coupling between excited states, ion chemistry,
or production of O(¹D) by O(¹S) radiative cascade. Those source terms would
need to be added explicitly before quantitative brightness predictions.
