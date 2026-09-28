# Electron-impact cross sections

## Local IAA dataset (recommended import path)

Download `elmolcs-v0.0.1.zip` from the [IAA Zenodo record](https://zenodo.org/records/8190461)
yourself and keep it outside this repository. The published archive has MD5
`e731acd2c0f203975679ce029b22b187`. After checking the hash, extract it to
a local directory. The resulting layout must contain
`elmolcs/Data/cs/Cross section.txt` and the three species files under
`elmolcs/Data/cs/{O,O2,N2}/`.

```cpp
#include "src/io/IaaCrossSectionLoader.hpp"

// The explicit 2000 eV limit enables documented model tails beyond tables.
auto crossSections = IaaCrossSectionLoader::load("/absolute/path/to/extracted/archive",
                                                 2000.0);
crossSections.requireEnergyCoverage("N2", 1000.0, 2000.0);
```

The loader imports the tabulated `EXCITATION` and `IONIZATION` processes for
O, O₂ and N₂ from the combined file, plus exactly one recommended **integral**
elastic section for each gas from its species file. It ignores `TOTAL`,
`EFFECTIVE`, momentum-transfer elastic sections, and alternative elastic
datasets; adding those to the event rates would double-count collisions.
It validates units, thresholds, nonnegative ordered tables and duplicate IDs.
The local-data integration test runs when `AURORA_IAA_DATA_ROOT` points at the
extraction directory:

```sh
AURORA_IAA_DATA_ROOT=/absolute/path/to/extracted/archive ctest --test-dir build -R iaa_cross_section_loader --output-on-failure
```

Without the second argument, loading is strict: most electronic and
ionization tables stop at 1000 eV; vibrational resonance tables stop at much
lower energies. To run the existing 1–2 keV source, pass `2000.0` explicitly.
Beyond a table's last point, channels ending below 100 eV are assigned zero;
other channels use `sigma(Emax) * Emax / E`, only up to the requested limit.
These are **model assumptions**, not measured or validated cross sections.
Every affected channel records the tail rule in its `source` metadata, and
queries above the chosen limit still fail. The zero-tail assumption for
low-energy channels is particularly coarse and should be revisited before
quantitative predictions.

`ElectronTransportMC` now supports elastic, excitation and ionization events:
elastic events change direction without recoil, excitation removes its stated
threshold, and ionization splits the excess energy between primary and
secondary electrons. The current angular distribution is isotropic in the
laboratory frame and the split fraction is uniform; neither is supplied by
the integral cross-section tables. The secondary is returned in the event but
is not yet transported. Rotational and attachment channels are not imported.
Consequently this completes the **stage-17 transport algorithm and tests**,
but not a fully validated physical auroral cascade: secondary transport and
the energy ledger belong to stage 18, and angular/energy-sharing data need a
better model for quantitative brightness predictions.

The repository includes the two **direct atomic oxygen excitation** channels
`O(3P) -> O(1D)` and `O(3P) -> O(1S)` in
`src/atmosphere/ReferenceCrossSections.hpp`. They are transcribed from Tables 1
and 3 of Laher and Gilmore, *J. Phys. Chem. Ref. Data* **19**, 277 (1990),
[doi:10.1063/1.555872](https://doi.org/10.1063/1.555872). The source tables
use `10^-18 cm^2`; the code converts them to `m^2` with factor `10^-22`.
The thresholds are 1.96 and 4.18 eV, respectively. The last tabulated
point is 200 eV; queries above this energy throw rather than extrapolate.

For N₂ and O₂, obtain the **IST-Lisbon** BOLSIG+ cross-section sets from
[LXCat](https://us.lxcat.net/contributors/) yourself. Preserve the downloaded
file, download date, contributor citation and the references in its header.
LXCat's [redistribution policy](https://us.lxcat.net/instructions/redistribution.php)
does not allow us to bundle its downloaded tables here. The loader in
`src/io/LxCatCrossSectionLoader.hpp` reads one explicitly selected
`EXCITATION` or `IONIZATION` block from a local LXCat/BOLSIG+ text file.
Selection uses the exact reaction line from the file, such as
`N2 -> N2(A)` **only if that exact label occurs in the downloaded set**.
The file must declare `Energy (eV) | Cross Section (m2)`. If the first
tabulated energy is above the stated threshold, the loader adds a zero at
the threshold; it never extrapolates above the last point.

Example after downloading and inspecting a file:

```cpp
#include "src/atmosphere/ReferenceCrossSections.hpp"
#include "src/io/LxCatCrossSectionLoader.hpp"

auto channels = atomicOxygenExcitationChannels();
channels.push_back(LxCatCrossSectionLoader::load(
    "path/to/your/ist-lisbon-n2.txt",
    {"EXACT REACTION LINE FROM FILE", "N2-selected-excitation", "N2",
     CollisionType::Excitation, "EXACT TARGET STATE",
     "IST-Lisbon / LXCat, download date, original channel reference"}));
CrossSectionDatabase db(std::move(channels));
```

Do **not** load `TOTAL` plus its constituent channels into the same database:
that would double count collisions. LXCat's elastic momentum-transfer
cross section is not the integral elastic event cross section; the current
loader intentionally rejects elastic blocks until a compatible scattering
model is defined. The database currently has no representation for
attachment, dissociation, secondary-electron energies or product branching.
Consequently, the included O channels and any individually imported N₂/O₂
channels are **not yet a complete auroral collision model**. In particular,
`totalSigma("O", E)` for the built-in O data sums only the two included
excitation channels, not all O collisions. The O₂ dissociative production of
excited O, which can affect visible emission, is not represented yet.

For independent comparison, see the evaluated reviews of
[O/O₂ (2026)](https://doi.org/10.1063/5.0287254) and
[N₂ (2023)](https://doi.org/10.1063/5.0150618). Before using a downloaded
set, inspect each channel's energy limit against the energies actually
produced by the particle source; `db.requireEnergyCoverage(species, min_eV,
max_eV)` checks this for each loaded channel. The existing electron-source
test samples **1000–2000 eV**, well above the 200 eV limit of the included
atomic-O tables. That configuration is therefore deliberately rejected for
an atmospheric collision calculation until a justified high-energy extension
is supplied. This check establishes numerical coverage, not physical
completeness of the collision model.

`src/atmosphere/ElectronTransportMC.hpp` samples optical depth and integrates
`sum(n_s * sigma_s,k)` with midpoint path cells. Each elastic or ionization
event changes the direction of subsequent segments. Tests check homogeneous
mean free path, channel frequencies, the no-collision probability, the energy
split, direction changes, atmosphere exit and a real local IAA import.
