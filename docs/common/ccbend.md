# CCBEND in AT

The three CCBEND PassMethods port elegant's straight-axis combined-function
bend, including its actual fringe maps and reference-trajectory optimization.

| PassMethod | Radiation loss | Quantum excitation |
| --- | --- | --- |
| `CCBendPass` | No | No |
| `CCBendRadPass` | Yes | Diffusion matrix only |
| `CCBendQuantPass` | Yes | Gaussian ISR |

## Usage

```matlab
B=atccbend('CB',1.2,-0.12,-0.8,'PolynomB',[0 -0.8 1.5], ...
    'FringeModel',1,'fringeInt1',[0 .001 0 0 0 0 .001 .0002], ...
    'fringeInt2',zeros(1,8));
rout=CCBendPass(B,zeros(6,1));
B.PassMethod='CCBendRadPass';B.Energy=3e9;
radiating=CCBendRadPass(B,zeros(6,1));
```

```python
import at
import numpy as np

bend = at.CCBend('CB', 1.2, -0.12, -0.8,
                PolynomB=[0, -0.8, 1.5], FringeModel=1,
                fringeInt1=[0, .001, 0, 0, 0, 0, .001, .0002], fringeInt2=[0]*8)
rout = at.element_pass(bend, np.zeros((6, 1)))
bend.PassMethod = 'CCBendQuantPass'
excited = at.element_pass(bend, np.zeros((6, 1000)), energy=3e9, seed=42)
```

`elegant2at` reads saved, flattened elegant lattices. The Python elegant importer
and exporter also support CCBEND. Unsupported nonzero options fail explicitly.

## Parameters and conventions

| AT | elegant | Default / meaning |
| --- | --- | --- |
| `Length` | `L` | Nominal arc length; also used for lattice length |
| `BendingAngle` | `ANGLE` | Signed, nonzero; absolute value less than pi |
| `PolynomB[n]` | `Kn/n!` | MATLAB index is n+1; nominal dipole coefficient is zero |
| `NumIntSteps` | `N_SLICES` | 4 |
| `IntegrationOrder` | `INTEGRATION_ORDER` | 4; accepts 2, 4, 6 |
| `FringeModel` | `FRINGEMODEL` | 0; -1=off, 0=original, 1=Lindberg integral map |
| `EdgeOrder` | `EDGE_ORDER` | 3; model 1 uses its second-order map for values ≥2 |
| `fringeInt1/2` | `FRINGE1/2K0`, `K2`, `K4`, `K5`, `K6`, `K7`, `I0`, `I1` | Eight values per end, zero by default |
| `FullGap`, `FringeInt1/2` | `2*HGAP`, `FINT1/2` | Scalars for model 0, zero by default |
| `OptimizeFSE`, `OptimizeDX` | `OPTIMIZE_FSE`, `OPTIMIZE_DX` | Both 1 |
| `CompensateKN` | `COMPENSATE_KN` | 0 |
| `ReferenceCorrection` | `REFERENCE_CORRECTION` | 1=path, 2=transverse, 3=both |
| `EdgeFlip`, `Yaw` | `EDGE_FLIP`, `YAW` | 0 |
| `FSE`, `FSEDipole`, `FSEQuadrupole` | `FSE`, `FSE_DIPOLE`, `FSE_QUADRUPOLE` | 0 |
| `Energy` / tracking energy | Beam energy | Total energy in eV; required for radiation |
| `GaussianCutoff` | Global `srGaussianLimit` | 3 |
| `QuantumSingleParticle` | `ISR1PART` | 1; 0 suppresses noise in one-particle calls |

Field names are case-sensitive: `fringeInt1/2` are arrays, whereas
`FringeInt1/2` are scalar FINT values. The intrinsic dipole field is added
internally; leave `PolynomB[0]` zero for the nominal element.

AT uses `[x,px,y,py,delta,ct]`. The adapter converts to elegant's
`[x,x',y,y',s,delta]` and back, returning ct=s−Length. Reference correction
does not alter `Length`. Both optimization flags off disable reference correction.
Optimization uses elegant's residuals with a bounded local Nelder-Mead solver,
so offsets need not be bitwise identical to elegant's `simplexMin` results.
Reset tracking caches after editing elements. Do not split or merge CCBEND.

## Radiation and scope

Radiation is applied once per positive-length slice using elegant's update.
Reference optimization excludes radiation and errors. `CCBendQuantPass` uses
Gaussian ISR with default ±3σ rejection, rather than discrete photon emission.
AT and elegant random sequences differ. Small physical-constant differences
can remain between AT and an elegant executable.

Use `CCBendRadPass` for optics and envelope calculations. The diffusion kernel
transports each slice's ISR impulse through the complete map and returns
untruncated theoretical diffusion. Numerical sensitivities are computed only
when diffusion is requested; ordinary tracking has no extra cost.

The port supports negative angles, skew polynomials, AT alignment/apertures and
both elegant fringe models. It excludes spin, SDDS multipole-error files,
partial/backward tracking, discrete photons and elegant's 3D misalignment.
Use AT `R1/R2/T1/T2` for alignment. Original fringe formulas retain elegant's
small-parameter cancellation. ISR without radiation loss is not a supported mode.

## Tests

```sh
python -m pytest pyat/test/test_ccbend.py
```

```matlab
assertSuccess(runtests('atmat/attests/testCCBend.m'));
```

Use an AT-enabled Python/MATLAB environment with compiled PassMethods.
Both suites compare all R/T coefficients with the saved elegant `ccbend1` Q5
output, using the example's nine-point grid and full coordinate conversion.
The shared `pyat/test/ccbend1_reference.json` contains the full-precision
reference and provenance; tests need neither SDDS nor a local elegant installation.
Separate radiation checks cover loss, diffusion and sampled ISR variance.

## Attribution

The fringe and plane-switch functions derive from elegant `src/ccbend.c`.
See `atintegrators/elegant_ccbend_fringe.h` for copyright and adaptation notices,
and `atintegrators/elegant.LICENSE` for the original license. The license is
also included in the Python integrator package.
