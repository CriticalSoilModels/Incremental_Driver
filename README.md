# Incremental Driver

A Modern Fortran library for driving and testing geotechnical constitutive models
(Abaqus-style UMATs). The driver applies prescribed stress/strain paths to a UMAT
and records the response — decoupled from any FEA solver so models can be calibrated
and validated independently.

Originally written by Andrzej Niemunis ([soilmodels.com/idriver](https://soilmodels.com/idriver/)).
This repo is an ongoing modernisation: refactoring legacy Fortran into modular,
testable, modern Fortran (2008+) while preserving mathematical correctness.

---

## Using as a library

Add Incremental Driver as a dependency in your project's `fpm.toml`:

```toml
[dependencies]
Incremental_Driver = { git = "https://github.com/CriticalSoilModels/Incremental_Driver" }
```

Then import the public API in your Fortran source:

```fortran
use incremental_driver
```

or selectively:

```fortran
use incremental_driver, only: integrate_step, material_state_t, umat_runner_t
```

### File-free API (recommended)

Build and run a load step entirely in memory — no input files required:

```fortran
use incremental_driver, only: step_config_t, material_state_t, &
                               umat_runner_t, integrate_step

type(umat_runner_t)                  :: runner
type(material_state_t)               :: state
type(step_config_t)                  :: config
type(material_state_t), allocatable  :: results(:)

! Wire up your UMAT
runner%proc   => YOUR_UMAT
runner%cmname =  'YOUR_MODEL'
runner%nstatv =  nstatv
allocate(runner%props(nprops), source=props)

! Set initial state
state%sig    = 0.0_dp
state%eps    = 0.0_dp
state%time   = 0.0_dp
state%dt     = 0.0_dp
state%temp   = 0.0_dp
state%F_start = identity
state%F_end   = identity
allocate(state%statev(nstatv), source=0.0_dp)

! Configure a strain-controlled step
config%load_type     = '*LinearLoad'
config%coord_sys     = '*Cartesian'
config%n_inc         = 10
config%ifstress      = 0          ! strain-controlled in all directions
config%delta_load    = 0.0_dp
config%delta_load(1) = 0.01_dp   ! total axial strain
config%delta_time    = 1.0_dp
! ... (see step_config_t for all fields)

! Run — results(:) contains one material_state_t per increment
call integrate_step(config, state, runner, results)
```

See `example/elastic_umat_demo.f90` for a complete working example.

### File-based API

For running the driver with `test.inp`, `parameters.dat`, and `initial_conditions.dat`
input files (the original iDriver workflow):

```fortran
use incremental_driver, only: run_model

call run_model(YOUR_UMAT)
```

---

## Key types

| Type | Description |
|---|---|
| `material_state_t` | Stress, strain, state variables, time, temperature, deformation gradient |
| `step_config_t` | Load step parameters: type, increments, mixed BCs, coordinate system |
| `umat_runner_t` | Wraps an Abaqus UMAT procedure pointer with props, cmname, nstatv |
| `model_runner_t` | Abstract base — extend to wrap non-UMAT models |
| `step_record_t` | Parsed step config + output write frequency (from `parse_test_file`) |

---

## Building from source

Requires gfortran and [fpm](https://github.com/fortran-lang/fpm).

```bash
# With conda (recommended)
conda env create -f environment.yml
conda activate fpm

fpm build       # build library and stub app
fpm test        # run unit tests
fpm run --example elastic_umat_demo   # run the example
```

---

## License

BSD-style — see [LICENSE](LICENSE).

Please cite Andrzej Niemunis if this code is used in your research.
