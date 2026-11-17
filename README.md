# C-Digital-Beamforming

A C++/CUDA research prototype for digital beamforming, developed in the context of the [Project 8 neutrino mass experiment](https://www.project8.org/home).

The original Phase III concept envisaged an antenna array for detecting cyclotron radiation from electrons. That detector concept has since been abandoned. This prototype was used to investigate reconstruction techniques suitable for online triggering, studying both computational performance and trigger discrimination through simulated signals and noise.

The code combines a simplified electron-event and antenna-response model with frequency-domain beamforming on the GPU. It includes spatial reconstruction, signal-to-noise studies, and threshold-trigger evaluation using ROC curves.

The prototype was discontinued because the studies required substantially more modelling complexity than it provided, leading to the development of a separate, more complete simulation model (see [CRESana](https://github.com/MCFlowMace/CRESana)).

**Status:** Unmaintained research prototype representing a 2020 development state. Targets the CUDA 11.1 toolchain.

## Dependencies and tested environment

### C++/CUDA

The project targets the following toolchain:

| Dependency | Version / configuration |
|---|---|
| CUDA Toolkit | 11.1.1 |
| GCC/G++ | 9 |
| OpenMP | GCC OpenMP support (`-fopenmp`, linked against `libgomp`) |
| Armadillo | Tested with 12.6.7; no explicit version constraint |
| Build tools | GNU Make and binutils |

GPU reconstruction requires a compatible NVIDIA GPU and driver. The Makefile defaults to `ARCH=sm_75` for the Tesla T4. The original development GPU was a Titan V (`sm_70`).

### Python

| Dependency | Version information |
|---|---|
| Python | tested with Python 3.8 |
| NumPy | tested with 1.19.5 |
| Cython | tested with 0.29.36 |
| setuptools | tested with 59.8.0 |
| Matplotlib | No explicit version constraint |
| SciPy | ≥ 1.6 |
| Seaborn | No explicit version constraint |
| Statsmodels | No explicit version constraint |

The listed extension versions form a tested combination, rather than a declaration that other versions are supported or unsupported.

### Runtime verification

C++ builds, GPU execution, the Cython extension, and plotting experiments were exercised on Ubuntu 24.04.4 with a Tesla T4, NVIDIA driver 580.82.07, CUDA Toolkit 11.1.1, GCC/G++ 9, and Armadillo 12.6.7. The Python extension used the versions listed above.

### Development container

[`.devcontainer/`](.devcontainer/) provides an Ubuntu 20.04 environment based on:

```text
nvidia/cuda:11.1.1-devel-ubuntu20.04
```

The image installs compiler, Armadillo, and Python dependencies from the Ubuntu repositories; those packages are not individually pinned. Its default SciPy may require updating for the ROC plotting script. The Dev Container configuration requests GPU access.

## Structure

| Path | Purpose |
|---|---|
| `include/beamforming/`, `lib/beamforming/` | Simulation model, antenna response, reconstruction, and trigger evaluation |
| `lib/beamforming/reconstruction_gpu.cu` | CUDA beamforming kernels |
| `lib/pybeamforming/` | Cython bindings |
| `lib/analysis/` | Plotting utilities |
| `test/` | Numerical experiments and plotting scripts |

## Build

```bash
make GCC=g++-9 NVCC="nvcc -ccbin g++-9"
```

The default target builds `build/lib/libbeamforming.a`. Example targets link against this library and place executables in `build/bin/`:

| Make target | Executable |
|---|---|
| `reconstruction` | `test_rec` |
| `response_map` | `response_map` |
| `measure-snr` | `measure-snr` |
| `threshold-trigger` | `threshold-trigger` |
| `roc-evaluation` | `roc-evaluation` |
| `event_generator` | `generate_events` |
| `e_loss` | `e_loss` |

For example:

```bash
make reconstruction response_map measure-snr threshold-trigger \
    roc-evaluation event_generator e_loss \
    GCC=g++-9 NVCC="nvcc -ccbin g++-9"
```

## Experiments

Commands below assume the repository root as the working directory. Plots are saved in the current directory.

These are exploratory numerical experiments with visual inspection of their outputs, rather than an automated regression suite. The reconstruction, response-map, SNR, and threshold-trigger experiments require GPU execution.

### Spatial reconstruction

```bash
python3 test/reconstruction/test_reconstruction.py
```

Runs `test_rec` and saves `beamforming_rec_test.pdf`, showing the reconstruction at the frequency bin selected by the beamformer.

The executable has no command-line configuration. Parameters are defined in `test/reconstruction/test_rec.cpp` and the defaults in `lib/beamforming/simulation.cpp`.

### Spatial response map

```bash
./build/bin/response_map 20 500 10000 > response.out
python3 test/response_map/plot_result.py response.out 5
```

Evaluates reconstruction response as the simulated source is moved across the grid. Saves the normalized logarithmic response as `response_map.pdf`.

Arguments: `grid_size n_samples snr`. The current implementation parses `snr` as an integer.

### SNR scaling

Sample-count sweep:

```bash
./build/bin/measure-snr 20 100 1001 30 31 > snr_samples.dat
python3 test/measure-snr/plot_snr.py snr_samples.dat True
```

Antenna-count sweep:

```bash
./build/bin/measure-snr 20 500 501 15 41 > snr_antennas.dat
python3 test/measure-snr/plot_snr.py snr_antennas.dat False
```

Each configuration is evaluated over 50 repetitions, comparing signal-plus-noise reconstruction with a noise-only reconstruction. Plots are saved as `<input>.pdf`.

Arguments:

```text
grid_size n_min n_max N_min N_max
```

Upper bounds are exclusive. Sample and antenna increments are 100 and 5 respectively.

Output columns:

```text
N  n_samples  snr  signal_err  noise_err  dsnr
```

### Trigger discrimination

```bash
python3 test/threshold-trigger/plot_roc.py
```

Sweeps classification thresholds and compares ROC curves for grid sizes 25, 50, and 100 at input SNR values 0.1 and 0.5. Saves `roc_1.pdf` and `roc_5.pdf`. Sweep parameters are defined in the script.

For direct invocation:

```bash
./build/bin/threshold-trigger \
    20 500 0.5 1000 42 3.0 0.0 26000000000 > roc.dat
```

Arguments:

```text
grid_size n_samples snr total_packets seed r phi f0
```

`f0` is in Hz. A negative `r` selects random event generation; otherwise the source is configured manually. `total_packets` must be a positive multiple of the fixed batch size, 100.

Output columns:

```text
false_positives  negative_labels  true_positives  positive_labels
```

### Event generation and energy-loss sampling

These examples do not execute GPU code.

```bash
./build/bin/generate_events 0.05 0.1 24.6 26.2 5.0

python3 test/event_generation/plot_event.py \
    event0_freq.dat 0.05 0.1 24.6 26.2 event0_pos.dat 5.0
```

Generates an event and plots its frequency evolution and trajectory as `event0_spec.pdf` and `event0_pos.png`. Generator arguments are `t_min t_max f_min f_max R`, with frequency limits in GHz and radius in centimetres.

```bash
./build/bin/e_loss > e_loss.dat
python3 test/event_generation/plot_e_loss.py
```

Generates 100000 energy-loss samples and compares their distribution with the reference model in `e_loss.pdf`.

```bash
python3 test/event_generation/scatter_dis.py
```

Runs a separate Python sampling experiment and saves `scatter_dis.pdf`.

```bash
./build/bin/roc-evaluation
```

Prints confusion-matrix results for hardcoded classification examples.

## Python bindings

```bash
make pybeamforming PYTHON=python3

PYTHONPATH="$PWD/build/lib" python3 -c "import pybeamforming"
```

`lib/pybeamforming/setup.py` invokes `nvcc` directly, independently of the Makefile's `NVCC` variable. Ensure that invocation selects the intended CUDA toolkit and GCC 9 host compiler. The setup script also contains a fixed CUDA library path, `/usr/local/cuda-11.1/lib64`.

The extension exposes `PySimulation_Settings`, `PyBeamformer`, and double-precision variants. Returned reconstruction arrays have axes:

```text
[packet, frequency, y, x]
```

Use `full_frequency=False` with the existing simulated-data path. The plotting experiments above invoke the C++ executables rather than the bindings.

## Implementation notes

- Timing diagnostics are emitted to stderr.
- `Reconstruction_GPU::get_img()` returns spatial axes in the opposite order to the Cython array interface. The reconstruction plotting script accounts for this with a transpose.
- The container uses Matplotlib's `Agg` backend; scripts generally save figures rather than display them.
- Older shell helpers under `test/` reference outdated executable paths; use the commands above.
- `plot_event.py` labels frequency as GHz but currently converts angular frequency to Hz. The GHz conversion requires division by `2*pi*1e9`.