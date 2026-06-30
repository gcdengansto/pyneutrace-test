# pyneutrace: Neutron Ray Tracing Python Package

A pure Python-based Monte Carlo neutron ray-tracing package for neutron instrument design, structural optimization, resolution calculation, and virtual experiments. 

---

## Features

- **High-Performance Tracing** — Monte Carlo neutron ray tracing simulation utilizing fully vectorized operations via `numpy` to trace millions of rays simultaneously without the need for C/C++ compilation.
- **Modular Pipeline** — Easy-to-use component assembly system to seamlessly chain neutron optics, sample environments, and detectors.
- **Built-in Crystallography** — Custom CIF parser natively built into the package to handle real-world crystal structures for powder and single-crystal scattering.
- **Advanced Ray Checkpointing** — Transparently export and merge live neutron states (with metadata validation) to overcome extreme attenuation and accelerate complex downstream simulations.
- **Optimization Engine** — SciPy-backed `InstrumentOptimizer` to tune beamline parameters (e.g., guide focus, chopper phase) for maximum flux/resolution.

---

## Directory Structure

```text
pyneutrace/
├── src/
│   └── pyneutrace/ 
│       ├── __init__.py 
│       ├── constants.py       # Physical constants and standard D-spacings
│       ├── utils.py           # Core utilities, including the custom CIF parser
│       ├── instrument/        # InstrumentAssemble pipeline manager
│       ├── components/        # Optics, samples, choppers, and detectors
│       ├── simulation/        # InstrumentOptimizer
│       └── visualization/     # InstrumentSimPlot for statistical visualization
├── test/                      # Comprehensive test suite covering all components
├── pyproject.toml             # Package build configuration
└── README.md
```

---

## Quick Start

1. **Installation:** Run `pip install .` in the project root path for the core engine, or `pip install ".[gui]"` to include 3D visualization dependencies.
2. **Build your instrument:** Instantiate components and chain them using the `InstrumentAssemble` class.
3. **Run the simulation:** Execute `pipeline.run(num_rays=1_000_000)` to trace the beam.
4. **Visualize 3D Geometry:** Call `pipeline.visualize_3d()` to render the physical instrument layout via PyVista.
5. **Analyze Results:** Use the `InstrumentSimPlot` class to visualize transmission drop-off and spatial/energy beam profiles at any stage.

---

## User Interface

### Available Components (`pyneutrace.components`)

| Component Class | Description |
|---|---|
| **NeutronSource** | Continuous neutron source (reactor) with a dual-temperature Maxwellian flux spectrum. |
| **Moderator** | Pulsed neutron source (spallation) with wavelength-dependent Ikeda-Carpenter pulse shaping. |
| **NeutronGuide** | Straight or tapered supermirror guide (supports m-value reflectivity decay). |
| **CurvedGuide** | Horizontally curved supermirror bender for filtering fast neutrons and gammas. |
| **Elliptic / Parabolic Guides**| Focusing optic guides (Single Elliptic, Double Elliptic, Double Parabolic). |
| **DiskChopper** | Standard rotating disk chopper for pulse shaping and order sorting. |
| **Straight/Curved FermiChopper**| High-speed rotating slit packages for precise energy selection. |
| **VelSelector** | Mechanical velocity selector (helical turbine) for monochromatic beam selection. |
| **Monochromator / Analyzer** | Single crystal (e.g., PG, Si, Ge) arrays to select specific neutron energies via Bragg's law. |
| **PowderSample** | Realistic sample environment supporting CIF loading and Debye-Scherrer cone diffraction. |
| **VRodSample** | Cylindrical sample demonstrating isotropic elastic scattering. |
| **Monitor / CylindMonitor** | Flat and curved Position Sensitive Detectors (PSD) for capturing 2D spatial and 1D energy histograms. |
| **Soller / RadCollimator** | Linear and radial collimation packages to define beam divergence. |
| **ExportRays / ImportRays** | I/O components to checkpoint live beams to disk and merge multiple runs. |
| **Spacer** | A virtual component used to advance the beam a specific distance in empty space. |

### Core Pipeline (`pyneutrace.instrument`)

| Class | Description |
|---|---|
| **InstrumentAssemble** | The central pipeline manager. Handles the addition of components, calculates global 4x4 coordinate transformations, and executes the sequential Monte Carlo propagation (`propagate_TOF`). |

### Visualization (`pyneutrace.visualization`)

| Class | Description |
|---|---|
| **InstrumentSimPlot** | Generates Matplotlib dashboards summarizing transmission efficiency and "Virtual Monitor" beam profiles without requiring heavy memory storage during runs. |

---

## Physics & Architecture

### Coordinate System

```text
  X  →  Horizontal (perpendicular to the beam propagating direction)
  Y  →  Vertical (up)
  Z  →  Beam propagating direction 
```

**Local vs. Global Frames**  
Each component in the instrument possesses its own local coordinate system. A component clearly defines its **entry window** (usually at $Z = 0$) and its **exit window** (usually at $Z = L$). 

`InstrumentAssemble` builds the instrument as a **seamless continuous pipeline**. This means all components, from the source to the detector, are connected back-to-back. The exit window of an upstream component becomes the exact mathematical entry window of the subsequent component. 

To achieve this:
1. All surviving neutrons exiting a component are mathematically translated (and rotated, if the component bends the beam like a `CurvedGuide` or `Monochromator`) so that their positions and velocities are expressed perfectly relative to the **center of the exit window**.
2. This allows the next component to simply assume the incoming beam is crossing its own $Z=0$ plane. 
3. If physical empty space is required between two optical elements, a `Spacer` component must be explicitly inserted into the pipeline.

**The Sample Exception**  
The only exception to the strict Entry/Exit window rule is the Sample environment. Because neutrons scattered from a sample diverge outward into the entire 3D sphere ($4\pi$ steradians), it is impossible to define a single planar "exit window". 

Instead, neutron coordinates scattered away from a sample use the **sample center** as their origin ($0,0,0$), with the incident primary beam defining the $+Z$ direction. Downstream components (such as a `RadCollimator` or `CylindMonitor`) are explicitly designed to account for this geometry, often utilizing cylindrical or spherical coordinates to capture the flying neutrons.

### Bragg Scattering

For the `Monochromator`, `Analyzer`, and `PowderSample`, the Bragg condition is applied to filter and scatter neutrons at specific wavelengths. The momentum transfer vector $\vec{Q}$ is calculated for reciprocal lattice vectors $\vec{G}$ and the incident wavevector $\vec{k}_{in}$, ensuring that only neutrons satisfying $|\vec{k}_{out}| = |\vec{k}_{in}|$ and $\vec{Q} = \vec{G}$ are scattered with high probability according to their structure factors ($|F|^2$).

---

## Known Limitations

- **Gravity:** The impact of gravity ($\vec{g}$) has not yet been taken into account. Neutrons currently travel in perfectly straight lines between component boundaries.
- **Detector Efficiency:** Monitors currently act as perfect 100% absorbing planes. Wavelength-dependent gas detection efficiency ($1 - \exp(-c\lambda)$) and finite depth tracking are not yet implemented.
- **Inelastic Scattering:** Sample kernels for inelastic energy transfer ($S(Q, \omega)$) to simulate phonons and magnons are not yet supported.
