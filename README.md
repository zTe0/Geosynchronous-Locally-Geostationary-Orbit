# Locally Geostationary Orbit (LGO) Analysis Tool

A computational tool for designing and analysing **Locally Geostationary Orbits (LGO)** — elliptic satellite orbits whose apogee velocity matches the Earth-surface velocity at a given latitude, maximising dwell time over a target region.

Based on the geometric framework from:

> Razoumny, Y. N. (2018). *Locally Geostationary Orbits: Optimal Geometry of Elliptic Orbit for Earth Coverage.* Journal of Spacecraft and Rockets, 56(4), 1017–1023. https://www.researchgate.net/publication/329948694_Locally_Geostationary_Orbits_Optimal_Geometry_of_Elliptic_Orbit_for_Earth_Coverage

Uses a tool from:

> Luciardello Lecardi M., *"NGSO interference induced on GSO ground station: extension and application of a tool"*, Politecnico di Milano, MSc Thesis, 2024. https://hdl.handle.net/10589/227983


Deployed at: https://zTe0.github.io/Geosynchronous-Locally-Geostationary-Orbit/

---

## Background

An LGO is an elliptic orbit where the apogee velocity equals the velocity of the end-point of the geocentric vector **r_A** rotating with Earth at a given inclination and sub-satellite latitude. This property maximises the time a satellite dwells in the apogee region (visibility zone), making LGOs attractive for high-latitude communication coverage.

The tool computes, as a function of apogee altitude H_A:

- **Circular (first cosmic) velocity** V_c(H_A)
- **Apogee velocity** V_A at fixed perigee altitudes (500 km, 5 000 km, 15 000 km)
- **Earth-surface velocity** V_i at the apogee point for given inclinations (special case: φ = i)
- **LGO apogee velocity and repetition factor m/n** for the special case (φ = i) and the general case (φ ≠ i)
- **Argument of perigee ω** derived from formula (23) of [Razoumny 2018], completing the orbital element set {a, e, i, ω} required by orbit propagators

J2 nodal correction is applied throughout.

---

## Development History

### 1 — MATLAB prototype

The analysis was originally scripted in MATLAB®. Curves replicating Fig. 6 of [Razoumny 2018] were produced using a **search-grid method**: the apogee radius r_A was found by scanning a dense altitude grid and identifying the index where the computed m/n matched each target value.

### 2 — C implementation (`LGO_C2.c` / `LGO_GUI.c`)

The MATLAB results were re-implemented in C for portability and to add a graphical front-end. The C code:

- Reproduces all curves of the MATLAB script (validated by direct plot and table comparison)
- Supports both the **special case** (φ = i, ω = 90°) and the **general case** (φ ≠ i, arbitrary ω)
- Writes a results table to `data/table_results.txt`
- Drives **Gnuplot** for live plotting
- Compiles for both Linux and Windows (MSYS2 / MinGW64)

### 3 — Newton-Raphson solver (`lgo_newton.c`)

A **Newton-Raphson root-finding module** was implemented as a drop-in replacement for the original grid-search `find_hindex()` function, with the following algorithm:

1. **Coarse log-spaced scan** (30 points) over the valid r_A range to bracket the root.
2. **Newton-Raphson step** with central-difference numerical derivative.
3. **Bisection fallback** (Brent-style): if the Newton step escapes the bracket, the interval is halved instead.
4. Convergence when |m/n(r_A) − target| < 10⁻⁹.

Validation confirmed **no significant difference** in computed apogee altitudes or velocities between the grid-search and Newton-Raphson approaches, while the Newton solver offers tighter convergence guarantees and avoids grid-resolution artefacts.

---

## Repository Contents

| File | Description |
|---|---|
| `LGO_GUI.c` | Main C source — computation + GTK3 GUI front-end |
| `lgo_newton.c` | Newton-Raphson solver module (drop-in for grid search) |

---

## Dependencies

| Dependency | Purpose |
|---|---|
| GCC (≥ 9) | Compilation |
| GTK3 (`libgtk-3-dev`) | GUI |
| Gnuplot | Plotting |
| libm | Math library |

---

## Build Instructions

### Linux

```bash
sudo apt install libgtk-3-dev gnuplot
gcc -O2 -Wno-deprecated-declarations \
    $(pkg-config --cflags gtk+-3.0) \
    -o lgo_gui LGO_GUI.c \
    $(pkg-config --libs gtk+-3.0) -lm
```

### Windows (MSYS2 MinGW64 shell)

```bash
pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-pkg-config \
          mingw-w64-x86_64-gtk3 mingw-w64-x86_64-gnuplot

gcc -O2 -Wno-deprecated-declarations \
    $(pkg-config --cflags gtk+-3.0) \
    -o lgo_gui.exe LGO_GUI.c \
    $(pkg-config --libs gtk+-3.0) -lm -mwindows
```

### Newton solver (standalone test)

```bash
gcc -O2 -DLGO_NEWTON_TEST -o lgo_newton_test lgo_newton.c -lm
./lgo_newton_test
```

---

## Run

```bash
./lgo_gui          # Linux
lgo_gui.exe        # Windows
```

---

## GUI Overview

The GTK3 front-end allows full configuration without recompiling:

- **Inclinations I [deg]** — comma-separated list (e.g. `0, 63.4`)
- **m/n targets** — comma-separated list (e.g. `1.0, 1.5, 2.0`)
- **Case selector** — Special (φ = i) or General
- **General case input** — argument of perigee ω [deg] or sub-satellite latitude φ [deg]
- **Run button** — executes the solver and launches a Gnuplot window

Output is written to `data/table_results.txt` with columns:

```
Idx  H[km]  Va[km/s]  m/n  i[deg]  phi[deg]  om[deg]  a[km]  e  Case
```

---

## Key Results

Apogee dwell time (true anomaly window θ ∈ [170°, 190°]) computed via Kepler's equation:

| Satellite | Tundra | Super Tundra (LGO) | Raindrop | Super Raindrop (LGO) |
|---|---|---|---|---|
| Δt [min] | 130.7 | 176 | 131.5 | 175.4 |
| Improvement | — | **+35%** | — | **+33%** |

LGO variants (Super Tundra, Super Raindrop) show a ~33–35% improvement in apogee dwell time relative to the classical Tundra and Draim's Raindrop orbits.

---

## Physical Constants Used

| Symbol | Value | Unit | Meaning |
|---|---|---|---|
| μ | 398 600.4 | km³/s² | Earth gravitational parameter |
| ω_E | 7.2921159 × 10⁻⁵ | s⁻¹ | Earth rotation rate |
| R_mean | 6 371.0 | km | Earth mean radius |
| R_E | 6 378.0 | km | Earth equatorial radius |
| J2 | 0.001 082 7 | — | Second zonal harmonic |

---

## Notation

| Symbol | Meaning |
|---|---|
| H_A | Apogee altitude [km] |
| r_A, r_P | Apogee / perigee geocentric radii [km] |
| V_A | Apogee velocity [km/s] |
| i | Orbital inclination [deg] |
| φ | Sub-satellite latitude at apogee [deg] |
| ω | Argument of perigee [deg] |
| m/n | Repetition factor (m nodal revolutions in n efficient astronomical days) |
| T_ef | Efficient astronomical day (Earth rotation period) |
