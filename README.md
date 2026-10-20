# Plasmaless Tokamak Circuits

A dynamical circuit model of the current response of ITER coils and vessel conductive elements to applied power-supply voltages.

## Installation

The original model lives at https://github.com/federico-pesamosca/plasmaless-tokamak-circuits. The IMAS interface, the MUSCLE3 actor and the DD4 machine description (`data/md_dd4`) are on branch `IMAS_Muscle3`; clone that version on ITER SDCC:

```bash
git clone -b IMAS_Muscle3 https://github.com/federico-pesamosca/plasmaless-tokamak-circuits.git
```

## Description

The model evolves the coupled circuit dynamics of the ITER active **circuits**, \(c\), and passive vessel conductive elements, \(e\). The 14 active coils are connected in 12 circuits: **CS1U and CS1L are in series** (one current Ic_CS1) and **VS3U and VS3L are in antiseries** (one current Ic_VS3, with Ia_VS3U = -Ic_VS3 and Ia_VS3L = +Ic_VS3, so a positive VS3 voltage gives a positive VS3L current). The other coils are one circuit each. With the coil currents Ia = T Ic (T is the 14x12 matrix of `plasmaless_circuit_topology`), the circuit matrices are

Mcc = T' Maa T, Mce = T' Mae, Rcc = T' Raa T (Mee, Ree unchanged)

and, for applied circuit voltages **Vc**, the circuit equations are

|&nbsp; Mcc&nbsp; Mce&nbsp; |&nbsp; |&nbsp; dIc/dt&nbsp; | + |&nbsp; Rcc&nbsp; 0&nbsp; |&nbsp; |&nbsp;Ic&nbsp;| = |&nbsp;Vc| \
|&nbsp; Mec&nbsp; Mee&nbsp; |&nbsp; |&nbsp; dIe/dt&nbsp; |&nbsp;&nbsp;&nbsp;   |&nbsp; 0&nbsp; Ree&nbsp; |&nbsp; |&nbsp;Ie&nbsp;|&nbsp;&nbsp;&nbsp; |&nbsp; 0&nbsp;|

Here, **M** is the mutual-inductance matrix, **R** contains the circuit resistances, and **Ic** and **Ie** are the active-circuit and passive-element currents (the state is `x = [Ic; Ie]`, 12 + 102 values). The input of the state space model are the **14 power-supply voltages**

Vps = {CS3U, CS2U, CS1U, CS1L, CS2L, CS3L, PF1, PF2, PF3, PF4, PF5, VS1, PF6, VS3}

mapped to the circuit voltages Vc = {CS3U, CS2U, CS1, CS2L, CS3L, PF1, PF2, PF3, PF4, PF5, PF6, VS3} by `Vc = map_powersupplies_to_circuits(Vps)` (matrix B of the state space model includes this map):

- Vc_CS1 = Vps_CS1U + Vps_CS1L (the two supplies are in series);
- the VS1 supply is in parallel to the supplies of PF2, PF3, PF4, PF5: Vc_PF2 = Vps_PF2 - Vps_VS1, Vc_PF3 = Vps_PF3 - Vps_VS1, Vc_PF4 = Vps_PF4 + Vps_VS1, Vc_PF5 = Vps_PF5 + Vps_VS1;
- Vc_VS3 = Vps_VS3; all other circuits have Vc = Vps.

The output of the state space model is unchanged: the **14 coil currents** Ia = {CS3U, CS2U, CS1U, CS1L, CS2L, CS3L, PF1..PF6, VS3U, VS3L} (Ia_CS1U = Ia_CS1L = Ic_CS1, Ia_VS3U = -Ic_VS3, Ia_VS3L = Ic_VS3) and the passive currents Ie. `map_coils_to_circuits(Ia)` gives the circuit currents from coil currents (for example for a warm start; it raises an error if the coil currents violate the series/antiseries constraints). Magnetic flux on the DINA grid is computed from the output currents as

psi = Mxi |&nbsp;I_a&nbsp;|
          |&nbsp;I_e&nbsp;|

where **Mxi** is the Green’s-function matrix mapping (single-coil) conductor currents to flux on that grid. For ITER, **Ra = 0** (superconductivity) is assumed for all coils **except the VS** (vertical stabilization) circuit, so Rcc is zero except Rcc_VS3 = R_VS3U + R_VS3L.

The Simulink model `plasmaless_timestep.slx` is included as a referenced model in `simulator_plasmaless_model_simulink.slx`. It evolves the circuit response to an input voltage time series Vps (it has no knowledge of the circuits: the matrices `A,B,C,D` carry the mapping, so the Simulink models did not need changes).

## DD4 machine description data

**The model uses IMAS Data Dictionary 4 (DD4) only.** The machine description (`em_coupling`, `pf_active`, `pf_passive`, `wall`) is read from the `.h5` files in `data/md_dd4`. The original DD3 IDSs (DINA `em_coupling`, `pf_active` and `equilibrium`, ITER_MD `pf_passive` and `wall`) are converted to DD 4.1.1 once, with `tools/convert_md_dd3_to_dd4.py`; no DD3 data is read when the model runs.

`data/md_dd4` is shipped with the repository. Regenerate it (for example for another DINA entry) with the conversion script, which needs the DD3 sources and imas-python:

```bash
module use /work/projects/pds/modules/all
module load IMAS-MUSCLE3/1.0.0-intel-2025b-pds
python tools/convert_md_dd3_to_dd4.py     # writes to data/md_dd4
# options: --dina-uri, --pf-passive-uri, --wall-uri, --output-uri "imas:hdf5?path=<dir>/md_dd4"
```

Default sources (DD3, read only):
- `em_coupling`, `pf_active`, `equilibrium` (the DINA grid): `imas:mdsplus?path=/work/projects/dina/SRO_JINTRAC/15MA_10perc/output`
- `pf_passive`: `imas:mdsplus?path=/work/imas/shared/imasdb/ITER_MD/3/115005/3`
- `wall`: `imas:mdsplus?path=/work/imas/shared/imasdb/ITER_MD/3/116000/2`

The script rereads what it wrote and checks it against the DD3 sources. Expected files in `data/md_dd4` (DD 4.1.1, HDF5):

| File | Content |
|---|---|
| `em_coupling.h5` | five `coupling_matrix` entries: `mutual_active_active` (14x14), `mutual_passive_active` (102x14), `mutual_passive_passive` (102x102), `mutual_grid_active` (8385x14), `mutual_grid_passive` (8385x102); and the (R,Z) mesh of the flux maps in `grid_ggd(1)` (8385 nodes = 65 R x 129 Z, Z varying fastest) |
| `pf_active.h5` | coil description (time-dependent data removed), coil resistances |
| `pf_passive.h5` | passive loops (102), names and resistances, geometry |
| `wall.h5` | limiter outlines used for plotting: `First Wall` (unit 1) and `Divertor` (unit 2) |
| `master.h5` | HDF5 master file linking the four IDSs above |

## MATLAB setup for simulation (SDCC)

`run_plasmaless_model` (Simulink and MATLAB modes), `run_plasmaless_imas`, the plotting and the tests run with MATLAB R2024a on SDCC:

```bash
module purge
module load IMAS-MATLAB/5.5.0-intel-2023b-DD-4.1.0
module unload MATLAB/2023b-r5-GCCcore-13.2.0
module load MATLAB/2024a-r6-GCCcore-13.2.0
matlab
```

Then, in MATLAB, from the repository root:

```matlab
add_plasmaless_paths     % adds the local subfolders (works whatever the name of the cloned folder)
run_plasmaless_model     % MATLAB or Simulink simulation
run_plasmaless_imas      % IMAS interface
```

## MUSCLE3 standalone run

```bash
muscle3/run_plasmaless_muscle3.sh     # loads its own modules
```

MATLAB R2025b hangs on figure rendering/export on the SDCC login nodes, so plotting is off by default (`do_plot = 0`); use the MATLAB R2024a setup above.

No other installation is required.

## MATLAB/Simulink basic simulation
- **Add the local paths**
run `add_plasmaless_paths.m` to add the local subfolders.
- **Run the simulation:** `run_plasmaless_model.m` loads the DD4 machine description (`load_plasmaless_md`), builds the model (`configure_plasmaless_model`, which returns `model.A/B/C/D` and the initial state `x0`) and runs `simulator_plasmaless_model_simulink.slx` (`A`, `B`, `C`, `D` are taken from the base workspace) or `simulator_plasmaless_model_matlab.m`. It saves `time`, `Ia`, `Ie` and the coil/loop names to `plasmaless_matlab_results.mat` in the repository root (same variables, units and names as `run_plasmaless_imas.m`), and plots the results with `plot_plasmaless_model.m` if `do_plot = 1` is set before the call (default 0). Adapt it to set the simulation parameters and input for your application.
  Set `mode` inside `run_plasmaless_model.m` to select the simulator:
  ```matlab
  mode = 'matlab';    % MATLAB simulation loop
  % or
  mode = 'simulink';  % Simulink simulation
  ```
  Only required parameters
  - dt: discretization time
  - tstart: start of simulation
  - tend: end of simulation
    
- **Prepare inputs:** `timeseries_plasmaless_model.m` illustrates the structure and format used to create the input voltage time series: a +/-10 V pulse on each of the 14 power supplies Vps.
- **Visualize results:** `plot_plasmaless_model.m` plots power-supply voltage and (all 14 coils and passive) current histories and animates the dynamic flux on the DINA grid. It reads the coil/loop geometry, the flux maps `Mxi`, the (R,Z) grid and the wall from `data/md_dd4` (`plasmaless_grid_from_ids`), and expects `out.simout` (timeseries of `[Ia Ie]`), `TS`, `tstart` and `tend` in the workspace; both run scripts provide them.
- **Run regression tests:** `tests/test_plasmaless_model.m` contains five tests:
  - a basic run of the scripts
  - the expected response to voltage steps in two superconductive circuits (CS3U, PF1; matrices Mcc, Mce, Mee)
  - the expected response to a voltage step on the resistive VS3 circuit (VS3U/VS3L in antiseries, settling to V/Rcc, positive VS3L current)
  - the circuit mappings (series, antiseries, `map_powersupplies_to_circuits`, `map_coils_to_circuits`)
  - MATLAB–Simulink output equivalence for the same input waveform and timing, checking matching timestamps and all active/passive current samples within numerical tolerance.

  `tests/test_vertical_instability.m` covers the vertical model (see `vertical_instability/README.md`).
  Two more scripts cover the IMAS interface: `tests/test_plasmaless_imas.m` (IMAS step loop vs Simulink, input/output contract errors, VS3 voltage step; run with `runtests` from the repository root) and `tests/check_dd3_dd4_equivalence.m` (the DD3 sources vs the DD4 entry `data/md_dd4`, two parts with two module stacks; see `howto.txt`).

Both simulation modes require MATLAB, Control System Toolbox (for `ss` and `c2d`), IMAS-MATLAB (DD 4), and the DD4 files in `data/md_dd4`. Simulink is additionally required for Simulink mode and for running the full regression test suite.

## IMAS interface (imas_model/)

The same model as plain MATLAB functions with IDS inputs and outputs, for use as a step in a coupled workflow. **Note:** the circuit/power-supply formulation above (`model.A/B/C/D` with 12 + 102 states and 14 power-supply inputs) has not yet been propagated to the IDS functions (`plasmaless_step_ids`, `simulator_plasmaless_model_imas`, `run_plasmaless_imas.m`, the MUSCLE3 actor), which still assume one state and one voltage input per coil; they need an update before use. The functions in `imas_model/` are:

| Function | Role |
|---|---|
| `load_plasmaless_md` | reads `em_coupling`, `pf_active`, `pf_passive` (and optionally `wall`) from `data/md_dd4` |
| `configure_plasmaless_model` | builds the discrete model from the IDSs (DD4 only); state = circuit currents, input = power-supply voltages |
| `plasmaless_circuit_topology` | coil/circuit/power-supply topology: matrices `T` (Ia = T Ic) and `P` (Vc = P Vps) |
| `map_powersupplies_to_circuits` | `Vc = map_powersupplies_to_circuits(Vps)`, 14 supplies to 12 circuits |
| `map_coils_to_circuits` | `Ic = map_coils_to_circuits(Ia)`, 14 coil currents to 12 circuit currents |
| `get_coupling_matrix` | reads and validates one `coupling_matrix` entry of `em_coupling` |
| `plasmaless_grid_from_ids` | returns the coil-to-grid flux map `Mxi` and the (R,Z) grid `RR`, `ZZ` |
| `plasmaless_step_ids` | one time step with IDS input and output |
| `simulator_plasmaless_model_imas` | runs the step loop for a voltage time series |

The model is built once, then stepped once per time step:

- Build once: `[em_coupling, pf_active, pf_passive] = load_plasmaless_md();` reads the DD4 machine description in `data/md_dd4`, then `[model, x0] = configure_plasmaless_model(em_coupling, pf_active, pf_passive, dt)` builds **M**, **R**, the continuous state space and its `c2d` discretization (`model.A/B/C/D`). Each coil and loop resistance is read from the IDS (`pf_active.coil{i}.resistance`, `pf_passive.loop{i}.resistance`). Missing resistances or inconsistent sizes raise an error. With these four inputs `x0` is zero (cold start). For a warm start, add a one-slice pair: `[model, x0] = configure_plasmaless_model(..., dt, pf_active_init, pf_passive_init)` takes the coil currents (matched by short name) and the loop currents (by position, names checked).
- Per step: `[pf_active_out, pf_passive_out, x] = plasmaless_step_ids(model, pf_active_in, x)` is the main function. It stands in for NICE in a coupled workflow and performs one time step, `x_new = A*x + B*u`, as `plasmaless_timestep.slx`.
  - Input: a one-slice `pf_active` (`homogeneous_time = 1`, a single time `t`) with one `coil{i}.voltage.data` value per coil. Coils are matched by short name (CS3U ... PF6, VS3U, VS3L). A missing coil, an extra coil or a missing voltage raises an error.
  - Output: `pf_active` with the coil currents and `pf_passive` with the loop currents, both at `t + dt`. The `pf_active` output also echoes the voltages applied over [t, t+dt].
- Flux on the DINA grid: `[Mxi, RR, ZZ] = plasmaless_grid_from_ids(em_coupling)` gives the map from `[Ia; Ie]` to the poloidal flux on the grid; `psi = reshape(Mxi*[Ia; Ie], size(RR))`.
- Example driver: `run_plasmaless_imas.m` (repo root, next to `run_plasmaless_model.m`; type `run_plasmaless_imas` from the repo root) uses the same input voltages as `run_plasmaless_model.m`, calls `simulator_plasmaless_model_imas` (one input slice per step), and can write the output slices with `ids_put_slice` (`write_output = 1`).
  It saves `time`, `Ia`, `Ie` and the coil/loop names to `results_file` (default `plasmaless_imas_results.mat` in the repository root; generated result files are ignored by git). Set `do_plot = 1` before the call to plot at the end with `plot_plasmaless_model.m`; `run_plasmaless_imas.m` builds the `out.simout` timeseries `[Ia Ie]` that this script expects.

## DD4 machine description

The `em_coupling` is the DD4 IDS with `coupling_matrix` entries named `mutual_active_active`, `mutual_passive_active`, `mutual_passive_passive`, `mutual_grid_active` and `mutual_grid_passive`; `get_coupling_matrix` checks the name, the quantity (magnetic flux) and that `rows_uri`/`columns_uri` refer to the expected nodes. The files in `data/md_dd4` (not pulse-specific; default of `run_plasmaless_imas.m`, `run_plasmaless_model.m`, the plots and the MUSCLE3 actor) are made from the DD3 IDSs by `tools/convert_md_dd3_to_dd4.py` (see Installation). DD3 IDSs are not read by the model: the DD4 IMAS-MATLAB API returns an empty `coupling_matrix` for a DD3 entry, so DD3 data must be converted first.

Limitations:
- The inputs are coil voltages only. There is no supply/circuit-to-coil mapping, so CS1U/CS1L and VS3U/VS3L are treated as independent coil voltages, as in the original model.
- The IMAS step functions return currents only; the flux on the DINA grid is computed in `plot_plasmaless_model.m` from `Mxi`.

## PDS integration

The model runs in IMAS-PDS as the `plasmaless` actor of the `plasmaless_controller` workflow (branch `feature/plasmaless_controller`), closed on the PCS magnetic controller. How to run it and the standalone checks: `instructions.txt`. How the integration was done: `pds_integration.txt`.

## Licence

This software is open source. ITER-related intellectual property remains &copy; 2026 ITER Organization.

## Contacts

Originally developed in 2026 by Federico Pesamosca (federico.pesamosca@iter.org)

IMAS interfacing and MUSCLE3 adaptation by Mireille Schneider (mireille.schneider@iter.org)
