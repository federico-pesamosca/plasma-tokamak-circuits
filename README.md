# Plasmaless Tokamak Circuits

A dynamical circuit model of the current response of ITER coils and vessel conductive elements to applied active-coil voltages.

## Installation

The original model lives at https://github.com/federico-pesamosca/plasma-tokamak-circuits. The IMAS interface, the MUSCLE3 actor, the DD4 machine description (`data/md_dd4`), the external linear model and the vertical-instability model are on branch `IMAS_Muscle3_with_plasma`; it also contains the plasma-free model. Branch `IMAS_Muscle3` holds the plasma-free version only. Clone the extended version on ITER SDCC:

```bash
git clone -b IMAS_Muscle3_with_plasma https://github.com/federico-pesamosca/plasma-tokamak-circuits.git
```

## Description

The model evolves the coupled circuit dynamics of the ITER active coils, \(a\), and passive vessel conductive elements, \(e\). For applied active-coil voltages **Va**, the circuit equations are

|&nbsp; Maa&nbsp; Mae&nbsp; |&nbsp; |&nbsp; dIa/dt&nbsp; | + |&nbsp; Ra&nbsp; 0&nbsp; |&nbsp; |&nbsp;Ia&nbsp;| = |&nbsp;Va| \
|&nbsp; Mea&nbsp; Mee&nbsp; |&nbsp; |&nbsp; dIe/dt&nbsp; |&nbsp;&nbsp;&nbsp;   |&nbsp; 0&nbsp; Re&nbsp; |&nbsp; |&nbsp;Ie&nbsp;|&nbsp;&nbsp;&nbsp; |&nbsp; 0&nbsp;|

Here, **M** is the mutual-inductance matrix, **R** contains the circuit resistances, and **Ia** and **Ie** are the active- and passive-element currents. Magnetic flux on the DINA grid is computed from the currents as

psi = Mxi |&nbsp;I_a&nbsp;|
          |&nbsp;I_e&nbsp;|

where **Mxi** is the Green’s-function matrix mapping conductor currents to flux on that grid. For ITER, **Ra = 0** (superconductivity) is assumed for all coils **except the VS** (vertical stabilization) circuit.

The Simulink model `plasmaless_timestep.slx` is included as a referenced model in `simulator_plasmaless_model_simulink.slx`. It evolves the circuit response to an input voltage time series Va.

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
    
- **Prepare inputs:** `timeseries_plasmaless_model.m` illustrates the structure and format used to create the input voltage time series.
- **Visualize results:** `plot_plasmaless_model.m` plots coil-voltage and current histories and animates the dynamic flux on the DINA grid. It reads the coil/loop geometry, the flux maps `Mxi`, the (R,Z) grid and the wall from `data/md_dd4` (`plasmaless_grid_from_ids`), and expects `out.simout` (timeseries of `[Ia Ie]`), `TS`, `tstart` and `tend` in the workspace; both run scripts provide them.
- **Run regression tests:** `tests/test_plasmaless_model.m` contains four tests:
  - a basic run of the scripts
  - the expected response to voltage steps in two superconductive coils
  - the expected response to voltage steps in two resistive coils
  - MATLAB–Simulink output equivalence for the same input waveform and timing, checking matching timestamps and all active/passive current samples within numerical tolerance.

  Two more scripts cover the IMAS interface: `tests/test_plasmaless_imas.m` (IMAS step loop vs Simulink, input/output contract errors, VS3 voltage step; run with `runtests` from the repository root) and `tests/check_dd3_dd4_equivalence.m` (the DD3 sources vs the DD4 entry `data/md_dd4`, two parts with two module stacks; see `howto.txt`).

Both simulation modes require MATLAB, Control System Toolbox (for `ss` and `c2d`), IMAS-MATLAB (DD 4), and the DD4 files in `data/md_dd4`. Simulink is additionally required for Simulink mode and for running the full regression test suite.

## IMAS interface (imas_model/)

The same model as plain MATLAB functions with IDS inputs and outputs, for use as a step in a coupled workflow. The functions in `imas_model/` are:

| Function | Role |
|---|---|
| `load_plasmaless_md` | reads `em_coupling`, `pf_active`, `pf_passive` (and optionally `wall`) from `data/md_dd4` |
| `configure_plasmaless_model` | builds the discrete model from the IDSs (DD4 only) |
| `get_coupling_matrix` | reads and validates one `coupling_matrix` entry of `em_coupling` |
| `plasmaless_grid_from_ids` | returns the coil-to-grid flux map `Mxi` and the (R,Z) grid `RR`, `ZZ` |
| `plasmaless_step_ids` | one time step with IDS input and output (circuit or linear model) |
| `simulator_plasmaless_model_imas` | runs the step loop for a voltage time series (repository root) |
| `plasmaless_coil_key` | short coil name of a `pf_active.coil` element (shared by all functions and the actor) |
| `configure_linear_model` | builds the discrete model from an external linear model (`linear_model.json` + matrices) |
| `read_npy` | minimal reader of numpy `.npy` files (`<f8`/`<i8`, 1-2 dimensions) |
| `extract_equilibrium_data` (`vertical_instability/`) | reads an equilibrium (any IMAS URI or entry directory, DD3 or DD4) and its conductor currents for the vertical-instability model |
| `configure_vertical_instability_model` (`vertical_instability/`) | builds the reduced vertical RZIp model (kind `'vertical'`) around that equilibrium |
| `run_vertical_instability`, `plot_vertical_instability` (repository root) | open-loop voltage-step demo and flux/LCFS animation of the vertical model |
| `plasmaless_observables_to_equilibrium` | writes the known observables of a linear model into a one-slice `equilibrium` IDS |
| `plasmaless_observables_from_equilibrium` | reads them back from an `equilibrium` time slice (operating point of the actor) |

The model is built once, then stepped once per time step:

- Build once: `[em_coupling, pf_active, pf_passive] = load_plasmaless_md();` reads the DD4 machine description in `data/md_dd4`, then `[model, x0] = configure_plasmaless_model(em_coupling, pf_active, pf_passive, dt)` builds **M**, **R**, the continuous state space and its `c2d` discretization (`model.A/B/C/D`). Each coil and loop resistance is read from the IDS (`pf_active.coil{i}.resistance`, `pf_passive.loop{i}.resistance`). Missing resistances or inconsistent sizes raise an error. With these four inputs `x0` is zero (cold start). For a warm start, add a one-slice pair: `[model, x0] = configure_plasmaless_model(..., dt, pf_active_init, pf_passive_init)` takes the coil currents (matched by short name) and the loop currents (by position, names checked).
- Per step: `[pf_active_out, pf_passive_out, x, y_out] = plasmaless_step_ids(model, pf_active_in, x)` is the main function. It stands in for NICE in a coupled workflow and performs one time step, `x_new = A*x + B*(u - op.voltage)`, as `plasmaless_timestep.slx` (`op` is the operating point, zero for the circuit model).
  - Input: a one-slice `pf_active` (`homogeneous_time = 1`, a single time `t`) with one `coil{i}.voltage.data` value per coil. Coils are matched by short name (CS3U ... PF6, VS3U, VS3L). A missing coil, an extra coil or a missing voltage raises an error.
  - Output: `pf_active` with the coil currents (`op.coil_current + x_to_coil*x_new`) and `pf_passive` with the loop currents (`op.loop_current + x_to_loop*x_new`), both at `t + dt`. The `pf_active` output also echoes the voltages applied over [t, t+dt]. `y_out` (optional 4th output) holds the observables of a linear or vertical model (`names`, `values = op.y + C_y*x_new`, `known`); it is empty for the circuit model.
- Model structure: the three builders (`configure_plasmaless_model`, `configure_linear_model`, `configure_vertical_instability_model`) return the same fields: `kind` (`'circuit'`, `'linear'` or `'vertical'`), `sysd`, `A`, `B`, `C`, `D`, `dt`, `n_states`, `state_names`, `coil_names` (machine-description order), `n_coils`, `x_to_coil` (coil currents = `x_to_coil*x`), `loop_names`, `n_loops`, `x_to_loop`, `y_names`, `C_y` (observables = `op.y + C_y*x`, 0 x n for the circuit model) and `op` (operating point: `coil_current`, `loop_current` (zero for the circuit and linear models, the equilibrium passive currents for the vertical model), `voltage`, `y`, `y_known`, plus the `*_given` flags of the descriptor entries). The circuit model keeps `M` and `R`; its selections are `x_to_coil = [eye(14) zeros(14,102)]`, `x_to_loop = [zeros(102,14) eye(102)]` and the step results are unchanged.
- Flux on the DINA grid: `[Mxi, RR, ZZ] = plasmaless_grid_from_ids(em_coupling)` gives the map from `[Ia; Ie]` to the poloidal flux on the grid; `psi = reshape(Mxi*[Ia; Ie], size(RR))`.
- Example driver: `run_plasmaless_imas.m` (repo root, next to `run_plasmaless_model.m`; type `run_plasmaless_imas` from the repo root) uses the same input voltages as `run_plasmaless_model.m`, calls `simulator_plasmaless_model_imas` (one input slice per step), and can write the output slices with `ids_put_slice` (`write_output = 1`).
  It saves `time`, `Ia`, `Ie`, `Y`, the coil/loop/observable names and `model_kind` to `results_file` (default `plasmaless_imas_results.mat` in the repository root; generated result files are ignored by git). Set `do_plot = 1` before the call to plot at the end with `plot_plasmaless_model.m`; `run_plasmaless_imas.m` builds the `out.simout` timeseries `[Ia Ie]` that this script expects.

## External linear model (NICE linearisation)

An externally supplied continuous linear model `dI/dt = A I + B U`, `Y = C I` (any state dimension) can replace the circuit model everywhere: IMAS step interface, run scripts, Simulink/MATLAB simulators and the MUSCLE3 actor. The received model `linearModelSimple_105073_t84.417598/` (NICE linearisation delivered under the label 105078, confirmed by the user on 2026-10-14 as pulse 105073 at t = 84.417598 s: 17 states, 14 generator voltages, 9 outputs) is the first instance; see its `README.md`.

**Descriptor** `linear_model.json` in the model directory (read by `configure_linear_model`, `jsondecode`):

| Field | Content |
|---|---|
| `description`, `shot`, `time` | free text and the linearisation point |
| `files` | `A`, `B`, `C`: basenames of the matrices, `.npy` (read by `read_npy`) or whitespace-separated `.txt` (`readmatrix`) |
| `states` | n state names (rows of A) |
| `inputs` | m input names, in the column order of B (the received model: the 14 generator voltages `CS3U, CS2U, CS1U, CS1L, CS2L, CS3L, PF1 .. PF5, VS1, PF6, VS3`) |
| `input_map` | object input name -> {coil short name: coefficient}: the input as a linear combination of the 14 `pf_active` coil voltages (coil sense). An input absent from `input_map`, or mapped to an empty object `{}`, gets a zero row: it cannot be driven through `pf_active` (warning `plasmaless:inputUndriven`, once at build time). Without `input_map`, `inputs` must be exactly the coil set and the map is the permutation to machine-description order (the previous format) |
| `coil_states` | object coil short name -> state name (several coils may share one state, e.g. CS1U and CS1L -> CS1); keys = machine-description coils |
| `passive_states` | state names of the passive conductors; every state is a coil state or a passive state, not both |
| `outputs` | p output names (rows of C), in the vocabulary below |
| `operating_point` (optional) | `coil_current` (coil -> A), `voltage` (coil -> V), `outputs` (output -> value): absolute values at the linearisation point, keyed by coil (coil space), in the wrapper's (IMAS) convention; see the MUSCLE3 setting `operating_point` for how they combine with the F_INIT data |
| `convention` (optional) | sign frame of the matrix files relative to the wrapper's (IMAS) convention: `state_sign` (+1/-1), `input_sign` (+1/-1), `output_sign` (object output name -> +1/-1); everything absent is +1. Meaning: `x_file = state_sign * x`, `u_file = input_sign * u`, `y_file(k) = output_sign(k) * y(k)` with x, u, y in the wrapper's convention. Applied before the input map and the discretisation: A unchanged, `B_conv = (input_sign/state_sign) * B`, `C_conv(k,:) = (state_sign/output_sign(k)) * C(k,:)`. Values other than +1/-1, unknown output names and unknown fields are `plasmaless:descriptor` errors. The raw file matrices are kept in `A_file`, `B_file`, `C_file`; `A_c`, `B_c`, `B_c_input`, `C_c` and the discrete matrices are in the wrapper's convention. The received descriptor: `state_sign -1`, `input_sign -1`, `output_sign` Ip -1 and +1 for the geometric outputs: A and B unchanged, the Ip row of C unchanged, the geometric rows of C negated |
| `notes`, `units`, `provenance`, `provider_note` | free text |

Checks (`plasmaless:descriptor`): input names unique, `numel(inputs)` = columns of B (`plasmaless:size`), every `input_map` key is an input, every coil referenced exists in the machine description, coefficients finite scalars. Other errors: `plasmaless:size` (matrix sizes, non-zero discrete D), `plasmaless:observable` (output name), `plasmaless:file` (missing matrix); `plasmaless:unstable` is a warning when A has an eigenvalue with real part > 1e-9.

**Generator inputs of the received model.** `G` (m x 14) is built from `input_map`; the continuous input matrix in coil space is `B_coil = B*G` and the discrete system is `c2d` of `(A, B_coil, C)`, so the step function and the simulators always take the 14 coil voltages. The model stores `input_names`, `input_map` (= G), `B_c` (continuous, coil space) and `B_c_input` (continuous, original input order). In the received descriptor each CS/PF generator maps to its coil with coefficient 1 (PF6 is input column 13); `VS3` (column 14, one voltage for the anti-series VS3U/VS3L pair) maps to `{VS3U: -1, VS3L: 1}` (generator voltage = V(VS3L) - V(VS3U) in coil sense; derived on 2026-10-14 from the provider's positive-voltage-positive-current rule applied to the delivered B, whose VS3 column drives the `-VS3` state (VS3L) positively and the `VS3` state (VS3U) negatively; it supersedes the earlier `{VS3U: 1, VS3L: -1}`, and an explicit confirmation of the loop sense is still pending); `VS1` (column 12, the star-point converter in series with PF2 .. PF5) maps to `{}`: no `pf_active` coil carries it, the column is not drivable (the PDS controller bridge drops the VS1 command anyway).

**Output vocabulary** and the DD4 (4.1.x) `equilibrium` fields they map to (`plasmaless_observables_to_equilibrium`, `time_slice{1}`):

| Name | Field |
|---|---|
| `rmag`, `zmag` | `global_quantities.magnetic_axis.r` / `.z` |
| `rX`, `zX` | `contour_tree.node{1}.r` / `.z` with `critical_type` 1 (saddle = X-point); DD4 has no `boundary_separatrix/x_point` any more |
| `Ip` | `global_quantities.ip` |
| `rc`, `zc` | `global_quantities.current_centre.r` / `.z` (current centre; vertical model) |
| `gap<k>` | `boundary.gap{k}.value`, `name` = `gap<k>` |
| `Dist2ndSep` | not written: DD4 removed `boundary_secondary_separatrix/distance_inner_outer` without replacement (kept in `y_out` and in the results files) |

**Semantics.** The model's states, inputs and outputs are deviations from the linearisation point: `du = u - op.voltage`, `x_new = A x + B du`, coil currents `op.coil_current + x_to_coil*x_new`, passive currents `x_to_loop*x_new` (one `pf_passive` loop per passive state, 4 for the received model: `PS_1, PS_2, VV_1, VV_2`), observables `op.y + C_y*x_new` (sent only where `op.y_known`). `x_to_coil` has identical rows for coils sharing a state (CS1U/CS1L); VS3U reads the state `VS3` and VS3L the state `-VS3`. The run scripts take the operating point from the descriptor (`op.voltage` = 0 for the received model, so the pulse train is applied as is) and save **absolute** values: `Ia` = `op.coil_current` + deviation, `Y` = `op.y` + deviation (`op.y` = 0 for the outputs the descriptor does not give). Without `operating_point` in the descriptor the operating point is zero and the results are deviations.

**How to run.**
- `run_plasmaless_imas.m` / `run_plasmaless_model.m`: set `linear_model_dir` in the script (default: the environment variable `PLASMALESS_LINEAR_MODEL_DIR` when set, else `''` = circuit model), e.g. `export PLASMALESS_LINEAR_MODEL_DIR=$PWD/linearModelSimple_105073_t84.417598`. Results: `plasmaless_linear_results.mat` / `plasmaless_matlab_linear_results.mat` (`time`, `Ia`, `Ie`, `Y`, `coil_names`, `loop_names`, `y_names`, `model_kind`). The Simulink/MATLAB simulators get `C_full = [x_to_coil; x_to_loop; C_y]` and feed `u` directly (no `op.voltage` subtraction, no `op` offsets: their results are deviations). With `do_plot = 1` the plots come from `plot_plasmaless_linear.m` (coil-current deviations, passive currents, one subplot per observable).
- MUSCLE3: setting `plasmaless.linear_model_dir` (absolute path, default `''`); `muscle3/run_plasmaless_muscle3.sh <model dir>` adds it to the generated `paths.ymmsl` and saves `plasmaless_muscle3_linear_results.mat`. The actor starts at t0 with x = 0; its operating point combines the descriptor and the F_INIT data at t0 (coil currents/voltages of the pf_active slice closest to t0, observables found in the equilibrium slice closest to t0) according to the optional setting `plasmaless.operating_point`: `descriptor` (default; descriptor values win, F_INIT fills what the descriptor does not give) or `f_init` (F_INIT values win, descriptor values only for what F_INIT does not have; any other value is the error `plasmaless_actor:settings`). Consequence: with `descriptor` the absolute coil currents and Ip sent to the controller start at the linearisation point, which differs from the F_INIT scenario values, so the controller sees a step error at t0; with `f_init` they start on the scenario. Each step uses `du = V - op.voltage`, so the operating point is held until the controller changes the voltages; `equilibrium_o_i` carries the known observables (a known Ip overrides the reference pass-through ip; `boundary.geometric_axis` stays the reference pass-through). See `muscle3/README.md`.
- Tests: `tests/test_plasmaless_linear.m`, eight tests (build, input map and state mapping, step loop vs direct recursion, operating point, descriptor errors, descriptor without `input_map`, equilibrium round trip, circuit model unchanged, one bipolar pulse per generator on the continuous file matrices: own-state response of each CS/PF generator, CS1U = CS1L, VS1 acting on PF2-PF5, VS3 acting on the VS3 pair with opposite signs).

**Status of the received model** (updated 2026-10-14, third update: see the end of this paragraph). Pulse correction (confirmed by the user on 2026-10-14): the provider's e-mail and the original directory name said 105078, but the data match pulse **105073** at t = 84.4176 s (DINA 105073 equilibrium slice at 84.4176 s: ip = -15.0632 MA, the provider's 1.50632e7 A to the last digit with the opposite sign; magnetic axis (6.2996, 0.6115) m vs the provider's (6.29704, 0.596501) m, the 79.2176 s slice has (6.2977, 0.5958) m; neighbouring slices 79.2176 s -15.0732 MA, 89.7176 s -15.0547 MA; entry `imas:hdf5?path=/work/projects/pds/pds-scenarios/105073/data/in`, DD 4.1.1), whereas DINA 105078 has -7.48 MA there. The descriptor's `shot` is 105073 (`provider_note` records the 105078 label); the directory was renamed from `linearModelSimple_105073_t84.417598` to `linearModelSimple_105073_t84.417598` on 2026-10-14. Sign frame (derived, to be confirmed by the provider): the DINA 105073 `pf_active` currents at 84.4376 s (closest sample: CS3U -22500, CS2U 9300, CS1U = CS1L 34771, CS2L 5428, CS3L -36059, PF1 -22192, PF2 29143, PF3 20117, PF4 40774, PF5 9766, PF6 -19143, VS3U -204, VS3L 204 A) are the provider's operating-point currents with the opposite sign and the same magnitude within about 3 % (CS2L 16 %): the matrices and the provider's operating point are written in a frame where all currents (coils, passive, Ip) have the opposite sign to the IMAS data (global negation of the current sense, as between a +Ip frame and the IMAS COCOS 17 data), geometric quantities unaffected. Hence the descriptor's `convention` block: `state_sign -1`, `input_sign -1` (ASSUMPTION: the generator voltages flip with the current sense, as power = V*I; to be confirmed), `output_sign` Ip -1 and +1 for rmag, zmag, rX, zX, gap1..3, Dist2ndSep; effect: A and B unchanged, the Ip row of C unchanged, the geometric rows of C negated, so the sign of the Ip output row is consistent with the IMAS data by construction; the VS3 generator-sign question is unchanged by the frame flip (B is unchanged when both x and u flip). Operating point in the descriptor (IMAS convention): the negated provider currents (CS3U -23248.2, CS2U 9514.24, CS1U = CS1L 34331.9, CS2L 6293.26, CS3L -37480.8, PF1 -21873.6, PF2 29058.6, PF3 19884.8, PF4 40412.3, PF5 10061.6, PF6 -19263.3, VS3U -203.564, VS3L 203.564 A), U = 0, rmag 6.29704 m, zmag 0.596501 m, Ip -15063200 A; its origin is thereby explained: the NICE equilibrium of 105073 at 84.4176 s, within about 3 % of the DINA 105073 coil currents (opposite frame). SUPERSEDED (2026-10-14, kept below for the record): the reading of the model as pulse 105078, the "factor two / typo" reading of the provider's 15.06 MA against the DINA 105078 Ip of -7.48 MA, and the user decision to set the operating-point Ip to -7478000 A: the pulse number was wrong, not the value. Still to confirm with the provider: the sign frame (incl. `input_sign`), the VS3 loop sense explicitly, the VS3 pair dynamics (pulse 105073 confirmed by the user; the provider's second answer and the VS3 map flip are in the third update at the end of this paragraph). Second run of the pulse train with the convention block (dt 0.001, 28 s): 35.9 s wall time; all responses finite, bounded and decaying (end values < 1 % of the peaks); max |coil current deviation| CS3U 7.6, CS2U 7.9, CS1U = CS1L 2.8, CS2L 7.9, CS3L 7.5, PF1 8.1, PF2 13.6, PF3 3.1, PF4 3.8, PF5 3.9, PF6 2.4, VS3U 773, VS3L 711 A; PS_1 745, PS_2 1045, VV_1 617, VV_2 198 A; max |observable deviation| rmag 8.0e-5, zmag 2.8e-3, rX 5.4e-4, zX 1.2e-3 m, Ip 9.5e3 A, gap1..3 2.3e-4, 9.5e-5, 7.3e-5 m, Dist2ndSep 6.0e-4 m (saved absolute values: coil currents = negated operating point + deviation, Ip = -15.0632 MA + deviation). Signed comparison with the previous run (file matrices used as is): the coil, passive and Ip deviations are identical (difference 0; VS3 1e-13 A), the geometric observable deviations (rmag, zmag, rX, zX, gaps, Dist2ndSep) have the opposite sign, exactly as the convention formulas prescribe (B and the Ip row of C unchanged, geometric rows negated); all magnitudes unchanged. Tests: `tests/test_plasmaless_linear.m` 7/7, `tests/test_plasmaless_imas.m` 3/3, `tests/test_plasmaless_model.m` 4/4 (MATLAB R2024a, IMAS-MATLAB 5.5.0 DD 4.1.0). Third update (2026-10-14, provider's second answer, relayed by the user, translated): (1) on the sign frame, "the generators and the currents are linked by the circuit equations; normally a positive voltage must give a positive current, unless COCOS changes signs, which would have to be tested" (voltages and currents share the sense in the provider's frame, consistent with `input_sign = state_sign`; the frame itself is not settled, a test is suggested); (2) "VS3 and -VS3 mean VS3U and VS3L, which are currents driven by a single generator and placed in a circuit such that their signs are opposite"; (3) on the non-mirror rows, "perhaps with the linearised model the currents end up decorrelated, and that is what bothers him" (not considered normal; no fix). Own check on the delivered B (provider frame, inputs in the confirmed generator order): the own-generator responses are CS3U 1.522, CS2U 1.640, CS1 0.578 (per CS1U or CS1L volt), CS2L 1.631, CS3L 1.516, PF1 1.606, PF2 2.814, PF3 0.647, PF4 0.814, PF5 0.811, PF6 0.492 A/s per V (all positive, consistent with the rule); state `VS3` (VS3U) responds with -3192.8 A/s per V of the VS3 generator, state `-VS3` (VS3L) with +1734.7, so the positive sense of the VS3 generator is the VS3L current direction and the generator voltage in coil-sense `pf_active` voltages is V(VS3L) - V(VS3U): `input_map.VS3` changed to `{VS3U: -1, VS3L: 1}` (derived; the previous `{VS3U: 1, VS3L: -1}` is superseded; an explicit confirmation of the loop sense is still pending). The two responses are not opposite-equal (ratio -1.84), the decorrelation the provider mentions; unresolved. `tests/test_plasmaless_linear.m` updated (B_c(:,VS3U) = -B(:,14), B_c(:,VS3L) = +B(:,14)) and re-run: 7/7 passed (MATLAB R2024a, IMAS-MATLAB 5.5.0 DD 4.1.0).

**PDS closed-loop test (2026-10-14).** The linear-model path was run in a PDS case with the PCS magnetic controller: IMAS-PDS worktree `IMAS-PDS-plasmaless` (branch `feature/plasmaless_controller`), case `plasmaless_controller_105073` (override `cases/overrides/plasmaless_controller_105073.ymmsl`: `source.t_min` 81.0 / `t_max` 100.0 on the `metis_nice_inverse_from_dina_105073` inverse output of the same day, t0 = 81.4576 s, `plasmaless.t_end` 83.4576, dt 0.005, `plasmaless.linear_model_dir` = this model directory, `plasmaless.operating_point` descriptor), Slurm run `cases/runs/plasmaless_controller_105073_20261014_103330`: finished without error, 400 exchanges, all outputs finite. It is an integration test (settings, coil order, ports, outputs), not a validation of the physics: the CS/PF commands sit at the clip from the first step as in every `plasmaless_controller` case and the VS3 currents reach -100 kA; numbers and caveats in `pds_integration.txt` 4.7. The PDS local install must be refreshed (or `PLASMALESS_REPO` exported at submission) before the linear model is available there.

Earlier status text (2026-10-14, first update; the 105078 and -7.478 MA statements are SUPERSEDED, see above). Verified by the model provider (2026-10-14): the 17 states are the coil currents `CS3U, CS2U, CS1, CS2L, CS3L, PF1 .. PF6, VS3, -VS3` (state 13 = the negated VS3 current, i.e. the VS3L current), then 2 passive-structure and 2 vacuum-vessel currents; the 14 inputs are the generator voltages `CS3U, CS2U, CS1U, CS1L, CS2L, CS3L, PF1 .. PF5, VS1, PF6, VS3` (VS1 = star-point converter in series with PF2 .. PF5, VS3 = one voltage for the anti-series pair; this explains the PF6 self term under column 13 and column 12 acting on PF2 .. PF5 and the passives); I, U, Y are deviations from the equilibrium at 84.4176 s; operating point: the coil currents of the descriptor (CS3U 23248.2, CS2U -9514.24, CS1 -34331.9, CS2L -6293.26, CS3L 37480.8, PF1 21873.6, PF2 -29058.6, PF3 -19884.8, PF4 -40412.3, PF5 -10061.6, PF6 19263.3, VS3 203.564 A), U = 0, rmag 6.29704 m, zmag 0.596501 m, Ip 1.50632e7 A as sent by the provider; rX, zX, gaps and Dist2ndSep are not used for now and their definitions are not needed. The output order is from `linModOrder.txt`. Operating-point Ip in the descriptor (user decision, 2026-10-14, not a provider confirmation): the provider's 15063200 A is treated as a typo and `operating_point.outputs.Ip` is set to -7478000 A, the DINA 105078 plasma current at 84.417598 s in the IMAS sign convention (linear interpolation of the DINA equilibrium ip, flat between 73.6 s (-7.474 MA) and 95.0 s (-7.498 MA); entry `imas:hdf5?path=/work/projects/pds/pds-scenarios/105078/data/in`, DD 4.1.1); the sign of the model's Ip output row (C row 5) is kept as delivered (assumption: the model's Ip deviations are in the IMAS sign convention); value and sign still to be confirmed with the provider, so the statements on the 15.06 MA value below are superseded in the descriptor but the question stays open. Still UNVERIFIED: the VS3 voltage map `{VS3U: 1, VS3L: -1}` (derived from the anti-series connection and the anti-symmetric VS3U/VS3L voltages of the DINA 105078 data, -5.96/+5.96 V); VS1 undriven; the sign convention of Ip (and possibly of the coil currents) relative to IMAS: the provider's Ip is +15.06 MA whereas the IMAS data of this pulse (`imas:hdf5?path=/work/projects/pds/pds-scenarios/105078/data/in`, DD 4.1.1, 40 samples) have Ip = -7.48 MA, so a sign check against the circuit model is not meaningful yet; the origin of the operating point: it is not the DINA 105078 state near 84.4 s (closest DINA sample 87.88 s: Ip -7.48 MA, magnetic axis (6.33, 0.54) m, CS3U -8.9 kA, CS1U = CS1L +11.4 kA, VS3U/VS3L +-2668 A), presumably a NICE equilibrium of another run; and the VS3 pair dynamics (A columns 12/13 are exact negatives but rows 12/13 of A and B are not, so the VS3U and VS3L deviations of a run are not mirror images). Run of the pulse train with the new map (dt 0.001, 28 s, 42 s wall time): all responses finite, bounded and decaying (end values < 1 % of the peaks); coil-current deviations up to 13.6 A (CS/PF; CS1U = CS1L 2.8 A) and 773/711 A (VS3U/VS3L), passive currents up to 1.05 kA, Ip deviation up to 9.5 kA, positions and gaps below 3 mm (run made with the provider's Ip 15.06 MA in the descriptor; the deviations do not depend on the operating-point value, only the absolute Ip saved by `run_plasmaless_imas.m` and sent by the actor does).

## Vertical-instability model

`vertical_instability/` holds a reduced RZIp model of the rigid vertical displacement of the plasma around a given equilibrium (theory, equations and references in `vertical_instability/README.md`). It extends the circuit model (same 14 coils and 102 passive loops, same order, voltages of the 14 coils as inputs) with the instantaneous vertical force balance (singular perturbation): the states are the current deviations `[dIa; dIe]` from the equilibrium currents `Ic0`, the outputs are `[zIp; Ia; Ie]` with `zIp = Ip0*delta_z` [A m]. The plasma current is held fixed and radial motion is excluded. For an ITER equilibrium the model has exactly one unstable continuous pole (the vertical instability).

**Loading the equilibrium** (`extract_equilibrium_data(uri, time, opts)`): any IMAS URI readable by the loaded IMAS-MATLAB module, or a plain path to an IMAS HDF5 entry directory (or to its `equilibrium.h5`), converted to `imas:hdf5?path=<dir>`. The equilibrium needs `profiles_2d` (r, z, psi and `j_phi`, or the Grad-Shafranov profiles `dpressure_dpsi` / `f_df_dpsi` of `profiles_1d`), `global_quantities.ip`, the LCFS and `pf_active` / `pf_passive` currents. Tested with:

| Entry | DD | Content | Growth rate | `zc0` vs entry `current_centre.z` |
|---|---|---|---|---|
| `imas:mdsplus?path=/work/projects/dina/imasdb/iter/3/105108/2`, 100 s (the model author's default) | 3.39.0 (read with the DD-4.1.0 module) | no `j_phi` (Grad-Shafranov reconstruction), 102 passive currents | 3.172 1/s | 0.521388 vs 0.524346 m (-3.0 mm) |
| `/work/projects/pds/pds-scenarios/105073/data/in`, 84.4176 s (PDS scenario entry) | 4.0.0 | `j_phi` (DINA), no passive currents, stored LCFS outline unusable | 7.617 1/s (`passive_currents` `'zero'`) | 0.590943 vs 0.590917 m (+0.03 mm) |

- psi sign: the machine-description flux maps need `psi0 = -psi` for a DD3 entry (COCOS 11, the author's convention) and `psi0 = +psi` for a DD4 entry (COCOS 17). The sign comes from the entry's DD version and is checked with the author's Ampere-law test against the signed Ip (`vertical:fluxSign`): DD3 105108 gives -7.375 MA for Ip -7.474 MA, DD4 105073 gives -15.053 MA for Ip -15.063 MA (with the opposite sign, +15.053 MA: rejected).
- LCFS: `boundary.outline` when it is a closed curve inside the grid; otherwise the closed contour of the slice's psi at `global_quantities.psi_boundary` around the magnetic axis (`equilibrium.boundary_source` says which). The 105073 entry needs this: its interpolated slice has no outline (linear interpolation drops it, neighbouring samples have different sizes) and the stored outlines run to the grid corner (open contour at the X-point flux).
- Passive currents: with `opts.passive_currents = 'equilibrium'` (default) the entry's `pf_passive` must hold one current per machine-description loop (same names and order), else error `vertical:passiveCurrents`; `'zero'` sets them to 0 explicitly, recorded in `equilibrium.passive_source` (`'zero (not in entry)'`). The PDS scenario entries have no passive currents (2 loops, no time). NICE inverse outputs (`nice_out`, PDS `initial_state`) have no `profiles_2d` (only `ggd`) and cannot feed this model.
- Other options: `interpolation` (`ids_get_slice` mode, default 3), `pf_passive_md` (machine-description `pf_passive`, default `load_plasmaless_md()`).

**IMAS interface.** The model has the common fields (kind `'vertical'`, `x_to_coil`, `x_to_loop`, ...) and is stepped by `plasmaless_step_ids` like the other two: coil currents `op.coil_current + x_to_coil*x` and loop currents `op.loop_current + x_to_loop*x` with `op.coil_current`/`op.loop_current` = `Ic0`; observables `y_names = {'Ip','rc','zc'}`, `C_y = [0; 0; C(ind_zIp,:)/Ip0]` (discrete C row of zIp divided by Ip0), `op.y = [Ip0; rc0; zc0]`, so Ip and rc stay constant and `zc = zc0 + zIp/Ip0`. `zc0`, `rc0` (`model.zc0`, `model.rc0`) are the centre of the equilibrium current distribution, `sum(J.*area.*Z)/Ip0` with the same area weights as the force balance. `op.voltage` = 0 and not given: the equilibrium has no voltage reference. The author's fields (`sys`, `Ac`..`Dc`, `M11`..`M22`, `Czr`, `Ip0`, `Ic0`, `ind_*`) are unchanged, so `run_vertical_instability.m` and `plot_vertical_instability.m` work as before.

**Why `current_centre`.** The model output is the displacement of the current distribution's centre (rigid shift of J), not of the magnetic axis or of the boundary: the DD4 field for it is `equilibrium.time_slice{1}.global_quantities.current_centre.z` ("vertical coordinate of current centre", m), with `.r` and `.velocity_z`; DINA fills it in the scenario entries (105073 at 84.4176 s: z 0.5909 m, r 6.1206 m), and the model's `zc0` reproduces it to 0.03 mm there.

**MUSCLE3 actor** (`muscle3/muscle_plasmaless_actor.m`): settings `vertical_equilibrium_uri` (default `''` = not used), `vertical_equilibrium_time` (default t0), `vertical_passive_currents` (`equilibrium` default, or `zero`); setting it together with `linear_model_dir` is an error. `operating_point: descriptor` means the model's own reference (here: the equilibrium's coil and passive currents, Ip, rc, zc; voltages, which the equilibrium does not give, from the F_INIT pf_active at t0, so du = 0 until the controller changes them); `f_init`: coil currents and voltages (and Ip / current centre when present) from the F_INIT data at t0, passive currents from the model. Each step writes Ip and `current_centre.r/z` into `equilibrium_o_i`; `boundary.geometric_axis` stays the reference pass-through.

**How to run.**
- MATLAB: `[model, x0] = configure_vertical_instability_model(dt, extract_equilibrium_data(uri, time, opts), em_coupling, pf_active, pf_passive)`; `run_vertical_instability` (author's open-loop demo, default entry). `run_plasmaless_imas` / `run_plasmaless_model` with the environment variables `PLASMALESS_VERTICAL_EQUILIBRIUM_URI` (+ optional `PLASMALESS_VERTICAL_EQUILIBRIUM_TIME`, `PLASMALESS_VERTICAL_PASSIVE_CURRENTS`): results `plasmaless_vertical_results.mat` / `plasmaless_matlab_vertical_results.mat` (`Y` = Ip, rc, zc; absolute values in the IMAS run, deviations in the Simulink/MATLAB run). The 28 s pulse train makes the open-loop unstable model diverge; shorten `tend` to look at it.
- MUSCLE3 standalone: `PLASMALESS_VERTICAL_EQUILIBRIUM_URI=/work/projects/pds/pds-scenarios/105073/data/in PLASMALESS_VERTICAL_EQUILIBRIUM_TIME=84.4176 PLASMALESS_VERTICAL_PASSIVE_CURRENTS=zero PLASMALESS_M3_DT=0.001 PLASMALESS_M3_T_END=1.0 muscle3/run_plasmaless_muscle3.sh` (results `plasmaless_muscle3_vertical_results.mat`). Run on 2026-10-15: exit 0, 1000 exchanges, zc constant until the first CS3U pulse edge at 0.5 s (-20 V deviation), then 0.97 mm away from `zc0` at 1.0 s; currents and zc equal a direct `plasmaless_step_ids` loop (max differences 4e-12 A, 1e-16 m).
- Tests: `tests/test_vertical_instability.m` (author's, 2 tests) and `tests/test_plasmaless_vertical.m` (7 tests: build from the DD4 entry with zero passive currents and from the DD3 default, step loop vs direct recursion, Ip/rc/zc round trip through the equilibrium IDS, explicit zero-passive option, entry-path conversion, circuit and linear models unchanged).

**Limitations.**
- Ip fixed, no radial motion, linear around one equilibrium; only the vacuum field of the equilibrium conductor currents enters the force-balance stiffness (author's README).
- Passive currents: the PDS scenario entries have none; `zero` is an explicit assumption that changes `Ic0` and so the stiffness (DD3 105108: growth rate 3.172 1/s with the entry's passive currents, 3.166 1/s with zero).
- The PDS PCS controller does not close the vertical loop: in `plasmaless_controller` the bridge drops the RZIp `VS1_volt_cmd` and commands VS3 only with the `V = R*I` fallback (magnetic-control reference `pcs_pcssp_in_pds.md` §3 and §5; `controllers/KCURR_RZIp/muscle_NICE_input.m:29-36` in the IMAS-PDS `feature/plasmaless_controller` worktree at 32849ec), and the controller reads Ip and `boundary.geometric_axis`, not `current_centre`. A PDS closed-loop run with this model will diverge (growth rate 3-8 1/s) until a vertical-stabilisation controller is wired.

## DD4 machine description

The `em_coupling` is the DD4 IDS with `coupling_matrix` entries named `mutual_active_active`, `mutual_passive_active`, `mutual_passive_passive`, `mutual_grid_active` and `mutual_grid_passive`; `get_coupling_matrix` checks the name, the quantity (magnetic flux) and that `rows_uri`/`columns_uri` refer to the expected nodes. The files in `data/md_dd4` (not pulse-specific; default of `run_plasmaless_imas.m`, `run_plasmaless_model.m`, the plots and the MUSCLE3 actor) are made from the DD3 IDSs by `tools/convert_md_dd3_to_dd4.py` (see Installation). DD3 IDSs are not read by the model: the DD4 IMAS-MATLAB API returns an empty `coupling_matrix` for a DD3 entry, so DD3 data must be converted first.

Limitations:
- The inputs of the step interface are coil voltages only. There is no supply/circuit-to-coil mapping in the circuit model, so CS1U/CS1L and VS3U/VS3L are treated as independent coil voltages, as in the original model (the linear model sums CS1U and CS1L into one state through its B matrix, and its generator inputs are mapped to the coil voltages by the descriptor's `input_map`; the VS1 generator has no coil and is not drivable).
- The circuit model returns currents only; the flux on the DINA grid is computed in `plot_plasmaless_model.m` from `Mxi`. The linear model returns its observables too, but no flux map, and `plot_plasmaless_model.m` (102 loops, DINA grid) does not apply to it.
- The linear model is a deviation model around one linearisation point; `Dist2ndSep` has no DD4 equilibrium field and is not sent to the controller.

## PDS integration

The model runs in IMAS-PDS as the `plasmaless` actor of the `plasmaless_controller` workflow (branch `feature/plasmaless_controller`), closed on the PCS magnetic controller. How to run it and the standalone checks: `instructions.txt`. How the integration was done: `pds_integration.txt`.

## Licence

This software is open source. ITER-related intellectual property remains &copy; 2026 ITER Organization.

## Contacts

Originally developed in 2026 by Federico Pesamosca (federico.pesamosca@iter.org)

IMAS interfacing and MUSCLE3 adaptation by Mireille Schneider (mireille.schneider@iter.org)
