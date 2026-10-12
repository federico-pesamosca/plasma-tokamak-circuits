# MUSCLE3 actor: plasmaless coil+vessel model

`muscle_plasmaless_actor.m` wraps the IMAS interface of the model (`../imas_model/`) as a
MUSCLE3 (libmuscle 0.10, Python in-process from MATLAB) actor. It stands in for NICE
direct evolutive (`nice_evo_rd`) in controller tests. In IMAS-PDS it is the program
`plasmaless` (`workflows/lib/local_programs.ymmsl`), used by the workflow
`plasmaless_controller`; PDS installs this repository into
`local_installs/plasmaless-tokamak-circuits` with `setup_files/setup_plasmaless.sh`.

## Ports

| Operator | Port | Content |
|---|---|---|
| F_INIT | `equilibrium_f_init` | multi-time equilibrium (reference trajectory); its message timestamp is the start time t0 |
| F_INIT | `pf_active_f_init` | multi-time pf_active, `homogeneous_time` 1, coil currents and voltages |
| S | `pf_active_s` | one-slice pf_active with the coil voltages at t (coil names not needed) |
| O_I | `equilibrium_o_i` | one-slice equilibrium at t (see below) |
| O_I | `pf_active_o_i` | one-slice coil currents at t (and the voltages applied over the step) |
| O_I | `pf_passive_o_i` | one-slice passive loop (vessel) currents at t |

Flag `KEEPS_NO_STATE_FOR_NEXT_USE`.

## Settings

- `dt` (mandatory): model time step [s] (discretisation of the model).
- `t_interval` (mandatory): exchange period [s]; must equal `dt`.
- `t_end` (mandatory): last output time [s]; t0 < `t_end` <= last time of the F_INIT
  equilibrium.
- `md_uri` (optional): IMAS URI of the DD4 machine description (`em_coupling`,
  `pf_active`, `pf_passive`). Default: `data/md_dd4` of this repository (derived from
  the actor file location), made with `tools/convert_md_dd3_to_dd4.py`.

## Timing

As `nice_imas_evo_rd_muscle3` (NICE 3.0.0.dev446, `main_imas_evo_rd_muscle3.cc`):
t0 = timestamp of `equilibrium_f_init`; the F_INIT pf_active is sliced at t0 (closest
sample). The first step uses those scenario voltages; each step gives the outputs at
t = t_prev + dt, sent on the three O_I ports with timestamp t and next_timestamp t + dt,
or None when t + dt > t_end + 1e-9 (then the actor stops). Otherwise the actor receives
`pf_active_s` (expected at timestamp t) and uses its voltages for the next step. The stop
decision comes from `t_end` only: the PDS controller always sends next_timestamp None.

## Assumptions

- Vessel (pf_passive loop) currents are zero at t0: the scenario data has none.
- The coils of the F_INIT pf_active and of the machine description are the same coils,
  in the same order (the controller works by index): the short names (CS3U ... PF6,
  VS3U, VS3L) must match exactly, otherwise the actor fails.
- The coil resistances come from the model machine description; the resistances sent
  by the controller are ignored.
- `equilibrium_o_i` has no plasma: `global_quantities.ip` and
  `boundary.geometric_axis.r/z` are the F_INIT equilibrium values linearly interpolated
  at t (so the controller's Ip/R/Z errors are zero), `code.output_flag` = 0, and
  `ids_properties.comment` says so.
- After the actor returns, the PDS launcher ends MATLAB with `os._exit(0)`: MATLAB
  R2025b with IMAS-MATLAB 5.6.0 segfaults on exit after reading an `imas:hdf5` entry.

## Test

Tested with a dummy controller in place of the PCS one (shot 105073 F_INIT, t0 = 84.3176 s,
with the former VSU/VSL name mapping,
20 steps): 20 exchanges, first output at t0 + 0.01, last with next_timestamp None, and
coil currents identical to direct `plasmaless_step_ids` calls.

Earlier test (historical, when the PDS workflow still ran `waveform_editor` inside the
workflow) for shot 105084: source and waveform_editor, dummy controller,
t0 = 136.2276 s, dt 0.005, t_end = t0 + 0.02: 4 exchanges; the F_INIT equilibrium and
pf_active received by the controller equal the NICE inverse output (`nice_out`) exactly.

In the current PDS `plasmaless_controller` workflow (IMAS-PDS,
`workflows/plasmaless_controller/`), the Waveform-Editor runs once at case creation
(`preprocess.sh`, writes `<case>/preprocess/initial_state` from `nice_out`); `source`
reads that entry and sends the F_INIT equilibrium and pf_active directly to
`magnetic_controller` and `plasmaless`, and a `recorder_plasmaless` component
(`visualization/plasmaless_currents.py`) shows live coil/passive-current figures in
m3dash. Slurm run with the PCS controller, shot 105084, t_end 160
(`cases/runs/plasmaless_controller_105084_20261008_200304`): 4754 exchanges, finished
without error.

## Standalone test

`plasmaless_standalone.ymmsl` couples the actor to `voltage_driver.py` only (no PDS, no
PCS controller): the driver sends a no-plasma F_INIT equilibrium and a cold-start
pf_active (machine-description coils, currents 0), then answers every output with the
voltage pulse train of `timeseries_plasmaless_model.m`, as `run_plasmaless_imas.m`
applies it (u(t_k) over [t_k, t_k+dt]). Settings: `dt` = `t_interval` = 0.01,
`t_end` = 28; the run script adds `driver.md_uri` and `driver.results_file` with absolute
repository paths (`paths.ymmsl` in the run directory); the actor uses its `md_uri`
default. Run (2 cores, ~3 min):

    ./run_plasmaless_muscle3.sh

Run directory (MUSCLE3 logs only, printed at the end):
`/scratch/users/$USER/plasmaless/m3_standalone_<date>`. The results file
`plasmaless_muscle3_results.mat` (repository root, ignored by git) has the
`run_plasmaless_imas.m` layout. To plot it with `plot_plasmaless_model.m`, load it in
MATLAB from the repository root and build the variables that script uses:
`load('plasmaless_muscle3_results.mat'); out.simout = timeseries([Ia Ie], time);`
`tstart = time(1); tend = time(end); TS = timeseries_plasmaless_model(tstart,tend,time(2)-time(1)); plot_plasmaless_model`.
With the same dt, the currents equal those of `run_plasmaless_imas.m`.
