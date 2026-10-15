# MUSCLE3 actor: plasmaless coil+vessel model

`muscle_plasmaless_actor.m` wraps the IMAS interface of the model (`../imas_model/`) as a
MUSCLE3 (libmuscle 0.10, Python in-process from MATLAB) actor. It stands in for NICE
direct evolutive (`nice_evo_rd`) in controller tests. In IMAS-PDS it is the program
`plasmaless` (`workflows/lib/local_programs.ymmsl`), used by the workflow
`plasmaless_controller`; PDS installs this repository into
`local_installs/plasma-tokamak-circuits` with `setup_files/setup_plasmaless.sh`.

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
- `linear_model_dir` (optional, default `''` = circuit model): absolute path of a
  directory holding `linear_model.json` and the matrices of an external linear model
  (`imas_model/configure_linear_model.m`, repository README "External linear model").
  The machine description is still read for the coil names and order.
- `vertical_equilibrium_uri` (optional, default `''` = not used): IMAS URI or entry
  directory of an equilibrium with `profiles_2d`; the actor then runs the
  vertical-instability model (`../vertical_instability/`, repository README
  "Vertical-instability model"). Setting it together with `linear_model_dir` is the
  error `plasmaless_actor:settings`.
- `vertical_equilibrium_time` (optional, default t0): time of that equilibrium.
- `vertical_passive_currents` (optional, default `equilibrium`): `equilibrium` (the
  entry's passive currents, error `vertical:passiveCurrents` when it has none) or `zero`.
- `operating_point` (optional, linear and vertical models, default `descriptor`): which
  source wins for the operating point, `descriptor` (the model's own reference: the
  descriptor of a linear model, the equilibrium of the vertical model; the value name is
  kept for compatibility) or `f_init` (see below); any other value is the error
  `plasmaless_actor:settings`.

## External linear model: operating point

The linear model's states, inputs and outputs are deviations from its linearisation
point. The actor starts at t0 with x = 0 and builds the operating point from two
sources: the descriptor (`operating_point`: coil currents, voltages, outputs) and the
F_INIT data at t0 (`op.coil_current` and `op.voltage` = currents and voltages of the
F_INIT pf_active slice closest to t0 in model coil order, `op.y` = the observables found
in the F_INIT equilibrium slice closest to t0, `plasmaless_observables_from_equilibrium`).
Setting `operating_point`:

- `descriptor` (default): descriptor values win, F_INIT fills what the descriptor does
  not give. The absolute coil currents and Ip sent to the controller then start at the
  linearisation point, which differs from the F_INIT scenario values: the controller sees
  a step error at t0 (for the received model, Ip +15.06 MA against the scenario's).
- `f_init`: the coil currents, voltages and the observables found in the F_INIT data win;
  descriptor values are used only for what F_INIT does not have. The outputs start on
  the scenario, and the model's deviations are added to it.

An observable known by neither source stays unknown and is not sent. Each step uses
du = V - `op.voltage`, so the operating point is held until the controller changes the
voltages; the actor prints one line naming the setting, how many coil entries came from
F_INIT, the observables taken from the F_INIT equilibrium, those sent as absolute
values and those not sent.
Outputs: `pf_active_o_i` = `op.coil_current` + deviations; `pf_passive_o_i` has as many
loops as the model has passive states (4 for the received model, `PS_1, PS_2, VV_1,
VV_2`), named as in the
descriptor; `equilibrium_o_i` = the reference pass-through plus the known observables
(`plasmaless_observables_to_equilibrium`: a known Ip overrides the pass-through ip;
magnetic axis, X-point (`contour_tree.node{1}`, `critical_type` 1) and gaps are added;
`boundary.geometric_axis` stays the reference pass-through because the model does not
give it; `Dist2ndSep` has no DD4 field and is not sent), with an `ids_properties.comment`
naming the model and the observables it provides.

## Vertical-instability model

States are the deviations from the equilibrium currents `Ic0`; x = 0 at t0. Operating
point: `op.coil_current` and `op.loop_current` = `Ic0`, `op.y` = Ip0, rc0, zc0 (all
given by the equilibrium); `op.voltage` is not given (the equilibrium has no voltage
reference), so with `descriptor` the voltages come from the F_INIT pf_active at t0 and
du = 0 until the controller changes them. With `f_init` the coil currents and voltages
(and Ip / `current_centre` when the F_INIT equilibrium has them) come from F_INIT; the
passive currents always stay the model's. `pf_passive_o_i` has the 102
machine-description loops. `equilibrium_o_i` carries Ip (`global_quantities.ip`, fixed)
and the current centre (`global_quantities.current_centre.r` fixed, `.z` = zc0 +
zIp/Ip0); `boundary.geometric_axis` stays the reference pass-through. The PCS controller
in PDS reads Ip and the geometric axis, not `current_centre`, and does not command VS1/VS3
(repository README, Limitations): a closed-loop PDS run with this model diverges until a
vertical-stabilisation controller is wired. The actor prints the growth rate, Ip0, rc0,
zc0 and the LCFS source at build time.

Standalone run (2026-10-15): `run_plasmaless_muscle3.sh` with
`PLASMALESS_VERTICAL_EQUILIBRIUM_URI=/work/projects/pds/pds-scenarios/105073/data/in`,
`PLASMALESS_VERTICAL_EQUILIBRIUM_TIME=84.4176`, `PLASMALESS_VERTICAL_PASSIVE_CURRENTS=zero`,
`PLASMALESS_M3_DT=0.001`, `PLASMALESS_M3_T_END=1.0` (the script writes them as
`plasmaless.vertical_*`, `dt`, `t_interval`, `t_end` into `paths.ymmsl`, which wins over
`plasmaless_standalone.ymmsl`): exit 0, 1000 exchanges, growth rate 7.6167 1/s; zc stays
at zc0 = 0.590943 m until the CS3U pulse edge at 0.5 s (-20 V from the F_INIT voltage),
then moves 0.97 mm by 1.0 s; currents and zc equal a direct `plasmaless_step_ids` loop
(4e-12 A, 1e-16 m). Results `plasmaless_muscle3_vertical_results.mat`; `Y` has the
columns `rc`, `zc` in addition (the driver reads `current_centre` too).

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

Slurm run with the PCS controller and the external linear model (2026-10-14), shot 105073,
case `plasmaless_controller_105073` (override `cases/overrides/plasmaless_controller_105073.ymmsl`:
`source.t_min` 81.0 / `t_max` 100.0 on the inverse `nice_out` of
`metis_nice_inverse_from_dina_105073_20261014_102619`, t0 = 81.4576 s, `plasmaless.t_end`
83.4576, dt = `t_interval` 0.005, `plasmaless.linear_model_dir` =
`linearModelSimple_105073_t84.417598`, `plasmaless.operating_point` descriptor; `PLASMALESS_REPO`
exported at submission because the PDS local install predates the linear-model code), run dir
`cases/runs/plasmaless_controller_105073_20261014_103330`: finished without error, 6 instances
exit 0, 400 exchanges (81.4626 .. 83.4576 s), 17 states / 14 coils / 4 passive loops / 9
observables, operating point from the descriptor (0/0 coil entries from F_INIT), one expected
`plasmaless:inputUndriven` warning (VS1). Integration test only; the physics (CS/PF commands at
the clip from the first step, VS3 currents to -100 kA, rmag +16 cm) is discussed in
`../pds_integration.txt` 4.7.

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
`run_plasmaless_imas.m` layout (`time`, `Ia`, `Ie`, `Y`, `coil_names`, `loop_names`,
`y_names`; `Y` holds the observables read back from `equilibrium_in`, NaN where absent).
With a linear model directory as first argument (`./run_plasmaless_muscle3.sh <dir>`)
the script sets `plasmaless.linear_model_dir` and writes
`plasmaless_muscle3_linear_results.mat`. With the default `operating_point: descriptor`
and the received model, the operating point is the descriptor's (coil currents of the
linearisation point, voltages 0, rmag/zmag/Ip known) and the saved `Ia` and the
`rmag`/`zmag`/`Ip` columns of `Y` are absolute values (operating point + deviation), the
same as `run_plasmaless_imas.m` with dt 0.01 (except the first row, t0, which the driver
fills with its own cold-start 0 and F_INIT ip 0); the driver's cold-start pf_active
(currents 0, pulse-train voltages at t0) and F_INIT equilibrium (ip 0) are not used
because the descriptor gives everything they carry. With `operating_point: f_init`
(uncomment it in `plasmaless_standalone.ymmsl` or add it to `paths.ymmsl`) the
operating point becomes the driver's: currents 0, `op.voltage` = u(t0) (10 V on CS3U, so
the first pulse is seen as 0 V and the later ones as -10 V), Ip known as 0, rmag/zmag
from the descriptor. To plot it with `plot_plasmaless_model.m`, load it in
MATLAB from the repository root and build the variables that script uses:
`load('plasmaless_muscle3_results.mat'); out.simout = timeseries([Ia Ie], time);`
`tstart = time(1); tend = time(end); TS = timeseries_plasmaless_model(tstart,tend,time(2)-time(1)); plot_plasmaless_model`.
With the same dt, the currents equal those of `run_plasmaless_imas.m`.
