#!/bin/bash
# Standalone MUSCLE3 test of the plasmaless actor: voltage_driver.py + actor,
# settings in plasmaless_standalone.ymmsl (dt 0.01, t = 0 .. 28 s: ~3 min).
# Optional first argument: directory of an external linear model (absolute path,
# linear_model.json + matrices); the actor then runs that model instead of the
# circuit model (setting plasmaless.linear_model_dir) and the results go to
# plasmaless_muscle3_linear_results.mat.
# Vertical-instability model instead (vertical_instability/): export
#   PLASMALESS_VERTICAL_EQUILIBRIUM_URI   IMAS URI or entry directory of the equilibrium
#   PLASMALESS_VERTICAL_EQUILIBRIUM_TIME  (optional) its time, default t0 (= 0 here)
#   PLASMALESS_VERTICAL_PASSIVE_CURRENTS  (optional) equilibrium (default) | zero
# (settings plasmaless.vertical_equilibrium_*, vertical_passive_currents); results
# in plasmaless_muscle3_vertical_results.mat. Not together with a linear model.
# Optional PLASMALESS_M3_DT and PLASMALESS_M3_T_END override dt (= t_interval) and
# t_end of plasmaless_standalone.ymmsl (the vertical model is unstable: use a short
# t_end, e.g. PLASMALESS_M3_DT=0.001 PLASMALESS_M3_T_END=1.0).
# Optional PLASMALESS_M3_COMMAND (supply | coil) sets driver.command: the driver sends
# the pulse train as the 14 power-supply voltages (pf_active.supply, default) or as
# the 14 coil voltages (the PDS controller's format; the VS1 pulse is then lost).
# Uses 2 cores (driver + one single-threaded MATLAB); on the shared login node,
# run it in Slurm for longer cases, e.g.
#   srun --partition=all_debug --time=00:30:00 --cpus-per-task=2 --mem=8G ./run_plasmaless_muscle3.sh

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PLASMALESS_REPO="$(dirname "$here")"
export OMP_NUM_THREADS=1
linear_model_dir="${1:-}"
if [ -n "$linear_model_dir" ] && [ ! -f "$linear_model_dir/linear_model.json" ]; then
    echo "run_plasmaless_muscle3.sh: $linear_model_dir has no linear_model.json" >&2
    exit 2
fi
vertical_uri="${PLASMALESS_VERTICAL_EQUILIBRIUM_URI:-}"
if [ -n "$linear_model_dir" ] && [ -n "$vertical_uri" ]; then
    echo "run_plasmaless_muscle3.sh: a linear model and PLASMALESS_VERTICAL_EQUILIBRIUM_URI are both given" >&2
    exit 2
fi
# the run directory only holds the MUSCLE3 logs; the results go to the repo root
run_dir=/scratch/users/$USER/plasmaless/m3_standalone_$(date +%Y%m%d_%H%M%S)
mkdir -p "$run_dir"
if [ -n "$linear_model_dir" ]; then
    results_file="$PLASMALESS_REPO/plasmaless_muscle3_linear_results.mat"
elif [ -n "$vertical_uri" ]; then
    results_file="$PLASMALESS_REPO/plasmaless_muscle3_vertical_results.mat"
else
    results_file="$PLASMALESS_REPO/plasmaless_muscle3_results.mat"
fi
# absolute paths for the driver (a relative path would resolve in its instance work dir)
cat > "$run_dir/paths.ymmsl" <<EOF
ymmsl_version: v0.2
settings:
  driver.md_uri: "imas:hdf5?path=$PLASMALESS_REPO/data/md_dd4"
  driver.results_file: "$results_file"
EOF
if [ -n "$linear_model_dir" ]; then
    echo "  plasmaless.linear_model_dir: \"$linear_model_dir\"" >> "$run_dir/paths.ymmsl"
fi
if [ -n "$vertical_uri" ]; then
    echo "  plasmaless.vertical_equilibrium_uri: \"$vertical_uri\"" >> "$run_dir/paths.ymmsl"
    if [ -n "${PLASMALESS_VERTICAL_EQUILIBRIUM_TIME:-}" ]; then
        echo "  plasmaless.vertical_equilibrium_time: $PLASMALESS_VERTICAL_EQUILIBRIUM_TIME" >> "$run_dir/paths.ymmsl"
    fi
    if [ -n "${PLASMALESS_VERTICAL_PASSIVE_CURRENTS:-}" ]; then
        echo "  plasmaless.vertical_passive_currents: $PLASMALESS_VERTICAL_PASSIVE_CURRENTS" >> "$run_dir/paths.ymmsl"
    fi
fi
# dt, t_interval and t_end overrides (a later settings file wins)
if [ -n "${PLASMALESS_M3_DT:-}" ]; then
    echo "  dt: $PLASMALESS_M3_DT" >> "$run_dir/paths.ymmsl"
    echo "  t_interval: $PLASMALESS_M3_DT" >> "$run_dir/paths.ymmsl"
fi
if [ -n "${PLASMALESS_M3_T_END:-}" ]; then
    echo "  t_end: $PLASMALESS_M3_T_END" >> "$run_dir/paths.ymmsl"
fi
# voltage command format of the driver (default supply, in plasmaless_standalone.ymmsl)
if [ -n "${PLASMALESS_M3_COMMAND:-}" ]; then
    case "$PLASMALESS_M3_COMMAND" in
        supply|coil) ;;
        *) echo "run_plasmaless_muscle3.sh: PLASMALESS_M3_COMMAND must be supply or coil (got '$PLASMALESS_M3_COMMAND')" >&2
           exit 2 ;;
    esac
    echo "  driver.command: $PLASMALESS_M3_COMMAND" >> "$run_dir/paths.ymmsl"
fi

module use /work/projects/pds/modules/all
module load IMAS-MUSCLE3/1.0.0-intel-2025b-pds

muscle_manager --start-all --run-dir "$run_dir" "$here/plasmaless_standalone.ymmsl" "$run_dir/paths.ymmsl"
rc=$?

echo "muscle_manager exit code: $rc"
echo "Logs:    $run_dir (muscle3_manager.log, instances/*/stdout.txt, stderr.txt)"
echo "Results: $results_file (time, Ia, Ie, Y, Vps, coil_names, loop_names, y_names, input_names)"
exit $rc
