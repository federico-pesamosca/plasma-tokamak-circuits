#!/bin/bash
# Standalone MUSCLE3 test of the plasmaless actor: voltage_driver.py + actor,
# settings in plasmaless_standalone.ymmsl (dt 0.01, t = 0 .. 28 s: ~3 min).
# Uses 2 cores (driver + one single-threaded MATLAB); on the shared login node,
# run it in Slurm for longer cases, e.g.
#   srun --partition=all_debug --time=00:30:00 --cpus-per-task=2 --mem=8G ./run_plasmaless_muscle3.sh

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PLASMALESS_REPO="$(dirname "$here")"
export OMP_NUM_THREADS=1
# the run directory only holds the MUSCLE3 logs; the results go to the repo root
run_dir=/scratch/users/$USER/plasmaless/m3_standalone_$(date +%Y%m%d_%H%M%S)
mkdir -p "$run_dir"
results_file="$PLASMALESS_REPO/plasmaless_muscle3_results.mat"
# absolute paths for the driver (a relative path would resolve in its instance work dir)
cat > "$run_dir/paths.ymmsl" <<EOF
ymmsl_version: v0.2
settings:
  driver.md_uri: "imas:hdf5?path=$PLASMALESS_REPO/data/md_dd4"
  driver.results_file: "$results_file"
EOF

module use /work/projects/pds/modules/all
module load IMAS-MUSCLE3/1.0.0-intel-2025b-pds

muscle_manager --start-all --run-dir "$run_dir" "$here/plasmaless_standalone.ymmsl" "$run_dir/paths.ymmsl"
rc=$?

echo "muscle_manager exit code: $rc"
echo "Logs:    $run_dir (muscle3_manager.log, instances/*/stdout.txt, stderr.txt)"
echo "Results: $results_file (time, Ia, Ie, coil_names, loop_names)"
exit $rc
