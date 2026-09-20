#!/bin/bash
#SBATCH -C cpu
#SBATCH -q shared
#SBATCH --module=cvmfs
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -J mcp_geant
#SBATCH -o output/slurm-%x_%A_%a.out
#SBATCH -e output/slurm-%x_%A_%a.err

module load python
PYTHON=$(command -v python3)
GCC14_LIBDIR=$(dirname "$(g++ -print-file-name=libstdc++.so.6)")

set -euo pipefail

PROJECT=${PROJECT:-${SLURM_SUBMIT_DIR:?}}
SAMPLER="$PROJECT/scripts/sample.py"
FILTER="$PROJECT/scripts/filter.py"
GEANT=${GEANT:-$PROJECT/build/g4Sim}
JOB_LIST=${JOB_LIST:-job.csv}
JOB_ID=${SLURM_ARRAY_TASK_ID:?Submit this script as a Slurm array}
export MCP_FAST_BLIPS=${MCP_FAST_BLIPS:-1}

case "$MCP_FAST_BLIPS" in
    0|1) ;;
    *) echo "MCP_FAST_BLIPS must be 0 or 1" >&2; exit 2 ;;
esac

cd "$PROJECT"

# Resolve a relative profile path before Geant4 changes to its temporary directory.
if [[ "$MCP_FAST_BLIPS" == 1 && -n "${MCP_BLIP_PROFILE:-}" ]]; then
    case "$MCP_BLIP_PROFILE" in
        /*) ;;
        *) MCP_BLIP_PROFILE="$PWD/$MCP_BLIP_PROFILE" ;;
    esac
    export MCP_BLIP_PROFILE
    if [[ ! -r "$MCP_BLIP_PROFILE" ]]; then
        echo "Cannot read MCP_BLIP_PROFILE: $MCP_BLIP_PROFILE" >&2
        exit 2
    fi
fi

row=$("$PYTHON" - "$JOB_LIST" "$JOB_ID" <<'PY'
import csv
import sys

with open(sys.argv[1]) as source:
    for row in csv.DictReader(source):
        if row["job_id"] == sys.argv[2]:
            print(f'{row["emitter"]}|{row["mcp_mass_GeV"]}')
            break
    else:
        raise SystemExit(f"Job ID {sys.argv[2]} not found")
PY
)

IFS='|' read -r emitter mass <<< "$row"

case "$emitter" in
    111) product=pi0 ;;
    221) product=eta ;;
    331) product=etap ;;
    113) product=rho0 ;;
    223) product=omega ;;
    333) product=phi ;;
    443) product=jpsi ;;
    0)   product=dy ;;
    *)
        echo "Unknown emitter: $emitter" >&2
        exit 2
        ;;
esac

mass_tag=$(printf "%.6f" "$mass" | tr '.' 'p')
WORK="/pscratch/sd/h/hogcai/mcp-geant4/tmp/${SLURM_ARRAY_JOB_ID}_${JOB_ID}"
INPUT_CSV="$WORK/mcp_input.csv"
RUN_MAC="$WORK/run.mac"
FILTERED_ROOT="$WORK/selected.root"
if [[ "$MCP_FAST_BLIPS" == 1 ]]; then
    OUTPUT_DIR=${OUTPUT_DIR:-/pscratch/sd/h/hogcai/mcp-geant4/output-fast}
    OUTPUT="$OUTPUT_DIR/mcp_fast_${product}_${mass_tag}.root"
else
    OUTPUT_DIR=${OUTPUT_DIR:-/pscratch/sd/h/hogcai/mcp-geant4/output}
    OUTPUT="$OUTPUT_DIR/mcp_${product}_${mass_tag}.root"
fi

mkdir -p "$WORK" "$OUTPUT_DIR"

if [[ -s "$OUTPUT" && "${OVERWRITE:-0}" != 1 ]]; then
    echo "Output exists, skipping: $OUTPUT"
    exit 0
fi

"$PYTHON" "$SAMPLER" "$mass" "$emitter" "$INPUT_CSV"

N_EVENTS=$(($(wc -l < "$INPUT_CSV") - 1))
printf '/run/beamOn %d\n' "$N_EVENTS" > "$RUN_MAC"

export MCP_INPUT_CSV="$INPUT_CSV"
export MCP_MASS_GEV="$mass"

(
    set +e
    set +u
    unset PERL5LIB PERLLIB
    source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
    setup geant4 v4_11_4_p02 -q e28:prof
    GEANT_SETUP_STATUS=$?
    export LD_LIBRARY_PATH="$GCC14_LIBDIR:${LD_LIBRARY_PATH:-}"
    set -euo pipefail

    if [[ "$GEANT_SETUP_STATUS" -ne 0 ]]; then
        echo "Geant4 setup failed" >&2
        exit "$GEANT_SETUP_STATUS"
    fi

    cd "$WORK"
    "$GEANT" "$RUN_MAC" 1
)

if [[ "$MCP_FAST_BLIPS" == 1 ]]; then
    if [[ ! -s "$WORK/mcp_fast.root" ]]; then
        echo "Missing fast output; check that GEANT points to the rebuilt executable." >&2
        exit 1
    fi
    mv "$WORK/mcp_fast.root" "$OUTPUT"
else
    "$PYTHON" "$FILTER" "$WORK/mcp_02.root" "$FILTERED_ROOT"
    mv "$FILTERED_ROOT" "$OUTPUT"
    rm -f "$WORK/mcp_02.root"
fi
rm -f "$INPUT_CSV" "$RUN_MAC"
rmdir "$WORK" 2>/dev/null || true

echo "emitter=$emitter"
echo "product=$product"
echo "mass_GeV=$mass"
echo "events=$N_EVENTS"
echo "fast_blips=$MCP_FAST_BLIPS"
if [[ "$MCP_FAST_BLIPS" == 1 ]]; then
    echo "blip_profile=${MCP_BLIP_PROFILE:-constant calibrated rate from N1/N0}"
fi
echo "output=$OUTPUT"
