#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "${SCRIPT_DIR}"

mkdir -p figs logs

SWE_REFERENCE_FILE=${SWE_REFERENCE_FILE:-/share/home/dq135/openbench/Reference/Grid/HigRes/Snow/Snow_Water_Equivalent/ERA5_Land/ERA5_Land_SWE_2001-2017.nc}
ERA5_SNOW_COVER_REFERENCE_FILE=${ERA5_SNOW_COVER_REFERENCE_FILE:-/share/home/dq135/openbench/Reference/Grid/HigRes/Snow/Snow_Cover_Fraction/ERA5_Land/ERA5_Land_SnowCover_2001-2017.nc}
MODIS_SCF_REFERENCE_FILE=${MODIS_SCF_REFERENCE_FILE:-/share/home/dq135/openbench/Reference/Grid/HigRes/Snow/Snow_Cover_Fraction/MODIS_MOD10CM/MOD10CM_SCF_2001_2017.nc}
SNOW_COVER_THRESHOLD=${SNOW_COVER_THRESHOLD:-0.01}

for required_file in \
    "${SWE_REFERENCE_FILE}" \
    "${ERA5_SNOW_COVER_REFERENCE_FILE}" \
    "${MODIS_SCF_REFERENCE_FILE}"
do
    if [[ ! -f "${required_file}" ]]; then
        echo "[ERROR] Required snow reference is missing: ${required_file}" >&2
        exit 1
    fi
done

run_logged() {
    local log_file=$1
    shift
    echo "[RUN] $*"
    "$@" 2>&1 | tee "${log_file}"
}

timestamp=$(date +%Y%m%d_%H%M%S)

snow_args=(
    --swe-reference-file "${SWE_REFERENCE_FILE}"
    --swe-reference-variable sd
    --era5-snow-cover-reference-file "${ERA5_SNOW_COVER_REFERENCE_FILE}"
    --era5-snow-cover-reference-variable snowc
    --scf-reference-file "${MODIS_SCF_REFERENCE_FILE}"
    --scf-reference-variable Snow_Cover_Monthly_CMG
    --snow-cover-threshold "${SNOW_COVER_THRESHOLD}"
)

echo "=== Input validation ==="
run_logged "logs/land_snow_check_${timestamp}.log" \
    python land_snow_response.py "${snow_args[@]}" --check-only
run_logged "logs/land_et_check_${timestamp}.log" \
    python land_evapotranspiration_response.py --check-only
run_logged "logs/land_runoff_check_${timestamp}.log" \
    python land_runoff_response.py --check-only
run_logged "logs/land_infiltration_check_${timestamp}.log" \
    python land_water_balance_response.py --check-only

echo "=== Figure generation ==="
run_logged "logs/land_snow_${timestamp}.log" \
    python land_snow_response.py "${snow_args[@]}"
run_logged "logs/land_et_${timestamp}.log" \
    python land_evapotranspiration_response.py
run_logged "logs/land_runoff_${timestamp}.log" \
    python land_runoff_response.py
run_logged "logs/land_infiltration_${timestamp}.log" \
    python land_water_balance_response.py

echo "[DONE] All supplementary land diagnostics completed successfully."
echo "[DONE] Figures: ${SCRIPT_DIR}/figs"
echo "[DONE] Logs use timestamp: ${timestamp}"
