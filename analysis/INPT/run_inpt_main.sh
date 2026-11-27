#!/bin/bash
# Part of the public analysis code for resident inpatient days in Korean districts, 2011-2024 (MIT licence). Inputs are read from environment variables; see README.md.
source ${ST_PROJECT}/env.sh
RS="$SC/rslot.sh"; RCLOG="$L/run_inpt_rc.log"; LG="$L/logs"; mkdir -p "$LG"
rcl(){ echo "$(date '+%F %T') rc=$2 $1" >> "$RCLOG"; }
chmod +x "$RS"; bash -n "$0" || { rcl precheck_bash 1; exit 1; }
Rscript -e "invisible(parse('$SC/81_fit_inla_inpt.R'))" || { rcl precheck_parse 1; exit 1; }
rcl precheck 0
export TAG=INPT ST_VER=std3 ST_THREADS=1 ST_GRAPH="$L/graph_593.graph" SK_RSLOTS=2
cd "$SC"; ST_MV_DEFER=1 ST_MV_PARTS=4 $RS Rscript 81_fit_inla_inpt.R INPT > "$LG/81_main_v1.log" 2>&1; rc=$?; rcl 81_main_v1 $rc
echo DONE rc=$rc; exit $rc
