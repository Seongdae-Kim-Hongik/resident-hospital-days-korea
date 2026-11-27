#!/bin/bash
# Part of the public analysis code for resident inpatient days in Korean districts, 2011-2024 (MIT licence). Inputs are read from environment variables; see README.md.
source ${ST_PROJECT}/env.sh
T=INPT; RS="$SC/rslot.sh"; GR="$L/graph_593.graph"; RCLOG="$L/run_inpt_rc.log"; LG="$L/logs"; mkdir -p "$LG"
rcl(){ echo "$(date '+%F %T') rc=$2 $1" >> "$RCLOG"; }
wl(){ while :; do l=$(sysctl -n vm.loadavg | awk '{print $2}'); n=$(ls -d ${TMPDIR:-.}/rslots/s* 2>/dev/null | wc -l | tr -d ' '); awk -v l="$l" -v n="$n" 'BEGIN{exit !(l<=12 || (n<=2 && l<=30))}' && break; sleep 60; done; }
PIDS=(); NAMES=(); bgadd(){ PIDS+=("$1"); NAMES+=("$2"); }
waitall(){ local i rc; for i in "${!PIDS[@]}"; do wait "${PIDS[$i]}"; rc=$?; rcl "${NAMES[$i]}" "$rc"; done; PIDS=(); NAMES=(); }
bash -n "$0" || { rcl rest_precheck_bash 1; exit 1; }
for f in 81_fit_inla_inpt.R 81R_inpt_refit_lik.R 85_window_maps.R 88g_sens_suite.R 88h_ppc_varpart.R 88k_struct_kh.R 88k_struct_kh_anem.R 88X_inpt_extra.R; do Rscript -e "invisible(parse('$SC/$f'))" || { rcl "precheck_parse_$f" 1; exit 1; }; done
test -s "$L/results_std3/S7_multiverse_specs.csv" || { rcl rest_precheck_S7specs 1; exit 1; }
rcl rest_precheck 0
export TAG=$T ST_VER=std3 ST_THREADS=1 ST_GRAPH="$GR" SK_RSLOTS=2 NCO=PERIO
cd "$SC"
VERS="s577 s2 sc s8 sexE fix15 x2021 supply"
for V in $VERS; do (cd "$SK" && ST_VER=$V python3 01_스크립트/80_build_disease_panel_v.py $T) > "$LG/80_$V.log" 2>&1; rcl "80_$V" $?; done
for R in 1-50 51-100 101-150 151-200; do wl; ST_MV_RANGE=$R $RS Rscript 81_fit_inla_inpt.R $T > "$LG/81_mv_$R.log" 2>&1 & bgadd $! "81_mv_$R"; sleep 3; done
MVP=("${PIDS[@]}"); MVN=("${NAMES[@]}"); PIDS=(); NAMES=()
wl; $RS Rscript 81R_inpt_refit_lik.R $T > "$LG/81R.log" 2>&1 & bgadd $! 81R_refit_lik
for V in $VERS; do wl; ( if [ $V = s577 ]; then unset ST_GRAPH; fi; ST_VER=$V ST_MV_N=2 ST_SENS_LITE=1 $RS Rscript 81_fit_inla_inpt.R $T ) > "$LG/81_$V.log" 2>&1 & bgadd $! "81_$V"; sleep 2; done
(cd "$SK" && python3 01_스크립트/87_window_grid.py $T std3) > "$LG/87.log" 2>&1; rcl 87_window_grid $?
wl; $RS Rscript 85_window_maps.R $T std3 > "$LG/85.log" 2>&1 & bgadd $! 85_window_maps
(cd "$SK" && python3 01_스크립트/86_yearly_maps_ehsa.py $T std3) > "$LG/86.log" 2>&1; rcl 86_yearly_maps_ehsa $?
(cd "$SK" && ST_VER=std3 python3 01_스크립트/86b_ehsa_modified_mk_inpt.py) > "$LG/86b.log" 2>&1; rcl 86b_ehsa_mk $?
(cd "$SK" && python3 01_스크립트/88g_sens_prep.py $T std3 PERIO) > "$LG/88gprep.log" 2>&1; rcl 88g_sens_prep $?
for X in A B C; do wl; EXG=$X $RS Rscript 88X_inpt_extra.R > "$LG/88X_$X.log" 2>&1 & bgadd $! "88X_$X"; sleep 2; done
for G in nb struct kh; do wl; $RS Rscript 88k_struct_kh.R $G > "$LG/88k_$G.log" 2>&1 & bgadd $! "88k_$G"; sleep 2; done
for i in "${!MVP[@]}"; do wait "${MVP[$i]}"; rcl "${MVN[$i]}" $?; done
ST_MV_PARTS=4 $RS Rscript 81_fit_inla_inpt.R $T > "$LG/81_merge.log" 2>&1; rc=$?; rcl 81_S7_merge $rc
test -s "$L/results_std3/S7_multiverse_summary.csv" || rcl S7_summary_missing 1
for G in 1 2 3 4 5 6; do wl; GROUP=$G $RS Rscript 88g_sens_suite.R > "$LG/88g_$G.log" 2>&1 & bgadd $! "88g_$G"; sleep 2; done
wl; $RS Rscript 88h_ppc_varpart.R > "$LG/88h.log" 2>&1 & bgadd $! 88h
wl; $RS Rscript 88k_struct_kh.R review > "$LG/88k_review.log" 2>&1 & bgadd $! 88k_review
waitall
if grep -q 'KH_typeII FAILED\|KH_typeIV FAILED' "$LG/88k_kh.log"; then KH_INIT=0 $RS Rscript 88k_struct_kh_anem.R kh > "$LG/88k_kh_init0.log" 2>&1; rcl 88k_kh_init0 $?; fi
(cd "$SK" && TAG=$T ST_VER=std3 python3 01_스크립트/89g_sens_combine.py $T std3) > "$LG/89g.log" 2>&1; rcl 89g_sens_combine $?
echo ALL_DONE
