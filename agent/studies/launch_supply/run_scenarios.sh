#!/bin/bash
# Single-seed forced-outcome paths (turn 6: no ensembles). Usage: run_scenarios.sh ARM WORLD SEED [...triples]
# Writes runs/<arm>__<world>__s<seed>.json; skips ones already present. Foreground, one at a time.
cd "$(dirname "$0")/../../.."
D=agent/studies/launch_supply
TR=cash,k_sold,t_sold,k_lost,t_lost,k_age,t_age,k_stock,t_stock,k_short,k_short_any,t_short,t_file_delay,own_used,cmo_used,cmo2_used,l2_used
while [ $# -ge 3 ]; do
  a=$1; w=$2; s=$3; shift 3
  out=$D/runs/${a}__${w}__s${s}.json
  [ -s $out ] && { echo "$a $w s$s cached"; continue; }
  julia --project=. agent/scripts/rdx.jl run $D/scenarios/${a}__${w}.rdj.json --seed=$s --trace=$TR > $out 2>/dev/null
  echo "$a $w s$s exit $?"
done
