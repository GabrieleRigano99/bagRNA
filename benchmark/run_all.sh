#!/usr/bin/env bash
# Run the benchmark species one after another (shortest first). Waits for each species' inputs
# (SHA256SUMS is written by fetch_inputs.sh when all reads are in). A failed species does not stop the queue.
# Usage: benchmark/run_all.sh [species ...]
cd "$(dirname "$0")/.."
INPUTS=${INPUTS:-/mnt/newvolume/bagRNA_repro/inputs}; LOGS=${LOGS:-/mnt/newvolume/bagRNA_repro/logs}
list=${*:-Candida_orthopsilosis Cryptococcus_depauperatus Aspergillus_aculeatinus Sporothrix_schenckii Syncephalastrum_racemosum}
for sp in $list; do
  until [ -s "$INPUTS/$sp/SHA256SUMS" ]; do sleep 120; done
  date +%s > "$LOGS/$sp.start"
  ./benchmark/run_species.sh "$sp" > "$LOGS/run_$sp.log" 2>&1; echo "EXIT=$?" >> "$LOGS/run_$sp.log"
  date +%s > "$LOGS/$sp.end"
done
