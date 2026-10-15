#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
mkdir -p results reports
${CC:-cc} -O3 analysis/daily_summary.c -lm -o results/daily_summary
results/daily_summary \
  03_Daily_Stocks/data/daily_stock.csv \
  results/daily_last5.csv results/daily_beta_sums.csv results/daily_market.csv
MPLCONFIGDIR=results/.mplcache python3 analysis/run.py

