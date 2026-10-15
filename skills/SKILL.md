---
name: momentum-replication
description: Reproduce and audit the Jegadeesh-Titman Tables I-III momentum study and extensions from the ECO6067 course CRSP inputs.
---

# Momentum replication

## Inputs

Read `data/raw/README.md` for the three licensed course CSV paths. Keep those CSVs local. Use `reference/F-F_Momentum_Factor.csv` for the public benchmark. Install `requirements.txt`, a C compiler, and TeX Live (`pdflatex`).

## Run

From the repository root, run `python3 code/run_all.py`. The program scans the daily file, computes Tables I–III and extensions, copies aggregate outputs to `outputs/tables/` and `outputs/figures/`, and builds `report.pdf`. After a successful full run, `python3 code/run_all.py --reuse-daily-summaries` can regenerate the Python results and PDF without scanning daily data again.

## Check

Read `results/validation.json` for sample boundaries and counts. Check `outputs/tables/table1.csv` for both immediate and five-trading-day-skip panels, `table2.csv` for 6/6 post-ranking beta and size, and `table3.csv` for raw returns and CAPM alphas. Compare the six-cell constructed Mom factor with the public benchmark using `outputs/tables/extension_mom_comparison.json`. Read `reports/02_replication.md` through `reports/06_validation_and_ai.md` for methods, inference, limitations, and AI disclosure.
