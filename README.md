# YufeiLIU_3036759762 — Jegadeesh–Titman momentum

**Yufei Liu · Student ID 3036759762 · ECO6067 Individual Project 04**

This project reproduces Tables I–III of Jegadeesh and Titman (1993) on the supplied 2000–2025 CRSP CIZ data. It separately studies (1) a constructed momentum factor versus Kenneth French's public Mom, (2) short- and long-horizon reversals, and (3) how the paper-style 6/6 strategy performs across prior market states.

## Repository map

| Path | Purpose |
| --- | --- |
| [report.pdf](report.pdf) | Final English research report. |
| [AI_USE_DISCLOSURE.md](AI_USE_DISCLOSURE.md) | Standalone disclosure by stage of work. |
| [data/raw/README.md](data/raw/README.md) | Manifest for restricted course inputs. |
| [data/processed/](data/processed/) | Derived aggregate monthly data. |
| [code/run_all.py](code/run_all.py) | Single command to regenerate analyses and deliverables. |
| [outputs/tables/](outputs/tables/) | Generated CSV tables and JSON summaries. |
| [outputs/figures/](outputs/figures/) | Generated figures. |
| [skills/SKILL.md](skills/SKILL.md) | Reusable reproduction instructions. |
| [reports/](reports/) | Detailed paper review, replication, extensions, and validation. |

## Inputs

Put the course's numbered folders in this directory without renaming the CSVs:

```
02_Monthly_Stocks_and_Factors/data/monthly_stock.csv
02_Monthly_Stocks_and_Factors/data/F-F factors and RF.csv
03_Daily_Stocks/data/daily_stock.csv
reference/F-F_Momentum_Factor.csv
```

The first three files are restricted course data and are intentionally ignored by Git. The final file is the public US monthly Mom factor, downloaded from the [Kenneth French Data Library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html) on 12 October 2026. Its source file says it uses the 202608 CRSP database. The project only compares dates through 2025-12.

## Reproduce

Requirements: Python 3 with `numpy`, `pandas`, and `matplotlib`, a C compiler (`clang` or `gcc`), and TeX Live (`pdflatex`). No WRDS account is required for the supplied inputs.

Install Python dependencies with `python3 -m pip install -r requirements.txt` if needed.

```sh
python3 code/run_all.py
```

The C program streams the 26 GB daily file twice and creates compact summaries for five-trading-day shifted rankings and prior-year daily betas. The Python program computes portfolios under `results/` and generates detailed result reports under `reports/`. The wrapper copies public aggregate files to `outputs/`, saves the processed monthly 6/6 series, and compiles `report.pdf`. The English paper review is a separate static document. Runtime depends on disk speed; daily extraction may take several minutes. After a complete run, `python3 code/run_all.py --reuse-daily-summaries` refreshes the Python analyses and PDF from existing local daily summaries. Generated security-level summaries remain local and are excluded from Git.

## Outputs

Read the [paper review](reports/01_paper_review.md) first, then the [Tables I–III replication](reports/02_replication.md). The CSV files behind the assigned tables are:

| Paper table | Result file | Contents |
| --- | --- | --- |
| Table I | [results/table1.csv](results/table1.csv) | All 32 J/K specifications: Panel A (immediate formation), Panel B (five-trading-observation gap), decile returns, P10−P1, and HAC t-statistics. |
| Table II | [results/table2.csv](results/table2.csv) | Post-ranking market beta, average formation market capitalization, and stock count for each 6/6 decile. Market capitalization remains in the supplied `MthCap` units. |
| Table III | [results/table3.csv](results/table3.csv) | Both Panel A raw returns and Panel B CAPM alphas, with HAC t-statistics, for the full sample and size/prior-beta groups. |
| Table III joint tests | [results/table3_joint_tests.csv](results/table3_joint_tests.csv) | Test statistics for equal raw decile means and zero CAPM alphas within each group. |

The extensions and supporting outputs are:

| Analysis | Report | Main result files |
| --- | --- | --- |
| Constructed Mom versus French Mom | [reports/03_extension_mom.md](reports/03_extension_mom.md) | [Monthly factor returns](results/extension_mom_monthly.csv), [comparison statistics](results/extension_mom_comparison.json), [cumulative figure](results/extension_mom_cumulative.png) |
| Short- and long-term reversal | [reports/04_extension_reversal.md](reports/04_extension_reversal.md) | [Short-term monthly returns](results/extension_short_reversal.csv), [long-term event path](results/extension_long_reversal_path.csv), [reversal summary](results/extension_reversal_summary.json) |
| 36-month follow-up of 6-month-ranked stocks | [reports/04_extension_reversal.md](reports/04_extension_reversal.md) | [Event-time path and HAC intervals](results/event_time_6m_path.csv), [cohort observations](results/event_time_6m_cohorts.csv), [figure](results/event_time_6m_path.png) |
| 6/6 momentum across prior market states | [reports/05_extension_market_state.md](reports/05_extension_market_state.md) | [Monthly 6/6 returns and market states](results/extension_market_state_monthly.csv), [conditional results](results/extension_market_state_summary.json), [figure](results/extension_market_state.png) |

Finish with [validation and AI-use disclosure](reports/06_validation_and_ai.md). The event-time follow-up is a separate analysis from the 36/36 long-term reversal strategy.

Method choices and deviations from the original papers are documented in the separate reports. The supplied course data must not be uploaded to a public GitHub repository.

