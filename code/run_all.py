"""Run the complete replication and package its public deliverables."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
TABLES = ROOT / "outputs/tables"
FIGURES = ROOT / "outputs/figures"


def value(path: str, **match: str) -> dict[str, str]:
    with (RESULTS / path).open(newline="") as stream:
        for row in csv.DictReader(stream):
            if all(row[key] == expected for key, expected in match.items()):
                return row
    raise ValueError(f"No matching row in {path}: {match}")


def write_metrics() -> None:
    t1 = value("table1.csv", panel="A_immediate", J="6", K="6", portfolio="winner_minus_loser")
    t2_l = value("table2.csv", portfolio="P1")
    t2_w = value("table2.csv", portfolio="P10")
    t3 = value("table3.csv", subsample="All", portfolio="P10-P1")
    mom = json.loads((RESULTS / "extension_mom_comparison.json").read_text())
    rev = json.loads((RESULTS / "extension_reversal_summary.json").read_text())
    state = json.loads((RESULTS / "extension_market_state_summary.json").read_text())
    event12 = value("event_time_6m_path.csv", event_month="12")
    event36 = value("event_time_6m_path.csv", event_month="36")
    macros = {
        "MainSpread": f'{float(t1["mean_monthly_pct"]):.3f}',
        "MainT": f'{float(t1["t_hac"]):.2f}',
        "LoserBeta": f'{float(t2_l["post_rank_beta"]):.3f}',
        "WinnerBeta": f'{float(t2_w["post_rank_beta"]):.3f}',
        "BetaSpread": f'{float(t2_w["post_rank_beta"])-float(t2_l["post_rank_beta"]):.3f}',
        "MainAlpha": f'{float(t3["panel_b_alpha_pct"]):.3f}',
        "AlphaT": f'{float(t3["panel_b_t_hac"]):.2f}',
        "MomCorrelation": f'{mom["correlation"]:.3f}',
        "MomTrackingError": f'{mom["tracking_error_monthly_pct"]:.3f}',
        "ShortReversal": f'{rev["short_mean_monthly_pct"]:.3f}',
        "LongReversal": f'{rev["long_car_36_pct"]:.2f}',
        "EventTwelve": f'{float(event12["avg_cumulative_spread_pct"]):.2f}',
        "EventThirtySix": f'{float(event36["avg_cumulative_spread_pct"]):.2f}',
        "StateUp": f'{state["states"][0]["mean_spread_pct"]:.3f}',
        "StateDown": f'{state["states"][1]["mean_spread_pct"]:.3f}',
        "StateDifference": f'{state["down_minus_nonnegative_mean_pct"]:.3f}',
        "StateP": f'{state["difference_two_sided_normal_p"]:.3f}',
    }
    (ROOT / "outputs/metrics.tex").write_text(
        "".join(f"\\newcommand{{\\{key}}}{{{number}}}\n" for key, number in macros.items())
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse-daily-summaries", action="store_true",
                        help="Skip the 26 GB daily scan when its three generated summaries already exist")
    args = parser.parse_args()
    if args.reuse_daily_summaries:
        for name in ("daily_last5.csv", "daily_beta_sums.csv", "daily_market.csv"):
            if not (RESULTS / name).is_file():
                parser.error(f"Missing results/{name}; run without --reuse-daily-summaries first")
        env = os.environ.copy()
        env["MPLCONFIGDIR"] = str(RESULTS / ".mplcache")
        subprocess.run(["python3", "analysis/run.py"], cwd=ROOT, env=env, check=True)
    else:
        subprocess.run(["bash", "analysis/run_all.sh"], cwd=ROOT, check=True)

    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    for path in RESULTS.glob("*"):
        if path.suffix in {".csv", ".json"} and path.name not in {
            "daily_last5.csv", "daily_beta_sums.csv", "daily_market.csv"
        }:
            shutil.copy2(path, TABLES / path.name)
        elif path.suffix == ".png":
            shutil.copy2(path, FIGURES / path.name)
    processed = ROOT / "data/processed"
    processed.mkdir(parents=True, exist_ok=True)
    shutil.copy2(RESULTS / "extension_market_state_monthly.csv", processed / "momentum_6_6_monthly.csv")
    write_metrics()
    if not shutil.which("pdflatex"):
        raise SystemExit("pdflatex is required to generate report.pdf. Install TeX Live and rerun.")
    for _ in range(2):
        subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                        "-output-directory", "outputs", "report.tex"], cwd=ROOT,
                       check=True, stdout=subprocess.DEVNULL)
    shutil.copy2(ROOT / "outputs/report.pdf", ROOT / "report.pdf")
    print("Generated report.pdf, outputs/tables/, outputs/figures/, and data/processed/.")


if __name__ == "__main__":
    main()
