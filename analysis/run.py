"""ECO6067: J&T (1993) modern-sample replication and three extensions.

Run from repository root with: python3 analysis/run.py
Raw course CSV files are read only. All generated files go to results/ and reports/.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MONTHLY = ROOT / "02_Monthly_Stocks_and_Factors/data/monthly_stock.csv"
FACTORS = ROOT / "02_Monthly_Stocks_and_Factors/data/F-F factors and RF.csv"
MOM_FILE = ROOT / "reference/F-F_Momentum_Factor.csv"
RESULTS = ROOT / "results"
REPORTS = ROOT / "reports"
J_SET = (3, 6, 9, 12)
K_SET = (3, 6, 9, 12)


def parse_monthly_factor(path: Path, name: str) -> pd.Series:
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        match = re.match(r"^\s*(\d{6})\s*,", line)
        if not match:
            continue
        cells = [c.strip() for c in line.split(",")]
        try:
            value = float(cells[-1] if name == "RF" else cells[1]) / 100.0
        except (ValueError, IndexError):
            continue
        if value <= -9:
            continue
        rows.append((int(match.group(1)), value))
    return pd.Series(dict(rows), name=name).sort_index()


def hac_mean(series: np.ndarray, lags: int) -> tuple[float, float, int]:
    x = np.asarray(series, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 3:
        return np.nan, np.nan, n
    mu = x.mean()
    e = x - mu
    s = e @ e / n
    for k in range(1, min(lags, n - 1) + 1):
        s += 2 * (1 - k / (lags + 1)) * (e[k:] @ e[:-k]) / n
    se = np.sqrt(max(s, 0) / n)
    return mu, mu / se if se > 0 else np.nan, n


def hac_alpha(y: np.ndarray, market: np.ndarray, rf: np.ndarray, lags: int,
              spread: bool = False) -> tuple[float, float, int]:
    mask = np.isfinite(y) & np.isfinite(market) & np.isfinite(rf)
    yy, mm, rr = y[mask], market[mask], rf[mask]
    n = len(yy)
    if n < 12:
        return np.nan, np.nan, n
    target = yy if spread else yy - rr
    X = np.column_stack([np.ones(n), mm - rr])
    coef = np.linalg.lstsq(X, target, rcond=None)[0]
    e = target - X @ coef
    z = X * e[:, None]
    meat = z.T @ z
    for k in range(1, min(lags, n - 1) + 1):
        cross = z[k:].T @ z[:-k]
        meat += (1 - k / (lags + 1)) * (cross + cross.T)
    inv = np.linalg.inv(X.T @ X)
    se = np.sqrt(max((inv @ meat @ inv)[0, 0], 0))
    return coef[0], coef[0] / se if se > 0 else np.nan, n


def hac_matrix(v: np.ndarray, lags: int) -> np.ndarray:
    """Long-run covariance of rows, without the division by sample size."""
    n = len(v)
    centered = v - v.mean(axis=0)
    out = centered.T @ centered / n
    for k in range(1, min(lags, n-1)+1):
        cross = centered[k:].T @ centered[:-k] / n
        out += (1-k/(lags+1)) * (cross + cross.T)
    return out


def table3_joint_tests(port: np.ndarray, market: np.ndarray, rf: np.ndarray, lags: int = 5):
    valid = np.isfinite(port).all(axis=1) & np.isfinite(market) & np.isfinite(rf)
    y, m, r = port[valid], market[valid], rf[valid]
    n = len(y)
    if n < 20:
        return np.nan, np.nan, n
    contrast = np.column_stack([-np.ones(9), np.eye(9)])
    diff = contrast @ y.mean(axis=0)
    cov_diff = contrast @ hac_matrix(y, lags) @ contrast.T / n
    raw_wald = float(diff @ np.linalg.pinv(cov_diff) @ diff)
    X = np.column_stack([np.ones(n), m-r])
    inv = np.linalg.inv(X.T @ X)
    coef = inv @ X.T @ (y-r[:, None])
    residual = y-r[:, None]-X@coef
    influence = (X @ inv[:, 0])[:, None] * residual
    alpha_cov = hac_matrix(influence, lags) * n
    alpha_wald = float(coef[0] @ np.linalg.pinv(alpha_cov) @ coef[0])
    return raw_wald, alpha_wald, n


def deciles(score: np.ndarray, eligible: np.ndarray) -> np.ndarray:
    out = np.full(len(score), -1, dtype=np.int8)
    idx = np.flatnonzero(eligible & np.isfinite(score))
    if len(idx) < 100:
        return out
    order = idx[np.argsort(score[idx], kind="mergesort")]
    for d, part in enumerate(np.array_split(order, 10)):
        out[part] = d
    return out


def weighted_group_return(r: np.ndarray, groups: np.ndarray,
                          weight: np.ndarray | None = None, count: int = 10) -> np.ndarray:
    valid = (groups >= 0) & np.isfinite(r)
    if weight is None:
        sums = np.bincount(groups[valid], weights=r[valid], minlength=count)
        n = np.bincount(groups[valid], minlength=count)
    else:
        valid &= np.isfinite(weight) & (weight > 0)
        sums = np.bincount(groups[valid], weights=r[valid] * weight[valid], minlength=count)
        n = np.bincount(groups[valid], weights=weight[valid], minlength=count)
    return np.divide(sums, n, out=np.full(count, np.nan), where=n > 0)


def aggregate_vintages(vintage: np.ndarray, k: int) -> np.ndarray:
    T, _, D = vintage.shape
    out = np.full((T, D), np.nan)
    for m in range(k - 1, T):
        legs = np.stack([vintage[m - a, a] for a in range(k)])
        valid = np.all(np.isfinite(legs), axis=0)
        out[m, valid] = legs[:, valid].mean(axis=0)
    return out


def make_vintage(ret: np.ndarray, score: np.ndarray, eligible: np.ndarray,
                 cap: np.ndarray, max_k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    T, _ = ret.shape
    vintage = np.full((T, max_k, 10), np.nan, dtype=np.float32)
    cap_means = np.full((T, 10), np.nan)
    counts = np.zeros((T, 10), dtype=int)
    for t in range(1, T):
        rank = deciles(score[t], eligible[t - 1])
        if np.all(rank < 0):
            continue
        for d in range(10):
            selected = rank == d
            counts[t, d] = selected.sum()
            cap_means[t, d] = np.nanmean(cap[t - 1, selected])
        for age in range(min(max_k, T - t)):
            vintage[t, age] = weighted_group_return(ret[t + age], rank)
    return vintage, cap_means, counts


def load_data():
    cols = ["PERMNO", "YYYYMM", "MthCalDt", "SecInfoStartDt", "SecInfoEndDt",
            "MthRet", "MthCap", "vwretd", "PrimaryExch",
            "SecurityType", "SecuritySubType", "ShareType", "USIncFlg", "IssuerType"]
    print("Reading monthly stock data...", flush=True)
    df = pd.read_csv(MONTHLY, usecols=cols, low_memory=False)
    dup_count = int(df.duplicated(["PERMNO", "YYYYMM"]).sum())
    if dup_count:
        conflicting_returns = df.groupby(["PERMNO", "YYYYMM"]).MthRet.nunique(dropna=False).gt(1).sum()
        if conflicting_returns:
            raise ValueError(f"{conflicting_returns} duplicate security-month keys have conflicting returns")
        df["valid_info"] = (df.MthCalDt >= df.SecInfoStartDt) & (df.MthCalDt <= df.SecInfoEndDt)
        df["common_priority"] = (df.SecurityType.eq("EQTY") & df.SecuritySubType.eq("COM") &
                                  df.ShareType.eq("NS") & df.USIncFlg.eq("Y") & df.IssuerType.eq("CORP"))
        df = df.sort_values(["PERMNO", "YYYYMM", "valid_info", "common_priority"],
                            ascending=[True, True, False, False]).drop_duplicates(["PERMNO", "YYYYMM"])
        print(f"Removed {dup_count} duplicate security-month event rows", flush=True)
    months = np.array([int(p.strftime("%Y%m")) for p in pd.period_range("2000-01", "2025-12", freq="M")])
    permnos = np.sort(df.PERMNO.unique())
    def pivot(column: str, dtype=float):
        return df.pivot(index="YYYYMM", columns="PERMNO", values=column).reindex(index=months, columns=permnos).to_numpy(dtype=dtype)
    ret = pivot("MthRet")
    cap = pivot("MthCap")
    common = (df.SecurityType.eq("EQTY") & df.SecuritySubType.eq("COM") &
              df.ShareType.eq("NS") & df.USIncFlg.eq("Y") & df.IssuerType.eq("CORP"))
    df["eligible_jt"] = common & df.PrimaryExch.isin(["N", "A"])
    df["eligible_mom"] = common & df.PrimaryExch.isin(["N", "A", "Q"])
    df["nyse"] = df.PrimaryExch.eq("N")
    jt = pivot("eligible_jt", dtype=float) == 1
    mom = pivot("eligible_mom", dtype=float) == 1
    nyse = pivot("nyse", dtype=float) == 1
    market = df.groupby("YYYYMM").vwretd.first().reindex(months).to_numpy(dtype=float)
    del df
    tail_path = RESULTS / "daily_last5.csv"
    if not tail_path.exists():
        raise FileNotFoundError("Run analysis/daily_summary.c before analysis/run.py")
    tail = pd.read_csv(tail_path)
    tail_gross = tail.pivot(index="YYYYMM", columns="PERMNO", values="last5_gross").reindex(index=months, columns=permnos).to_numpy(dtype=float)
    rf_series = parse_monthly_factor(FACTORS, "RF")
    rf = rf_series.reindex(months).to_numpy(dtype=float)
    official_series = parse_monthly_factor(MOM_FILE, "Mom")
    official = official_series.reindex(months).to_numpy(dtype=float)
    return months, permnos, ret, cap, jt, mom, nyse, market, rf, official, tail_gross


def signals(ret: np.ndarray, tail: np.ndarray) -> tuple[dict[int, np.ndarray], dict[int, np.ndarray]]:
    T, N = ret.shape
    logret = np.full_like(ret, np.nan)
    valid = np.isfinite(ret) & (ret > -1)
    logret[valid] = np.log1p(ret[valid])
    cum = np.vstack([np.zeros((1, N)), np.cumsum(np.nan_to_num(logret, nan=0), axis=0)])
    cnt = np.vstack([np.zeros((1, N), dtype=np.int16), np.cumsum(np.isfinite(logret), axis=0, dtype=np.int16)])
    signal_a, signal_b = {}, {}
    for j in J_SET:
        a = np.full_like(ret, np.nan)
        b = np.full_like(ret, np.nan)
        for t in range(j, T):
            full = cnt[t] - cnt[t-j] == j
            a[t, full] = np.expm1(cum[t, full] - cum[t-j, full])
            if t > j and np.any(full):
                l0, l1 = tail[t-j-1], tail[t-1]
                good = full & np.isfinite(l0) & np.isfinite(l1) & (l0 > 0) & (l1 > 0)
                b[t, good] = np.expm1(cum[t, good] - cum[t-j, good] + np.log(l0[good]) - np.log(l1[good]))
        signal_a[j], signal_b[j] = a, b
    return signal_a, signal_b


def make_table1(months, ret, cap, jt, sig_a, sig_b):
    rows = []
    series = {}
    for panel, signals_by_j in [("A_immediate", sig_a), ("B_skip_5_trading_days", sig_b)]:
        for j in J_SET:
            print(f"Table I: {panel}, J={j}", flush=True)
            vintage, cap_means, counts = make_vintage(ret, signals_by_j[j], jt, cap, 12)
            if panel == "A_immediate" and j == 6:
                base = (vintage, cap_means, counts)
            for k in K_SET:
                port = aggregate_vintages(vintage, k)
                port[:24] = np.nan  # common January 2002 start, full overlap for all J/K
                series[(panel, j, k)] = port
                for label, col in [("P1_loser", port[:, 0]), ("P10_winner", port[:, 9]), ("winner_minus_loser", port[:, 9] - port[:, 0])]:
                    mean, tstat, n = hac_mean(col, k - 1)
                    rows.append(dict(panel=panel, J=j, K=k, portfolio=label,
                                     mean_monthly_pct=mean * 100, t_hac=tstat, months=n))
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "table1.csv", index=False)
    return out, series, base


def make_table2(months, base, market):
    vintage, cap_means, counts = base
    port = aggregate_vintages(vintage, 6)
    port[:24] = np.nan
    rows = []
    for d in range(10):
        mask = np.isfinite(port[:, d]) & np.isfinite(market)
        beta = np.cov(port[mask, d], market[mask], ddof=1)[0, 1] / np.var(market[mask], ddof=1)
        rows.append(dict(portfolio=f"P{d+1}", post_rank_beta=beta,
                         avg_formation_market_cap=np.nanmean(cap_means[24:, d]),
                         avg_stocks_per_vintage=np.nanmean(np.where(counts[24:, d] > 0, counts[24:, d], np.nan)),
                         months=mask.sum()))
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "table2.csv", index=False)
    return out


def make_market_state(months, port, market):
    """Compare Table I's immediate 6/6 spread across predetermined market states."""
    n_months = len(months)
    if port.shape != (n_months, 10):
        raise ValueError("The market-state extension requires ten monthly 6/6 portfolios")
    prior_12m = np.full(n_months, np.nan)
    for t in range(12, n_months):
        history = market[t-12:t]
        if np.isfinite(history).all() and (history > -1).all():
            prior_12m[t] = np.prod(1 + history) - 1
    spread = port[:, 9] - port[:, 0]
    valid = np.isfinite(spread) & np.isfinite(prior_12m) & np.isfinite(market)
    if valid.sum() < 100:
        raise ValueError("Too few monthly 6/6 observations for the market-state extension")
    down = prior_12m[valid] < 0
    if min(down.sum(), (~down).sum()) < 20:
        raise ValueError("Too few observations in a market-state group")
    y = spread[valid]
    x = down.astype(float)
    design = np.column_stack([np.ones(len(y)), x])
    coefficients = np.linalg.lstsq(design, y, rcond=None)[0]
    residual = y - design @ coefficients
    scores = design * residual[:, None]
    meat = scores.T @ scores
    lags = 5  # the Table I 6/6 holding portfolios overlap for six months
    for lag in range(1, lags + 1):
        cross = scores[lag:].T @ scores[:-lag]
        meat += (1 - lag / (lags + 1)) * (cross + cross.T)
    inv = np.linalg.inv(design.T @ design)
    covariance = inv @ meat @ inv
    state_rows = []
    for label, mask, contrast in [
        ("prior_12m_nonnegative", ~down, np.array([1.0, 0.0])),
        ("prior_12m_negative", down, np.array([1.0, 1.0])),
    ]:
        mean = float(contrast @ coefficients)
        se = math.sqrt(max(float(contrast @ covariance @ contrast), 0))
        values = y[mask]
        state_rows.append(dict(state=label, months=int(mask.sum()),
                               mean_spread_pct=mean*100, hac_se_pct=se*100,
                               ci95_low_pct=(mean-1.96*se)*100,
                               ci95_high_pct=(mean+1.96*se)*100,
                               negative_spread_months=int((values < 0).sum()),
                               negative_spread_share=float((values < 0).mean()),
                               p10_spread_pct=float(np.quantile(values, .1)*100),
                               worst_spread_pct=float(values.min()*100)))
    difference = float(coefficients[1])
    difference_se = math.sqrt(max(float(covariance[1, 1]), 0))
    difference_t = difference / difference_se if difference_se > 0 else np.nan
    rebound = down & (market[valid] > 0)
    down_other = down & ~(market[valid] > 0)
    monthly = pd.DataFrame(dict(month=months[valid],
                                p1_loser_pct=port[valid, 0]*100,
                                p10_winner_pct=port[valid, 9]*100,
                                winner_minus_loser_pct=y*100,
                                market_return_pct=market[valid]*100,
                                prior_12m_market_return_pct=prior_12m[valid]*100,
                                prior_market_down=down,
                                same_month_market_rebound=rebound))
    monthly.to_csv(RESULTS / "extension_market_state_monthly.csv", index=False)
    worst = monthly.nsmallest(1, "winner_minus_loser_pct").iloc[0]
    without_worst = monthly[monthly.month != worst.month]
    state_means_without_worst = without_worst.groupby("prior_market_down").winner_minus_loser_pct.mean()
    summary = dict(sample_start=int(months[valid][0]), sample_end=int(months[valid][-1]),
                   months=int(valid.sum()), state_definition="Prior 12 complete monthly market returns through t-1 compound to a negative return",
                   hac_lags=lags, states=state_rows,
                   down_minus_nonnegative_mean_pct=difference*100,
                   difference_hac_se_pct=difference_se*100,
                   difference_hac_t=float(difference_t),
                   difference_two_sided_normal_p=float(math.erfc(abs(difference_t)/math.sqrt(2))),
                   difference_ci95_low_pct=(difference-1.96*difference_se)*100,
                   difference_ci95_high_pct=(difference+1.96*difference_se)*100,
                   prior_down_current_positive_months=int(rebound.sum()),
                   prior_down_current_positive_mean_pct=float(y[rebound].mean()*100) if rebound.any() else None,
                   prior_down_current_nonpositive_months=int(down_other.sum()),
                   prior_down_current_nonpositive_mean_pct=float(y[down_other].mean()*100) if down_other.any() else None,
                   worst_month=int(worst.month),
                   worst_month_spread_pct=float(worst.winner_minus_loser_pct),
                   down_minus_nonnegative_without_worst_pct=float(state_means_without_worst[True] - state_means_without_worst[False]),
                   five_worst_months=monthly.nsmallest(5, "winner_minus_loser_pct")[
                       ["month", "winner_minus_loser_pct", "market_return_pct", "prior_market_down"]
                   ].to_dict(orient="records"))
    (RESULTS / "extension_market_state_summary.json").write_text(json.dumps(summary, indent=2))
    fig, ax = plt.subplots(figsize=(7, 4.5))
    means = np.array([row["mean_spread_pct"] for row in state_rows])
    errors = np.array([[row["mean_spread_pct"]-row["ci95_low_pct"] for row in state_rows],
                       [row["ci95_high_pct"]-row["mean_spread_pct"] for row in state_rows]])
    ax.bar(["Prior 12m nonnegative", "Prior 12m negative"], means,
           color=["#47769a", "#c27248"])
    ax.errorbar([0, 1], means, yerr=errors, fmt="none", ecolor="black", capsize=5)
    ax.axhline(0, color="black", lw=.8)
    ax.set(title="6/6 momentum by prior market state", ylabel="Mean monthly P10−P1, percent")
    fig.tight_layout()
    fig.savefig(RESULTS / "extension_market_state.png", dpi=160)
    plt.close(fig)
    return summary


def load_prior_betas(months, permnos):
    print("Reading daily beta summaries...", flush=True)
    sums = pd.read_csv(RESULTS / "daily_beta_sums.csv")
    mkt = pd.read_csv(RESULTS / "daily_market.csv")
    mkt["year"] = mkt.date_key // 10000
    rho = {}
    for year, group in mkt.groupby("year"):
        x = group.vwretd.to_numpy()
        rho[year] = np.corrcoef(x[1:], x[:-1])[0, 1]
    slopes = []
    for j in range(3):
        n = sums[f"n{j}"].to_numpy()
        sx = sums[f"sx{j}"].to_numpy()
        sy = sums[f"sy{j}"].to_numpy()
        sxx = sums[f"sxx{j}"].to_numpy()
        sxy = sums[f"sxy{j}"].to_numpy()
        denom = sxx - sx*sx / np.maximum(n, 1)
        slopes.append(np.divide(sxy - sx*sy / np.maximum(n, 1), denom,
                                out=np.full(len(sums), np.nan), where=(n >= 100) & (denom > 0)))
    sums["beta_sw"] = np.nansum(slopes, axis=0) / (1 + 2 * sums.year.map(rho).to_numpy())
    sums.loc[~np.isfinite(np.stack(slopes)).all(axis=0), "beta_sw"] = np.nan
    sums.loc[~sums.beta_sw.between(-10, 10), "beta_sw"] = np.nan
    by_year = sums.pivot(index="year", columns="PERMNO", values="beta_sw").reindex(columns=permnos)
    out = np.full((len(months), len(permnos)), np.nan, dtype=np.float32)
    for t, month in enumerate(months):
        prior = month // 100 - 1
        if prior in by_year.index:
            out[t] = by_year.loc[prior].to_numpy(dtype=np.float32)
    return out, int(sums.beta_sw.notna().sum())


def make_table3(months, ret, cap, jt, sig6, market, rf, beta):
    T, _ = ret.shape
    names = ["All", "Size_small", "Size_mid", "Size_big", "Beta_low", "Beta_mid", "Beta_high"]
    vintages = {name: np.full((T, 6, 10), np.nan, dtype=np.float32) for name in names}
    for t in range(6, T):
        if t % 36 == 0:
            print(f"Table III: forming {months[t]}", flush=True)
        universe = jt[t-1] & np.isfinite(sig6[t]) & np.isfinite(cap[t-1]) & (cap[t-1] > 0)
        size_breaks = np.quantile(cap[t-1, universe], [1/3, 2/3]) if universe.sum() >= 100 else [np.nan, np.nan]
        beta_univ = universe & np.isfinite(beta[t])
        beta_breaks = np.quantile(beta[t, beta_univ], [1/3, 2/3]) if beta_univ.sum() >= 100 else [np.nan, np.nan]
        groups = {
            "All": universe,
            "Size_small": universe & (cap[t-1] <= size_breaks[0]),
            "Size_mid": universe & (cap[t-1] > size_breaks[0]) & (cap[t-1] <= size_breaks[1]),
            "Size_big": universe & (cap[t-1] > size_breaks[1]),
            "Beta_low": beta_univ & (beta[t] <= beta_breaks[0]),
            "Beta_mid": beta_univ & (beta[t] > beta_breaks[0]) & (beta[t] <= beta_breaks[1]),
            "Beta_high": beta_univ & (beta[t] > beta_breaks[1]),
        }
        for name, eligible in groups.items():
            rank = deciles(sig6[t], eligible)
            if np.all(rank < 0):
                continue
            for age in range(min(6, T-t)):
                vintages[name][t, age] = weighted_group_return(ret[t+age], rank)
    rows = []
    joint_rows = []
    for name in names:
        port = aggregate_vintages(vintages[name], 6)
        port[:24] = np.nan
        raw_wald, alpha_wald, joint_n = table3_joint_tests(port, market, rf)
        joint_rows.append(dict(subsample=name, panel_a_equal_means_chi2_9=raw_wald,
                               panel_b_zero_alphas_chi2_10=alpha_wald, months=joint_n))
        for d in range(11):
            y = port[:, d] if d < 10 else port[:, 9] - port[:, 0]
            mean, tstat, n = hac_mean(y, 5)
            alpha, t_alpha, _ = hac_alpha(y, market, rf, 5, spread=(d == 10))
            rows.append(dict(subsample=name, portfolio=f"P{d+1}" if d < 10 else "P10-P1",
                             panel_a_mean_pct=mean*100, panel_a_t_hac=tstat,
                             panel_b_alpha_pct=alpha*100, panel_b_t_hac=t_alpha, months=n))
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "table3.csv", index=False)
    pd.DataFrame(joint_rows).to_csv(RESULTS / "table3_joint_tests.csv", index=False)
    return out


def make_mom(months, ret, cap, mom_univ, nyse, official):
    T, N = ret.shape
    built = np.full(T, np.nan)
    details = []
    for t in range(13, T):
        prior = ret[t-12:t-1]  # months t-12 through t-2
        valid_hist = np.isfinite(prior) & (prior > -1)
        valid = valid_hist.all(axis=0)
        score = np.full(N, np.nan)
        score[valid] = np.expm1(np.log1p(prior[:, valid]).sum(axis=0))
        eligible = (mom_univ[t-1] & valid & np.isfinite(cap[t-1]) & (cap[t-1] > 0)
                    & np.isfinite(cap[t-13]) & (cap[t-13] > 0))
        reference = eligible & nyse[t-1]
        if reference.sum() < 100:
            continue
        med_cap = np.median(cap[t-1, reference])
        low, high = np.quantile(score[reference], [0.3, 0.7])
        size = np.where(cap[t-1] <= med_cap, 0, 1)
        momentum = np.where(score <= low, 0, np.where(score >= high, 2, 1))
        cell = size * 3 + momentum
        cell[~eligible] = -1
        r6 = weighted_group_return(ret[t], cell, cap[t-1], count=6)
        if np.all(np.isfinite(r6[[0, 2, 3, 5]])):
            built[t] = (r6[2] + r6[5] - r6[0] - r6[3]) / 2
        details.append(dict(month=months[t], eligible_stocks=int(eligible.sum()), nyse_breakpoint_stocks=int(reference.sum()),
                            small_low=r6[0], small_high=r6[2], big_low=r6[3], big_high=r6[5], mom_built=built[t]))
    pd.DataFrame(details).to_csv(RESULTS / "extension_mom_monthly.csv", index=False)
    mask = np.isfinite(built) & np.isfinite(official) & (months >= 200201)
    x, y = built[mask], official[mask]
    comparison = dict(start=int(months[mask][0]), end=int(months[mask][-1]), months=int(mask.sum()),
                      built_mean_pct=float(x.mean()*100), official_mean_pct=float(y.mean()*100),
                      built_vol_monthly_pct=float(x.std(ddof=1)*100), official_vol_monthly_pct=float(y.std(ddof=1)*100),
                      correlation=float(np.corrcoef(x,y)[0,1]), mean_difference_pct=float((x-y).mean()*100),
                      tracking_error_monthly_pct=float((x-y).std(ddof=1)*100))
    (RESULTS / "extension_mom_comparison.json").write_text(json.dumps(comparison, indent=2))
    fig, ax = plt.subplots(figsize=(9, 4.5))
    dates = pd.to_datetime(months[mask].astype(str), format="%Y%m")
    ax.plot(dates, np.cumprod(1+x)-1, label="Constructed Mom")
    ax.plot(dates, np.cumprod(1+y)-1, label="French Mom")
    ax.set(title="Momentum factor cumulative returns, common sample", ylabel="Cumulative return")
    ax.legend(); ax.grid(alpha=.25); fig.tight_layout()
    fig.savefig(RESULTS / "extension_mom_cumulative.png", dpi=160)
    plt.close(fig)
    return comparison


def make_reversal(months, ret, cap, jt, market):
    T, N = ret.shape
    short = np.full(T, np.nan)
    for t in range(1, T):
        eligible = jt[t-1] & np.isfinite(ret[t-1])
        rank = deciles(ret[t-1], eligible)
        r = weighted_group_return(ret[t], rank)
        short[t] = r[0] - r[9]
    short[:24] = np.nan
    mean, tstat, n = hac_mean(short, 1)
    pd.DataFrame({"month": months, "loser_minus_winner": short}).dropna().to_csv(RESULTS / "extension_short_reversal.csv", index=False)
    events = []
    for t, month in enumerate(months):
        if month % 100 != 12 or month < 200212 or t + 36 >= T or (month // 100 - 2002) % 3:
            continue
        window = ret[t-35:t+1]
        eligible = jt[t] & np.isfinite(window).all(axis=0)
        if eligible.sum() < 100:
            continue
        signal = np.nansum(window - market[t-35:t+1, None], axis=0)
        rank = deciles(signal, eligible)
        for age in range(1, 37):
            r = weighted_group_return(ret[t+age], rank)
            events.append(dict(formation_month=month, event_month=age,
                               loser_market_adjusted=r[0]-market[t+age],
                               winner_market_adjusted=r[9]-market[t+age],
                               loser_minus_winner=r[0]-r[9]))
    events = pd.DataFrame(events)
    events.to_csv(RESULTS / "extension_long_reversal_events.csv", index=False)
    event_avg = events.groupby("event_month")[["loser_market_adjusted", "winner_market_adjusted", "loser_minus_winner"]].mean()
    event_avg["cum_loser_minus_winner"] = event_avg.loser_minus_winner.cumsum()
    event_avg.to_csv(RESULTS / "extension_long_reversal_path.csv")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(event_avg.index, event_avg.cum_loser_minus_winner * 100)
    ax.axhline(0, color="black", lw=.7)
    ax.set(title="Long-term reversal: cumulative loser minus winner", xlabel="Months after formation", ylabel="Cumulative percentage points")
    ax.grid(alpha=.25); fig.tight_layout()
    fig.savefig(RESULTS / "extension_long_reversal.png", dpi=160)
    plt.close(fig)
    summary = dict(short_mean_monthly_pct=mean*100, short_t_hac=tstat, short_months=n,
                   long_cohorts=int(events.formation_month.nunique()),
                   long_formations=sorted(events.formation_month.unique().tolist()),
                   long_car_12_pct=float(event_avg.loc[12,"cum_loser_minus_winner"]*100),
                   long_car_24_pct=float(event_avg.loc[24,"cum_loser_minus_winner"]*100),
                   long_car_36_pct=float(event_avg.loc[36,"cum_loser_minus_winner"]*100))
    (RESULTS / "extension_reversal_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def make_event_time6(months, ret, cap, jt, signal6):
    """Follow the same six-month winner and loser cohorts for 36 event months."""
    vintage, _, counts = make_vintage(ret, signal6, jt, cap, 36)
    formation = np.flatnonzero((months >= 200201) & (np.arange(len(months)) + 36 <= len(months)))
    complete = (counts[formation, 0] > 0) & (counts[formation, 9] > 0)
    complete &= np.isfinite(vintage[formation, :, 0]).all(axis=1)
    complete &= np.isfinite(vintage[formation, :, 9]).all(axis=1)
    formation = formation[complete]
    if len(formation) < 100:
        raise ValueError("Too few complete 6-month formation cohorts for the 36-month event study")
    loser = vintage[formation, :, 0].astype(float)
    winner = vintage[formation, :, 9].astype(float)
    spread = winner - loser
    cumulative = np.cumsum(spread, axis=1)
    n = len(formation)
    pd.DataFrame({
        "formation_month": np.repeat(months[formation], 36),
        "event_month": np.tile(np.arange(1, 37), n),
        "p1_loser_return": loser.ravel(),
        "p10_winner_return": winner.ravel(),
        "winner_minus_loser": spread.ravel(),
        "cumulative_winner_minus_loser": cumulative.ravel(),
    }).to_csv(RESULTS / "event_time_6m_cohorts.csv", index=False)
    rows = []
    for h in range(1, 37):
        mean, tstat, observed = hac_mean(cumulative[:, h-1], h-1)
        se = abs(mean / tstat) if np.isfinite(tstat) and tstat != 0 else np.nan
        rows.append(dict(event_month=h, cohorts=observed,
                         avg_loser_monthly_pct=loser[:, h-1].mean()*100,
                         avg_winner_monthly_pct=winner[:, h-1].mean()*100,
                         avg_spread_monthly_pct=spread[:, h-1].mean()*100,
                         avg_cumulative_spread_pct=mean*100,
                         hac_se_cumulative_pct=se*100,
                         hac_t_cumulative=tstat,
                         ci95_low_pct=(mean-1.96*se)*100,
                         ci95_high_pct=(mean+1.96*se)*100))
    path = pd.DataFrame(rows)
    path.to_csv(RESULTS / "event_time_6m_path.csv", index=False)
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    ax.plot(path.event_month, path.avg_cumulative_spread_pct, label="Winner minus loser")
    ax.fill_between(path.event_month.to_numpy(), path.ci95_low_pct.to_numpy(),
                    path.ci95_high_pct.to_numpy(), alpha=.2, label="Approx. 95% HAC interval")
    ax.axhline(0, color="black", lw=.8)
    ax.set(title="6-month ranking: spread after formation", xlabel="Months after formation",
           ylabel="Cumulative percentage points")
    ax.legend(); ax.grid(alpha=.25); fig.tight_layout()
    fig.savefig(RESULTS / "event_time_6m_path.png", dpi=160)
    plt.close(fig)
    return dict(first_formation=int(months[formation[0]]),
                last_formation=int(months[formation[-1]]), cohorts=n,
                horizons={h: path.iloc[h-1].to_dict() for h in (12, 24, 36)})


def write_reports(table1, table2, table3, mom, reversal, event_time, market_state, validation):
    pick = table1[(table1.panel == "A_immediate") & (table1.J == 6) & (table1.K == 6) & (table1.portfolio == "winner_minus_loser")].iloc[0]
    all_spread = table3[(table3.subsample == "All") & (table3.portfolio == "P10-P1")].iloc[0]
    p1, p10 = table2.iloc[0], table2.iloc[-1]
    worst_rows = "\n".join(
        f"| {row['month']} | {row['winner_minus_loser_pct']:.2f}% | {row['market_return_pct']:.2f}% | {'Yes' if row['prior_market_down'] else 'No'} |"
        for row in market_state["five_worst_months"])
    (REPORTS / "02_replication.md").write_text(f"""# Tables I–III: modern-sample replication

## Question and sample

Do past winners outperform past losers over the next 3–12 months? This is a modern-sample reproduction of Jegadeesh and Titman (1993), using supplied CRSP CIZ data from 2000–2025. The common reported return sample starts in January 2002, allowing all J/K specifications a complete formation and overlapping-holding history. The original paper studies NYSE/AMEX stocks in 1965–1989, so numerical equality is not expected.

## Construction

- Universe: US-incorporated corporate common shares (`EQTY/COM/NS`) on NYSE or AMEX (`PrimaryExch` N/A) at formation. Returns with delisting are retained where supplied.
- At the start of month t, rank J-month compounded total returns through t−1 into ten equally sized groups, P1 lowest through P10 highest. An equal-weight vintage is formed each month and held K months. In each month, returns of the K active vintages are averaged. Monthly equal weights are reset within each vintage.
- Table I Panel A uses immediately preceding J months. Panel B shifts the J-month measurement interval five trading observations earlier, using the last-five-daily-return gross factors. This operational definition approximates the paper's one-week interval.
- Table II uses the 6/6 Panel A portfolios. The post-ranking beta is the time-series slope of each portfolio's monthly returns on the CRSP value-weighted market return. Size is average formation-date market capitalization (`MthCap`, source units).
- Table III uses 6/6 Panel A within size terciles and within prior-calendar-year daily Scholes–Williams-beta terciles. Panel A reports raw mean monthly returns. Panel B reports the CAPM intercept from portfolio excess returns regressed on CRSP value-weighted market excess returns. The risk-free rate is converted from percent to decimal.
- Reported t statistics use Newey–West standard errors with K−1 monthly lags, because overlapping holdings induce serial correlation. This is a documented inference choice and can differ from the paper's reported t values.

## Core result

The 6/6 immediate-formation P10−P1 portfolio averages **{pick.mean_monthly_pct:.3f}% per month** (HAC t = **{pick.t_hac:.2f}**, {int(pick.months)} months). The complete 32-strategy grid is in `results/table1.csv`; `results/table2.csv` and `results/table3.csv` contain the other assigned tables.

## Interpretation and limits

The original paper's 12/3 spread is 1.31% per month without a gap and 1.49% with a one-week gap. These are reference values from the 1965–1989 sample, not targets for the 2000–2025 data. CIZ monthly returns compound daily returns and can differ from the historical CRSP format. The five-trading-observation interval is a transparent approximation to a calendar week. Security universe, data revisions, sample period, and inference method also affect comparisons. Read the full tables before interpreting sign or significance of any specification.

Source: https://doi.org/10.1111/j.1540-6261.1993.tb04702.x
""")
    (REPORTS / "03_extension_mom.md").write_text(f"""# Extension 1: constructed Mom versus Kenneth French Mom

## Question

How closely does a momentum factor built from the supplied CRSP sample track the public Kenneth French US Mom factor?

## Method

At each month-end, require common-stock status and positive prior market equity. Use months t−12 through t−2 to measure prior return, omitting t−1. Set size and momentum breakpoints from NYSE securities: median size and 30th/70th prior-return percentiles. Assign NYSE, AMEX, and NASDAQ stocks to the six size-by-momentum groups and value-weight their month-t returns using t−1 market equity. Construct ½(Small High + Big High) − ½(Small Low + Big Low). This follows the public definition but may differ from the official series because of eligibility details, source vintages, missing-return conventions, and CRSP revisions.

## Results, common {mom['start']}–{mom['end']} sample

| Measure | Constructed | French official |
| --- | ---: | ---: |
| Average monthly return | {mom['built_mean_pct']:.3f}% | {mom['official_mean_pct']:.3f}% |
| Monthly volatility | {mom['built_vol_monthly_pct']:.3f}% | {mom['official_vol_monthly_pct']:.3f}% |

There are {mom['months']} common months. Monthly correlation is **{mom['correlation']:.3f}**; constructed minus official mean is **{mom['mean_difference_pct']:.3f} percentage points**, with **{mom['tracking_error_monthly_pct']:.3f} percentage-point** monthly tracking error. See `results/extension_mom_monthly.csv` and `results/extension_mom_cumulative.png`.

The official factor is a value-weighted 2×3 sort with a one-month skip, while the paper's Tables I–III use equal-weighted deciles with J/K overlapping holding periods. Their return levels should not be treated as the same portfolio.

Source: https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_mom_factor.html
""")
    (REPORTS / "04_extension_reversal.md").write_text(f"""# Extension 2: reversal and the subsequent path of momentum

## Question

How does return predictability vary across one-month reversal, three-year reversal, and the later performance of stocks selected by a six-month momentum signal?

## Short horizon

Within the same NYSE/AMEX common-stock universe, rank month t−1 returns into deciles and hold for month t. Define reversal as P1−P10. It averages **{reversal['short_mean_monthly_pct']:.3f}% per month** (HAC t = **{reversal['short_t_hac']:.2f}**, {reversal['short_months']} months). The monthly path is in `results/extension_short_reversal.csv`. This is a transparent one-month test motivated by Jegadeesh (1990), not a literal reproduction of his forecasting regressions.

## Long horizon

Every third December, rank stocks on the sum of market-adjusted returns over the prior 36 months. Form equal-weighted bottom and top deciles and follow their next 36 monthly returns. The cumulative loser-minus-winner spread is **{reversal['long_car_12_pct']:.2f}**, **{reversal['long_car_24_pct']:.2f}**, and **{reversal['long_car_36_pct']:.2f} percentage points** after 12, 24, and 36 months. There are only **{reversal['long_cohorts']} nonoverlapping cohorts** ({', '.join(map(str, reversal['long_formations']))}), so this test has limited precision. The paths are in `results/extension_long_reversal_path.csv` and `results/extension_long_reversal.png`.

The long-horizon method is inspired by De Bondt and Thaler (1985), but uses return deciles rather than their primary 35-stock portfolios. A delisted stock contributes its available delisting return and leaves subsequent monthly averages. January effects and changing market composition remain important interpretation limits.

Sources: https://doi.org/10.1111/j.1540-6261.1990.tb05110.x ; https://doi.org/10.1111/j.1540-6261.1985.tb05004.x
""")
    nonnegative, down = market_state["states"]
    (REPORTS / "05_extension_market_state.md").write_text(f"""# Extension 3: momentum returns across market states

## Independent question and economic motivation

Does the 6/6 winner-minus-loser strategy perform differently after a market decline? The economic concern is that momentum's average return may hide concentrated losses when previous losers recover. This is a modern-sample, exploratory test motivated by [Daniel and Moskowitz, *Momentum Crashes*](https://www.nber.org/papers/w20439), not a replication of their full conditional-risk model.

## Design

Use **exactly the Table I Panel A, J=6/K=6 portfolio**: equal-weighted P10 winners minus P1 losers, with six overlapping monthly vintages. The outcome is its realized month-t return. At the end of t−1, compound the CRSP value-weighted market returns for months t−12 through t−1. Label month t **prior-market-down** when that return is below zero; the other state is nonnegative. This state uses information available before month t. Both portfolio returns and market returns are measured in percentage points in the output CSV.

Estimate `spread_t = a + b × prior_market_down_t + error_t`. The coefficient `b` is the **difference in mean monthly spread** between the prior-down and prior-nonnegative states. Report Newey–West/HAC standard errors with five monthly lags, matching the six-month overlap of the 6/6 portfolio. The 95% intervals and two-sided p value use a normal approximation. This is a descriptive conditional-mean comparison, not a causal or out-of-sample forecast test.

## Results, {market_state['sample_start']}–{market_state['sample_end']}

| Prior 12-month market state | Months | Mean P10−P1 per month | Approx. 95% HAC interval | Share of negative spread months | Worst month |
| --- | ---: | ---: | ---: | ---: | ---: |
| Nonnegative | {nonnegative['months']} | {nonnegative['mean_spread_pct']:.3f}% | [{nonnegative['ci95_low_pct']:.3f}%, {nonnegative['ci95_high_pct']:.3f}%] | {nonnegative['negative_spread_share']:.1%} | {nonnegative['worst_spread_pct']:.2f}% |
| Negative | {down['months']} | {down['mean_spread_pct']:.3f}% | [{down['ci95_low_pct']:.3f}%, {down['ci95_high_pct']:.3f}%] | {down['negative_spread_share']:.1%} | {down['worst_spread_pct']:.2f}% |

The prior-down minus prior-nonnegative mean is **{market_state['down_minus_nonnegative_mean_pct']:+.3f} percentage points per month** (HAC t = **{market_state['difference_hac_t']:.2f}**, approximate two-sided p = **{market_state['difference_two_sided_normal_p']:.3f}**, 95% interval **[{market_state['difference_ci95_low_pct']:.3f}, {market_state['difference_ci95_high_pct']:.3f}] percentage points**). The state split covers **{market_state['months']} months**. The interval includes zero, so this sample does **not** establish a difference at the conventional two-sided 5% level. The [monthly CSV](../results/extension_market_state_monthly.csv), [summary JSON](../results/extension_market_state_summary.json), and [figure](../results/extension_market_state.png) contain the complete outputs.

## Same-month rebound: description only

Among the {down['months']} prior-down months, the market return **during month t** is positive in {market_state['prior_down_current_positive_months']} months. Their average momentum spread is **{market_state['prior_down_current_positive_mean_pct']:.3f}%**; in the other {market_state['prior_down_current_nonpositive_months']} months it is **{market_state['prior_down_current_nonpositive_mean_pct']:.3f}%**. This split is known only after month t ends. It can describe whether momentum losses coincide with market rebounds, but it cannot be used as an ex-ante trading rule.

## Extreme-month sensitivity

| Month | 6/6 P10−P1 | Same-month market return | Prior market down? |
| --- | ---: | ---: | --- |
{worst_rows}

The worst spread month is **{market_state['worst_month']}** ({market_state['worst_month_spread_pct']:.2f}%). If that single month is omitted, the prior-down minus nonnegative mean changes from **{market_state['down_minus_nonnegative_mean_pct']:.2f}** to **{market_state['down_minus_nonnegative_without_worst_pct']:.2f} percentage points**. This is a descriptive sensitivity check selected after observing outcomes, not a replacement estimate or a new significance test.

## Interpretation and limits

The market-state cut is a simple 12-month sign rule, chosen for this extension; it is not Daniel and Moskowitz's panic-state definition. The same stock portfolios appear in every comparison, so the extension changes the conditioning variable rather than the underlying momentum strategy. State frequencies may be uneven, extreme months can dominate a group mean, and the 2002–2025 sample does not provide independent confirmation of a pattern selected from these results. The HAC interval addresses monthly serial correlation but does not remove these design limits. No trading costs or shorting frictions are included.
""")
    (REPORTS / "06_validation_and_ai.md").write_text(f"""# Validation, limitations, and AI use

## Automated checks

{json.dumps(validation, indent=2)}

The complete results and input boundaries are generated by `python3 analysis/run.py` after `analysis/daily_summary.c` has produced its small daily summaries. The original licensed CSV files are read only and must not be committed to a public repository.

## Remaining research checks

Review all table rows for economically implausible returns, compare the daily-derived last-week factor to spot-checked security histories, and review the French Mom tracking difference before submitting. Check the market-state extension's worst months and group sizes before drawing conclusions. No claim of an exact 1965–1989 reproduction is made.

## AI-use disclosure

AI helped interpret the papers and data dictionaries, revise the code and check intermediate output. The generated work was checked against the input date ranges, portfolio counts, original-paper methods, and the public French Mom series. Two material implementation problems were found and corrected during this process:

1. Repeated security-day and security-month records from the source files would have double-counted observations. The workflow now deduplicates these records before aggregation and verifies that repeated monthly returns and daily market returns do not conflict.
2. The first Table III beta grouping incorrectly indexed the prior calendar year, leaving the beta subsamples empty. The year lookup was corrected and all three beta terciles were regenerated and checked for nonmissing results.
""")
    def markdown(headers, rows):
        return "| " + " | ".join(headers) + " |\n| " + " | ".join(["---"]*len(headers)) + " |\n" + "\n".join("| " + " | ".join(map(str, row)) + " |" for row in rows) + "\n"

    sections = []
    for panel in ("A_immediate", "B_skip_5_trading_days"):
        subset = table1[(table1.panel == panel) & (table1.portfolio == "winner_minus_loser")]
        rows = []
        for j in J_SET:
            cells = [subset[(subset.J == j) & (subset.K == k)].iloc[0] for k in K_SET]
            rows.append([j] + [f"{x.mean_monthly_pct:.3f} ({x.t_hac:.2f})" for x in cells])
        sections.append(f"### {panel}: P10−P1 monthly return %, HAC t in parentheses\n\n" +
                        markdown(["J / K", "3", "6", "9", "12"], rows))
    t2_rows = [[r.portfolio, f"{r.post_rank_beta:.3f}", f"{r.avg_formation_market_cap:,.0f}",
                f"{r.avg_stocks_per_vintage:.0f}"] for r in table2.itertuples()]
    t3_spread = table3[table3.portfolio == "P10-P1"]
    t3_rows = [[r.subsample, f"{r.panel_a_mean_pct:.3f} ({r.panel_a_t_hac:.2f})",
                f"{r.panel_b_alpha_pct:.3f} ({r.panel_b_t_hac:.2f})"] for r in t3_spread.itertuples()]
    joint = pd.read_csv(RESULTS / "table3_joint_tests.csv")
    joint_rows = [[r.subsample, f"{r.panel_a_equal_means_chi2_9:.2f}",
                   f"{r.panel_b_zero_alphas_chi2_10:.2f}"] for r in joint.itertuples()]
    replication_path = REPORTS / "02_replication.md"
    with replication_path.open("a") as f:
        f.write("\n## Table I: winner minus loser summary\n\n" + "\n".join(sections))
        f.write("\n## Table II: beta and capitalization\n\n" +
                markdown(["Portfolio", "Post-ranking beta", "Average market cap, source units", "Stocks/vintage"], t2_rows))
        f.write(f"""
### Table II interpretation: original paper versus replication

In the original 1965–1989 sample, the 6/6 loser portfolio P1 has a post-ranking beta of 1.36 and the winner portfolio P10 has a beta of 1.28. The winner-minus-loser beta is therefore −0.08. Both endpoint portfolios have lower average market capitalization than the middle portfolios, and losers are smaller than winners. The authors conclude that the positive spread is unlikely to be explained simply by higher market beta or small-stock exposure. See the original paper's [Table II and discussion](https://www.bauer.uh.edu/rsusmel/phd/jegadeesh-titman93.pdf).

The modern sample has P1 beta **{p1.post_rank_beta:.3f}** and P10 beta **{p10.post_rank_beta:.3f}**, implying a P10−P1 beta of **{p10.post_rank_beta-p1.post_rank_beta:.3f}**. P1 is also smaller than P10 on average (**{p1.avg_formation_market_cap:,.0f}** versus **{p10.avg_formation_market_cap:,.0f}** in the source market-capitalization units); the middle portfolios are larger than both endpoints. With a positive market risk premium, the negative beta of the long–short portfolio works against its positive raw return. This reproduces the original *direction* of the risk and size patterns, with a substantially larger beta difference in the modern sample. Market-capitalization levels across the two data vintages are not directly comparable.

### Original versus modern 6/6 results

| Measure | Original, 1965–1989 | Replication, 2002–2025 |
| --- | ---: | ---: |
| P10−P1 raw monthly return | 0.950% (t = 3.07) | {all_spread.panel_a_mean_pct:.3f}% (HAC t = {all_spread.panel_a_t_hac:.2f}) |
| P10−P1 post-ranking market beta | −0.080 | {p10.post_rank_beta-p1.post_rank_beta:.3f} |
| P10−P1 CAPM alpha, monthly | 1.000% (t = 3.23) | {all_spread.panel_b_alpha_pct:.3f}% (HAC t = {all_spread.panel_b_t_hac:.2f}) |

The original raw return and alpha are from Table III's all-stock rows, and its beta spread is from Table II. Positive market beta compensation cannot explain either sample's positive winner-minus-loser mean. In the modern sample, the raw spread is positive but its HAC t statistic is below a conventional two-sided 5% cutoff, whereas the CAPM alpha's t statistic is above it. A negative beta and a positive average market premium make the market adjustment raise the estimated intercept. The historical and modern t statistics use different samples and need not share identical standard-error calculations; compare their signs and magnitudes with that limitation in mind.
""")
        f.write("\n## Table III: 6/6 winner minus loser across subsamples\n\n" +
                markdown(["Subsample", "Panel A return % (t)", "Panel B alpha % (t)"], t3_rows))
        f.write("\n### Table III joint tests\n\n" +
                markdown(["Subsample", "Equal decile means, HAC Wald χ²(9)", "Zero CAPM alphas, HAC Wald χ²(10)"], joint_rows) +
                "\nThese are HAC Wald analogues of the paper's joint F tests, not numerically identical statistics. Full P1–P10 cells are in `results/table3.csv`.\n")
        f.write("\n## Conclusion\n\nAll 32 modern-sample P10−P1 mean returns in Table I are positive. The 6/6 spread is positive but its HAC t statistic is below conventional 5% thresholds. Table II shows that the loser portfolio has a higher market beta and smaller average market capitalization than the winner portfolio. Table III shows larger raw spreads in small and high-beta groups; risk-adjusted estimates are reported separately. These are descriptive associations, not evidence that momentum is a risk-free or implementable trading profit. The [reversal extension](04_extension_reversal.md) also follows the same six-month-ranked cohorts for 36 months after formation.\n")
    horizon_rows = []
    for h in (12, 24, 36):
        row = event_time["horizons"][h]
        horizon_rows.append([h, f"{row['avg_cumulative_spread_pct']:.2f}",
                             f"{row['hac_t_cumulative']:.2f}",
                             f"[{row['ci95_low_pct']:.2f}, {row['ci95_high_pct']:.2f}]"])
    with (REPORTS / "04_extension_reversal.md").open("a") as f:
        f.write(f"""
## Additional analysis: 6-month ranking followed for 36 months

### Question and design

Do the returns of the **same stocks** selected by a six-month momentum signal persist or reverse over the following three years? This mirrors the question in Jegadeesh and Titman's Table VII. It is separate from the project's 36/36 long-term reversal extension, which ranks stocks on a different, 36-month history.

Each month from {event_time['first_formation']} through {event_time['last_formation']}, eligible NYSE/AMEX common stocks are ranked on their compounded return over the preceding six months. The bottom and top deciles are frozen at formation. Their available monthly returns are followed for event months 1–36, with event month 1 being the first month after the ranking window. Within each frozen decile, each event month's return is equally weighted across stocks with an observed return. A delisted stock contributes its supplied terminal return where available and then leaves subsequent monthly averages. This convention can affect later event months because the remaining holdings change.

Only formations with all 36 follow-up months in the supplied data are included, giving **{event_time['cohorts']} monthly cohorts**. At each horizon, the cumulative spread is the **sum of monthly P10−P1 returns in percentage points**, averaged across cohorts; it is not a compounded buy-and-hold wealth return. This follows the paper's event-time emphasis but does not reproduce its historical data or every delisting and reweighting detail.

### Modern-sample results

| Months after formation | Cumulative P10−P1, percentage points | HAC t | Approx. 95% interval, percentage points |
| ---: | ---: | ---: | ---: |
""" + "\n".join("| " + " | ".join(map(str, row)) + " |" for row in horizon_rows) + f"""

Between months 12 and 36, the mean cumulative spread changes by **{event_time['horizons'][36]['avg_cumulative_spread_pct']-event_time['horizons'][12]['avg_cumulative_spread_pct']:+.2f} percentage points**. The intervals at these horizons include zero in the supplied sample, and this table does not directly test whether the month-12-to-36 change differs from zero. The path should therefore be read as descriptive evidence, not as proof of a later reversal.

The full event-month path is in `results/event_time_6m_path.csv` and `results/event_time_6m_path.png`; cohort-level observations are in `results/event_time_6m_cohorts.csv`. Intervals use Newey–West standard errors across monthly formation cohorts with **h−1 lags** for each h-month cumulative return. Adjacent cohorts share holding months, so they should not be treated as independent observations. The intervals use a normal approximation and are descriptive uncertainty measures.

In the original 1965–1989 sample, Table VII's cumulative spread was about **9.5 percentage points at month 12** and about **4 points at month 36**. Its decline motivated the paper's discussion of partial reversal. Those are historical reference values, not targets for the modern sample. Source: [Jegadeesh and Titman (1993), Table VII](https://www.bauer.uh.edu/rsusmel/phd/jegadeesh-titman93.pdf).
""")
def main():
    RESULTS.mkdir(exist_ok=True); REPORTS.mkdir(exist_ok=True)
    months, permnos, ret, cap, jt, mom_univ, nyse, market, rf, official, tail = load_data()
    sig_a, sig_b = signals(ret, tail)
    table1, monthly_ports, base = make_table1(months, ret, cap, jt, sig_a, sig_b)
    market_state = make_market_state(months, monthly_ports[("A_immediate", 6, 6)], market)
    table1_6_6 = table1[(table1.panel == "A_immediate") & (table1.J == 6) &
                        (table1.K == 6) & (table1.portfolio == "winner_minus_loser")].iloc[0]
    monthly_spread = monthly_ports[("A_immediate", 6, 6)][:, 9] - monthly_ports[("A_immediate", 6, 6)][:, 0]
    if not np.isclose(np.nanmean(monthly_spread)*100, table1_6_6.mean_monthly_pct, atol=1e-8):
        raise ValueError("Market-state monthly spread does not reproduce Table I's 6/6 mean")
    table2 = make_table2(months, base, market)
    beta, beta_count = load_prior_betas(months, permnos)
    table3 = make_table3(months, ret, cap, jt, sig_a[6], market, rf, beta)
    mom = make_mom(months, ret, cap, mom_univ, nyse, official)
    reversal = make_reversal(months, ret, cap, jt, market)
    event_time = make_event_time6(months, ret, cap, jt, sig_a[6])
    validation = dict(monthly_start=int(months[0]), monthly_end=int(months[-1]),
                      securities=int(len(permnos)), stock_observations=int(np.isfinite(ret).sum()),
                      eligible_jt_median=int(np.median(jt.sum(axis=1))),
                      eligible_mom_median=int(np.median(mom_univ.sum(axis=1))),
                      nonmissing_daily_beta_years=beta_count,
                      event_time_6m_cohorts=event_time["cohorts"],
                      event_time_6m_last_formation=event_time["last_formation"],
                      market_state_6_6_months=market_state["months"],
                      market_state_prior_down_months=market_state["states"][1]["months"],
                      riskfree_last_month=int(months[np.flatnonzero(np.isfinite(rf))[-1]]),
                      official_mom_last_common_month=int(months[np.flatnonzero(np.isfinite(official))[-1]]))
    (RESULTS / "validation.json").write_text(json.dumps(validation, indent=2))
    write_reports(table1, table2, table3, mom, reversal, event_time, market_state, validation)
    print("Finished. See results/ and reports/.", flush=True)


if __name__ == "__main__":
    main()
