# Extension 3: momentum returns across market states

## Independent question and economic motivation

Does the 6/6 winner-minus-loser strategy perform differently after a market decline? The economic concern is that momentum's average return may hide concentrated losses when previous losers recover. This is a modern-sample, exploratory test motivated by [Daniel and Moskowitz, *Momentum Crashes*](https://www.nber.org/papers/w20439), not a replication of their full conditional-risk model.

## Design

Use **exactly the Table I Panel A, J=6/K=6 portfolio**: equal-weighted P10 winners minus P1 losers, with six overlapping monthly vintages. The outcome is its realized month-t return. At the end of t−1, compound the CRSP value-weighted market returns for months t−12 through t−1. Label month t **prior-market-down** when that return is below zero; the other state is nonnegative. This state uses information available before month t. Both portfolio returns and market returns are measured in percentage points in the output CSV.

Estimate `spread_t = a + b × prior_market_down_t + error_t`. The coefficient `b` is the **difference in mean monthly spread** between the prior-down and prior-nonnegative states. Report Newey–West/HAC standard errors with five monthly lags, matching the six-month overlap of the 6/6 portfolio. The 95% intervals and two-sided p value use a normal approximation. This is a descriptive conditional-mean comparison, not a causal or out-of-sample forecast test.

## Results, 200201–202512

| Prior 12-month market state | Months | Mean P10−P1 per month | Approx. 95% HAC interval | Share of negative spread months | Worst month |
| --- | ---: | ---: | ---: | ---: | ---: |
| Nonnegative | 223 | 1.456% | [0.892%, 2.020%] | 30.0% | -27.58% |
| Negative | 65 | -2.025% | [-5.729%, 1.679%] | 44.6% | -49.35% |

The prior-down minus prior-nonnegative mean is **-3.481 percentage points per month** (HAC t = **-1.82**, approximate two-sided p = **0.069**, 95% interval **[-7.234, 0.273] percentage points**). The state split covers **288 months**. The interval includes zero, so this sample does **not** establish a difference at the conventional two-sided 5% level. The [monthly CSV](../results/extension_market_state_monthly.csv), [summary JSON](../results/extension_market_state_summary.json), and [figure](../results/extension_market_state.png) contain the complete outputs.

## Same-month rebound: description only

Among the 65 prior-down months, the market return **during month t** is positive in 37 months. Their average momentum spread is **-7.011%**; in the other 28 months it is **4.564%**. This split is known only after month t ends. It can describe whether momentum losses coincide with market rebounds, but it cannot be used as an ex-ante trading rule.

## Extreme-month sensitivity

| Month | 6/6 P10−P1 | Same-month market return | Prior market down? |
| --- | ---: | ---: | --- |
| 200904 | -49.35% | 10.95% | Yes |
| 200903 | -29.32% | 8.68% | Yes |
| 202301 | -27.93% | 7.18% | Yes |
| 202011 | -27.58% | 12.37% | No |
| 200211 | -26.13% | 6.13% | Yes |

The worst spread month is **200904** (-49.35%). If that single month is omitted, the prior-down minus nonnegative mean changes from **-3.48** to **-2.74 percentage points**. This is a descriptive sensitivity check selected after observing outcomes, not a replacement estimate or a new significance test.

## Interpretation and limits

The market-state cut is a simple 12-month sign rule, chosen for this extension; it is not Daniel and Moskowitz's panic-state definition. The same stock portfolios appear in every comparison, so the extension changes the conditioning variable rather than the underlying momentum strategy. State frequencies may be uneven, extreme months can dominate a group mean, and the 2002–2025 sample does not provide independent confirmation of a pattern selected from these results. The HAC interval addresses monthly serial correlation but does not remove these design limits. No trading costs or shorting frictions are included.
