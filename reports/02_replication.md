# Tables I–III: modern-sample replication

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

The 6/6 immediate-formation P10−P1 portfolio averages **0.670% per month** (HAC t = **1.32**, 288 months). The complete 32-strategy grid is in `results/table1.csv`; `results/table2.csv` and `results/table3.csv` contain the other assigned tables.

## Interpretation and limits

The original paper's 12/3 spread is 1.31% per month without a gap and 1.49% with a one-week gap. These are reference values from the 1965–1989 sample, not targets for the 2000–2025 data. CIZ monthly returns compound daily returns and can differ from the historical CRSP format. The five-trading-observation interval is a transparent approximation to a calendar week. Security universe, data revisions, sample period, and inference method also affect comparisons. Read the full tables before interpreting sign or significance of any specification.

Source: https://doi.org/10.1111/j.1540-6261.1993.tb04702.x

## Table I: winner minus loser summary

### A_immediate: P10−P1 monthly return %, HAC t in parentheses

| J / K | 3 | 6 | 9 | 12 |
| --- | --- | --- | --- | --- |
| 3 | 0.342 (0.96) | 0.485 (1.44) | 0.510 (1.39) | 0.402 (1.16) |
| 6 | 0.636 (1.32) | 0.670 (1.32) | 0.596 (1.16) | 0.425 (0.89) |
| 9 | 0.831 (1.45) | 0.698 (1.20) | 0.537 (0.94) | 0.360 (0.66) |
| 12 | 0.471 (0.82) | 0.386 (0.65) | 0.290 (0.48) | 0.159 (0.28) |

### B_skip_5_trading_days: P10−P1 monthly return %, HAC t in parentheses

| J / K | 3 | 6 | 9 | 12 |
| --- | --- | --- | --- | --- |
| 3 | 0.508 (1.48) | 0.526 (1.61) | 0.547 (1.53) | 0.434 (1.29) |
| 6 | 0.716 (1.53) | 0.677 (1.35) | 0.602 (1.20) | 0.437 (0.93) |
| 9 | 0.821 (1.47) | 0.673 (1.17) | 0.537 (0.94) | 0.372 (0.70) |
| 12 | 0.556 (0.98) | 0.404 (0.68) | 0.301 (0.50) | 0.176 (0.31) |

## Table II: beta and capitalization

| Portfolio | Post-ranking beta | Average market cap, source units | Stocks/vintage |
| --- | --- | --- | --- |
| P1 | 1.890 | 1,806,074 | 160 |
| P2 | 1.483 | 5,502,145 | 160 |
| P3 | 1.310 | 9,297,750 | 160 |
| P4 | 1.194 | 11,361,691 | 160 |
| P5 | 1.105 | 12,690,280 | 160 |
| P6 | 1.042 | 13,789,680 | 159 |
| P7 | 1.010 | 14,309,656 | 159 |
| P8 | 1.045 | 13,485,011 | 159 |
| P9 | 1.088 | 11,469,891 | 159 |
| P10 | 1.255 | 5,722,706 | 159 |

### Table II interpretation: original paper versus replication

In the original 1965–1989 sample, the 6/6 loser portfolio P1 has a post-ranking beta of 1.36 and the winner portfolio P10 has a beta of 1.28. The winner-minus-loser beta is therefore −0.08. Both endpoint portfolios have lower average market capitalization than the middle portfolios, and losers are smaller than winners. The authors conclude that the positive spread is unlikely to be explained simply by higher market beta or small-stock exposure. See the original paper's [Table II and discussion](https://www.bauer.uh.edu/rsusmel/phd/jegadeesh-titman93.pdf).

The modern sample has P1 beta **1.890** and P10 beta **1.255**, implying a P10−P1 beta of **-0.635**. P1 is also smaller than P10 on average (**1,806,074** versus **5,722,706** in the source market-capitalization units); the middle portfolios are larger than both endpoints. With a positive market risk premium, the negative beta of the long–short portfolio works against its positive raw return. This reproduces the original *direction* of the risk and size patterns, with a substantially larger beta difference in the modern sample. Market-capitalization levels across the two data vintages are not directly comparable.

### Original versus modern 6/6 results

| Measure | Original, 1965–1989 | Replication, 2002–2025 |
| --- | ---: | ---: |
| P10−P1 raw monthly return | 0.950% (t = 3.07) | 0.670% (HAC t = 1.32) |
| P10−P1 post-ranking market beta | −0.080 | -0.635 |
| P10−P1 CAPM alpha, monthly | 1.000% (t = 3.23) | 1.130% (HAC t = 2.85) |

The original raw return and alpha are from Table III's all-stock rows, and its beta spread is from Table II. Positive market beta compensation cannot explain either sample's positive winner-minus-loser mean. In the modern sample, the raw spread is positive but its HAC t statistic is below a conventional two-sided 5% cutoff, whereas the CAPM alpha's t statistic is above it. A negative beta and a positive average market premium make the market adjustment raise the estimated intercept. The historical and modern t statistics use different samples and need not share identical standard-error calculations; compare their signs and magnitudes with that limitation in mind.

## Table III: 6/6 winner minus loser across subsamples

| Subsample | Panel A return % (t) | Panel B alpha % (t) |
| --- | --- | --- |
| All | 0.670 (1.32) | 1.130 (2.85) |
| Size_small | 0.988 (1.56) | 1.536 (3.03) |
| Size_mid | 0.577 (1.25) | 0.967 (2.50) |
| Size_big | 0.315 (1.04) | 0.592 (2.24) |
| Beta_low | 0.379 (0.95) | 0.663 (2.00) |
| Beta_mid | 0.404 (0.98) | 0.749 (2.32) |
| Beta_high | 1.227 (1.90) | 1.772 (3.35) |

### Table III joint tests

| Subsample | Equal decile means, HAC Wald χ²(9) | Zero CAPM alphas, HAC Wald χ²(10) |
| --- | --- | --- |
| All | 14.17 | 24.95 |
| Size_small | 6.51 | 29.00 |
| Size_mid | 13.81 | 21.90 |
| Size_big | 9.46 | 21.86 |
| Beta_low | 20.77 | 47.89 |
| Beta_mid | 10.56 | 18.82 |
| Beta_high | 26.87 | 34.31 |

These are HAC Wald analogues of the paper's joint F tests, not numerically identical statistics. Full P1–P10 cells are in `results/table3.csv`.

## Conclusion

All 32 modern-sample P10−P1 mean returns in Table I are positive. The 6/6 spread is positive but its HAC t statistic is below conventional 5% thresholds. Table II shows that the loser portfolio has a higher market beta and smaller average market capitalization than the winner portfolio. Table III shows larger raw spreads in small and high-beta groups; risk-adjusted estimates are reported separately. These are descriptive associations, not evidence that momentum is a risk-free or implementable trading profit. The [reversal extension](04_extension_reversal.md) also follows the same six-month-ranked cohorts for 36 months after formation.
