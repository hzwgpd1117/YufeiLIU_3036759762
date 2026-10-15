# Extension 2: reversal and the subsequent path of momentum

## Question

How does return predictability vary across one-month reversal, three-year reversal, and the later performance of stocks selected by a six-month momentum signal?

## Short horizon

Within the same NYSE/AMEX common-stock universe, rank month t−1 returns into deciles and hold for month t. Define reversal as P1−P10. It averages **0.404% per month** (HAC t = **1.07**, 288 months). The monthly path is in `results/extension_short_reversal.csv`. This is a transparent one-month test motivated by Jegadeesh (1990), not a literal reproduction of his forecasting regressions.

## Long horizon

Every third December, rank stocks on the sum of market-adjusted returns over the prior 36 months. Form equal-weighted bottom and top deciles and follow their next 36 monthly returns. The cumulative loser-minus-winner spread is **17.04**, **23.84**, and **19.48 percentage points** after 12, 24, and 36 months. There are only **7 nonoverlapping cohorts** (200212, 200512, 200812, 201112, 201412, 201712, 202012), so this test has limited precision. The paths are in `results/extension_long_reversal_path.csv` and `results/extension_long_reversal.png`.

The long-horizon method is inspired by De Bondt and Thaler (1985), but uses return deciles rather than their primary 35-stock portfolios. A delisted stock contributes its available delisting return and leaves subsequent monthly averages. January effects and changing market composition remain important interpretation limits.

Sources: https://doi.org/10.1111/j.1540-6261.1990.tb05110.x ; https://doi.org/10.1111/j.1540-6261.1985.tb05004.x

## Additional analysis: 6-month ranking followed for 36 months

### Question and design

Do the returns of the **same stocks** selected by a six-month momentum signal persist or reverse over the following three years? This mirrors the question in Jegadeesh and Titman's Table VII. It is separate from the project's 36/36 long-term reversal extension, which ranks stocks on a different, 36-month history.

Each month from 200201 through 202301, eligible NYSE/AMEX common stocks are ranked on their compounded return over the preceding six months. The bottom and top deciles are frozen at formation. Their available monthly returns are followed for event months 1–36, with event month 1 being the first month after the ranking window. Within each frozen decile, each event month's return is equally weighted across stocks with an observed return. A delisted stock contributes its supplied terminal return where available and then leaves subsequent monthly averages. This convention can affect later event months because the remaining holdings change.

Only formations with all 36 follow-up months in the supplied data are included, giving **253 monthly cohorts**. At each horizon, the cumulative spread is the **sum of monthly P10−P1 returns in percentage points**, averaged across cohorts; it is not a compounded buy-and-hold wealth return. This follows the paper's event-time emphasis but does not reproduce its historical data or every delisting and reweighting detail.

### Modern-sample results

| Months after formation | Cumulative P10−P1, percentage points | HAC t | Approx. 95% interval, percentage points |
| ---: | ---: | ---: | ---: |
| 12 | 2.86 | 0.49 | [-8.67, 14.39] |
| 24 | 0.09 | 0.01 | [-17.20, 17.37] |
| 36 | 2.30 | 0.24 | [-16.64, 21.25] |

Between months 12 and 36, the mean cumulative spread changes by **-0.56 percentage points**. The intervals at these horizons include zero in the supplied sample, and this table does not directly test whether the month-12-to-36 change differs from zero. The path should therefore be read as descriptive evidence, not as proof of a later reversal.

The full event-month path is in `results/event_time_6m_path.csv` and `results/event_time_6m_path.png`; cohort-level observations are in `results/event_time_6m_cohorts.csv`. Intervals use Newey–West standard errors across monthly formation cohorts with **h−1 lags** for each h-month cumulative return. Adjacent cohorts share holding months, so they should not be treated as independent observations. The intervals use a normal approximation and are descriptive uncertainty measures.

In the original 1965–1989 sample, Table VII's cumulative spread was about **9.5 percentage points at month 12** and about **4 points at month 36**. Its decline motivated the paper's discussion of partial reversal. Those are historical reference values, not targets for the modern sample. Source: [Jegadeesh and Titman (1993), Table VII](https://www.bauer.uh.edu/rsusmel/phd/jegadeesh-titman93.pdf).
