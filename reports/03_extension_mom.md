# Extension 1: constructed Mom versus Kenneth French Mom

## Question

How closely does a momentum factor built from the supplied CRSP sample track the public Kenneth French US Mom factor?

## Method

At each month-end, require common-stock status and positive prior market equity. Use months t−12 through t−2 to measure prior return, omitting t−1. Set size and momentum breakpoints from NYSE securities: median size and 30th/70th prior-return percentiles. Assign NYSE, AMEX, and NASDAQ stocks to the six size-by-momentum groups and value-weight their month-t returns using t−1 market equity. Construct ½(Small High + Big High) − ½(Small Low + Big Low). This follows the public definition but may differ from the official series because of eligibility details, source vintages, missing-return conventions, and CRSP revisions.

## Results, common 200201–202512 sample

| Measure | Constructed | French official |
| --- | ---: | ---: |
| Average monthly return | 0.137% | 0.134% |
| Monthly volatility | 4.511% | 4.483% |

There are 288 common months. Monthly correlation is **0.999**; constructed minus official mean is **0.002 percentage points**, with **0.214 percentage-point** monthly tracking error. See `results/extension_mom_monthly.csv` and `results/extension_mom_cumulative.png`.

The official factor is a value-weighted 2×3 sort with a one-month skip, while the paper's Tables I–III use equal-weighted deciles with J/K overlapping holding periods. Their return levels should not be treated as the same portfolio.

Source: https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_mom_factor.html
