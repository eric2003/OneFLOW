# Intermediate numerical verification assets

Scripts and CSV tables for formula-level checks that sit between pure gtest
and full CFD cases.

| Asset | Role |
|-------|------|
| `sa_closure_curves.py` | Theory curves for SA fv1(chi), fw(r) |
| `data/sa_*.csv` | Committed reference tables for review/diff |
| C++ `TurbSaPolicy.*Curve*` tests | Production pure API vs closed form |

Do not treat theory CSV alone as proof of solver correctness.
