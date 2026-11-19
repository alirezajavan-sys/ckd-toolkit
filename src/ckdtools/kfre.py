"""4-variable Kidney Failure Risk Equation (KFRE).

Predicts the probability of treated kidney failure (dialysis or transplant) at 2 and 5 years
in adults with CKD G3–G5 (eGFR < 60 mL/min/1.73 m²).

    risk = 1 - S0 ^ exp( -0.2201 × (age/10 - 7.036)
                         + 0.2467 × (male - 0.5642)
                         - 0.5567 × (eGFR/5 - 7.222)
                         + 0.4510 × (ln ACR[mg/g] - 5.137) )

Baseline survival S0:
    North American:      2 y 0.9750, 5 y 0.9240
    Non-North American:  2 y 0.9832, 5 y 0.9365  (multinational recalibration, 2016)

References
----------
* Tangri N, Stevens LA, Griffith J, et al. A predictive model for progression of chronic
  kidney disease to kidney failure. JAMA. 2011;305(15):1553-1559.
* Tangri N, Grams ME, Levey AS, et al. Multinational assessment of accuracy of equations
  for predicting risk of kidney failure: a meta-analysis. JAMA. 2016;315(2):164-174.
"""

from __future__ import annotations

import math
import warnings

from .units import acr_to_mgg

BASELINE_SURVIVAL = {
    "north_american": {2: 0.9750, 5: 0.9240},
    "non_north_american": {2: 0.9832, 5: 0.9365},
}


def linear_predictor(age: float, sex: str, egfr: float, acr_mgg: float) -> float:
    male = 1.0 if sex.strip().lower() in {"male", "m"} else 0.0
    return (
        -0.2201 * (age / 10 - 7.036)
        + 0.2467 * (male - 0.5642)
        - 0.5567 * (egfr / 5 - 7.222)
        + 0.4510 * (math.log(acr_mgg) - 5.137)
    )


def kfre(
    age: float,
    sex: str,
    egfr: float,
    acr: float,
    acr_unit: str = "mg/g",
    years: int = 5,
    region: str = "non_north_american",
) -> float:
    """Probability (0–1) of kidney failure within `years` (2 or 5).

    `region` selects the baseline-survival calibration: 'non_north_american' (default; Europe,
    Middle East, Asia, Oceania, Latin America) or 'north_american'.
    """
    if years not in (2, 5):
        raise ValueError("KFRE is defined for 2- or 5-year horizons only.")
    if region not in BASELINE_SURVIVAL:
        raise ValueError(f"region must be one of {list(BASELINE_SURVIVAL)}.")
    if sex.strip().lower() not in {"female", "f", "male", "m"}:
        raise ValueError("sex must be 'female' or 'male'.")
    acr_mgg = acr_to_mgg(acr, acr_unit)
    if acr_mgg <= 0:
        raise ValueError("ACR must be positive (log-transformed in the equation).")
    if egfr >= 60:
        warnings.warn("KFRE was developed and validated for eGFR < 60 mL/min/1.73 m².", stacklevel=2)
    s0 = BASELINE_SURVIVAL[region][years]
    return 1 - s0 ** math.exp(linear_predictor(age, sex, egfr, acr_mgg))
