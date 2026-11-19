"""Race-free CKD-EPI equations for estimated GFR in adults (mL/min/1.73 m²).

References
----------
* Inker LA, Eneanya ND, Coresh J, et al. New creatinine- and cystatin C-based equations
  to estimate GFR without race. N Engl J Med. 2021;385:1737-1749.
* Inker LA, Schmid CH, Tighiouart H, et al. Estimating glomerular filtration rate from
  serum creatinine and cystatin C. N Engl J Med. 2012;367:20-29.

All equations assume IDMS-traceable creatinine and IFCC-standardised cystatin C,
and are validated for adults (age >= 18 years).
"""

from __future__ import annotations

from .units import creatinine_to_mgdl


def _check_common(age: float, sex: str) -> bool:
    """Validate inputs and return True if female."""
    if age < 18:
        raise ValueError("CKD-EPI equations are validated for adults (age >= 18) only.")
    s = sex.strip().lower()
    if s not in {"female", "f", "male", "m"}:
        raise ValueError("sex must be 'female' or 'male'.")
    return s in {"female", "f"}


def _positive(name: str, value: float) -> float:
    if value is None or value <= 0:
        raise ValueError(f"{name} must be a positive number.")
    return float(value)


def egfr_creatinine(scr: float, age: float, sex: str, unit: str = "mg/dL") -> float:
    """2021 CKD-EPI creatinine equation (eGFRcr).

    eGFRcr = 142 × min(Scr/κ, 1)^α × max(Scr/κ, 1)^-1.200 × 0.9938^age × 1.012 [if female]
    κ = 0.7 (F) / 0.9 (M);  α = -0.241 (F) / -0.302 (M)

    >>> round(egfr_creatinine(1.0, 60, "female"), 1)
    64.5
    """
    female = _check_common(age, sex)
    scr = _positive("Serum creatinine", creatinine_to_mgdl(scr, unit))
    kappa, alpha = (0.7, -0.241) if female else (0.9, -0.302)
    ratio = scr / kappa
    egfr = 142 * min(ratio, 1) ** alpha * max(ratio, 1) ** -1.200 * 0.9938**age
    return egfr * 1.012 if female else egfr


def egfr_cystatin(scys: float, age: float, sex: str) -> float:
    """2012 CKD-EPI cystatin C equation (eGFRcys); cystatin C in mg/L.

    eGFRcys = 133 × min(Scys/0.8, 1)^-0.499 × max(Scys/0.8, 1)^-1.328 × 0.996^age × 0.932 [if female]
    """
    female = _check_common(age, sex)
    scys = _positive("Serum cystatin C", scys)
    ratio = scys / 0.8
    egfr = 133 * min(ratio, 1) ** -0.499 * max(ratio, 1) ** -1.328 * 0.996**age
    return egfr * 0.932 if female else egfr


def egfr_creatinine_cystatin(
    scr: float, scys: float, age: float, sex: str, unit: str = "mg/dL"
) -> float:
    """2021 CKD-EPI creatinine–cystatin C equation (eGFRcr-cys). Preferred when both are available.

    eGFRcr-cys = 135 × min(Scr/κ, 1)^α × max(Scr/κ, 1)^-0.544
                     × min(Scys/0.8, 1)^-0.323 × max(Scys/0.8, 1)^-0.778
                     × 0.9961^age × 0.963 [if female]
    κ = 0.7 (F) / 0.9 (M);  α = -0.219 (F) / -0.144 (M)
    """
    female = _check_common(age, sex)
    scr = _positive("Serum creatinine", creatinine_to_mgdl(scr, unit))
    scys = _positive("Serum cystatin C", scys)
    kappa, alpha = (0.7, -0.219) if female else (0.9, -0.144)
    r_cr, r_cys = scr / kappa, scys / 0.8
    egfr = (
        135
        * min(r_cr, 1) ** alpha
        * max(r_cr, 1) ** -0.544
        * min(r_cys, 1) ** -0.323
        * max(r_cys, 1) ** -0.778
        * 0.9961**age
    )
    return egfr * 0.963 if female else egfr
