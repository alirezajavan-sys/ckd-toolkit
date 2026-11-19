"""Unit conversions for common kidney laboratory values."""

from __future__ import annotations

# Creatinine: 1 mg/dL = 88.42 µmol/L
CREATININE_UMOL_PER_MGDL = 88.42

# Albumin-to-creatinine ratio: 1 mg/mmol = 8.84 mg/g
ACR_MGG_PER_MGMMOL = 8.84


def creatinine_to_mgdl(value: float, unit: str = "mg/dL") -> float:
    """Convert serum creatinine to mg/dL. Accepts 'mg/dL' or 'umol/L' / 'µmol/L'."""
    unit = unit.replace("µ", "u").lower()
    if unit == "mg/dl":
        return float(value)
    if unit == "umol/l":
        return float(value) / CREATININE_UMOL_PER_MGDL
    raise ValueError(f"Unknown creatinine unit: {unit!r}. Use 'mg/dL' or 'umol/L'.")


def acr_to_mgg(value: float, unit: str = "mg/g") -> float:
    """Convert urine albumin-to-creatinine ratio to mg/g. Accepts 'mg/g' or 'mg/mmol'."""
    unit = unit.lower()
    if unit == "mg/g":
        return float(value)
    if unit == "mg/mmol":
        return float(value) * ACR_MGG_PER_MGMMOL
    raise ValueError(f"Unknown ACR unit: {unit!r}. Use 'mg/g' or 'mg/mmol'.")
