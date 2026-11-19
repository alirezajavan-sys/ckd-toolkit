"""KDIGO CKD classification by GFR (G) and albuminuria (A) category, and the KDIGO risk heat map.

Reference: KDIGO 2024 Clinical Practice Guideline for the Evaluation and Management of
Chronic Kidney Disease. Kidney Int. 2024;105(4S):S117-S314.
"""

from __future__ import annotations

from .units import acr_to_mgg

GFR_CATEGORIES = {
    "G1": "Normal or high (≥90)",
    "G2": "Mildly decreased (60–89)",
    "G3a": "Mildly to moderately decreased (45–59)",
    "G3b": "Moderately to severely decreased (30–44)",
    "G4": "Severely decreased (15–29)",
    "G5": "Kidney failure (<15)",
}

ALBUMINURIA_CATEGORIES = {
    "A1": "Normal to mildly increased (<30 mg/g; <3 mg/mmol)",
    "A2": "Moderately increased (30–300 mg/g; 3–30 mg/mmol)",
    "A3": "Severely increased (>300 mg/g; >30 mg/mmol)",
}

# KDIGO prognosis heat map: risk of adverse outcomes by G and A category.
_RISK = {
    "G1": ("Low", "Moderate", "High"),
    "G2": ("Low", "Moderate", "High"),
    "G3a": ("Moderate", "High", "Very high"),
    "G3b": ("High", "Very high", "Very high"),
    "G4": ("Very high", "Very high", "Very high"),
    "G5": ("Very high", "Very high", "Very high"),
}


def gfr_category(egfr: float) -> str:
    """KDIGO GFR category (G1–G5) from eGFR in mL/min/1.73 m²."""
    if egfr < 0:
        raise ValueError("eGFR cannot be negative.")
    if egfr >= 90:
        return "G1"
    if egfr >= 60:
        return "G2"
    if egfr >= 45:
        return "G3a"
    if egfr >= 30:
        return "G3b"
    if egfr >= 15:
        return "G4"
    return "G5"


def albuminuria_category(acr: float, unit: str = "mg/g") -> str:
    """KDIGO albuminuria category (A1–A3) from urine ACR."""
    acr_mgg = acr_to_mgg(acr, unit)
    if acr_mgg < 0:
        raise ValueError("ACR cannot be negative.")
    if acr_mgg < 30:
        return "A1"
    if acr_mgg <= 300:
        return "A2"
    return "A3"


def kdigo_risk(g: str, a: str) -> str:
    """KDIGO heat-map risk: 'Low', 'Moderate', 'High' or 'Very high'."""
    if g not in _RISK:
        raise ValueError(f"Unknown GFR category {g!r}.")
    if a not in ALBUMINURIA_CATEGORIES:
        raise ValueError(f"Unknown albuminuria category {a!r}.")
    return _RISK[g][int(a[1]) - 1]


def classify(egfr: float, acr: float, acr_unit: str = "mg/g") -> dict[str, str]:
    """Full KDIGO classification.

    Note: a diagnosis of CKD also requires abnormalities persisting for > 3 months; G1/G2 with A1
    and no other marker of kidney damage is not CKD.
    """
    g = gfr_category(egfr)
    a = albuminuria_category(acr, acr_unit)
    return {
        "gfr_category": g,
        "gfr_description": GFR_CATEGORIES[g],
        "albuminuria_category": a,
        "albuminuria_description": ALBUMINURIA_CATEGORIES[a],
        "stage": f"{g}{a}",
        "risk": kdigo_risk(g, a),
    }
