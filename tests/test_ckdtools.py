import math
import warnings

import pytest

from ckdtools import (
    albuminuria_category,
    classify,
    egfr_creatinine,
    egfr_creatinine_cystatin,
    egfr_cystatin,
    gfr_category,
    kdigo_risk,
    kfre,
)
from ckdtools.cli import main

# --------------------------------------------------------------------------- eGFR


def test_egfr_creatinine_hand_calculated_female():
    # 60-year-old woman, Scr 1.0 mg/dL: 142 × (1/0.7)^-1.2 × 0.9938^60 × 1.012
    expected = 142 * (1.0 / 0.7) ** -1.2 * 0.9938**60 * 1.012
    assert egfr_creatinine(1.0, 60, "female") == pytest.approx(expected)
    assert round(egfr_creatinine(1.0, 60, "female"), 1) == 64.5


def test_egfr_creatinine_hand_calculated_male_low_scr():
    # 40-year-old man, Scr 0.7 mg/dL (below κ = 0.9): 142 × (0.7/0.9)^-0.302 × 0.9938^40
    expected = 142 * (0.7 / 0.9) ** -0.302 * 0.9938**40
    assert egfr_creatinine(0.7, 40, "male") == pytest.approx(expected)


def test_egfr_creatinine_continuous_at_knot():
    # At Scr = κ both spline arms give the same value.
    below = egfr_creatinine(0.9 - 1e-9, 50, "male")
    above = egfr_creatinine(0.9 + 1e-9, 50, "male")
    assert below == pytest.approx(above, rel=1e-6)


def test_egfr_creatinine_umol_equals_mgdl():
    assert egfr_creatinine(88.42, 55, "male", unit="umol/L") == pytest.approx(egfr_creatinine(1.0, 55, "male"))


def test_egfr_decreases_with_creatinine_and_age():
    assert egfr_creatinine(2.0, 60, "male") < egfr_creatinine(1.0, 60, "male")
    assert egfr_creatinine(1.0, 80, "male") < egfr_creatinine(1.0, 40, "male")


def test_egfr_cystatin_hand_calculated():
    expected = 133 * (1.2 / 0.8) ** -1.328 * 0.996**65 * 0.932
    assert egfr_cystatin(1.2, 65, "female") == pytest.approx(expected)


def test_egfr_combined_hand_calculated():
    # Man, 70 y, Scr 1.5 mg/dL, cystatin C 1.4 mg/L (both above their knots)
    expected = 135 * (1.5 / 0.9) ** -0.544 * (1.4 / 0.8) ** -0.778 * 0.9961**70
    assert egfr_creatinine_cystatin(1.5, 1.4, 70, "male") == pytest.approx(expected)


@pytest.mark.parametrize("bad", [dict(age=17), dict(sex="other"), dict(scr=0), dict(scr=-1)])
def test_egfr_rejects_invalid_input(bad):
    kwargs = dict(scr=1.0, age=50, sex="female") | bad
    with pytest.raises(ValueError):
        egfr_creatinine(**kwargs)


# --------------------------------------------------------------------------- KDIGO


@pytest.mark.parametrize(
    "egfr,cat",
    [(120, "G1"), (90, "G1"), (89.9, "G2"), (60, "G2"), (59.9, "G3a"), (45, "G3a"),
     (44.9, "G3b"), (30, "G3b"), (29.9, "G4"), (15, "G4"), (14.9, "G5"), (0, "G5")],
)
def test_gfr_category_boundaries(egfr, cat):
    assert gfr_category(egfr) == cat


@pytest.mark.parametrize(
    "acr,unit,cat",
    [(10, "mg/g", "A1"), (29.9, "mg/g", "A1"), (30, "mg/g", "A2"), (300, "mg/g", "A2"),
     (301, "mg/g", "A3"), (2.9, "mg/mmol", "A1"), (3.5, "mg/mmol", "A2"), (35, "mg/mmol", "A3")],
)
def test_albuminuria_category(acr, unit, cat):
    assert albuminuria_category(acr, unit) == cat


def test_kdigo_heat_map_corners():
    assert kdigo_risk("G1", "A1") == "Low"
    assert kdigo_risk("G3a", "A1") == "Moderate"
    assert kdigo_risk("G3a", "A2") == "High"
    assert kdigo_risk("G3b", "A2") == "Very high"
    assert kdigo_risk("G5", "A1") == "Very high"


def test_classify():
    c = classify(38, 450)
    assert c["stage"] == "G3bA3"
    assert c["risk"] == "Very high"


# --------------------------------------------------------------------------- KFRE


def test_kfre_at_cohort_means():
    # At the development-cohort means the linear predictor is 0 (except the sex term).
    age, egfr, acr = 70.36, 36.11, math.exp(5.137)
    lp_female = 0.2467 * (0 - 0.5642)
    assert kfre(age, "female", egfr, acr, years=5) == pytest.approx(1 - 0.9365 ** math.exp(lp_female))
    assert kfre(age, "female", egfr, acr, years=2, region="north_american") == pytest.approx(
        1 - 0.9750 ** math.exp(lp_female)
    )


def test_kfre_direction_of_effects():
    base = dict(age=65, sex="male", egfr=30, acr=300)
    r = kfre(**base)
    assert kfre(**(base | {"egfr": 20})) > r  # lower eGFR -> higher risk
    assert kfre(**(base | {"acr": 1000})) > r  # more albuminuria -> higher risk
    assert kfre(**(base | {"age": 80})) < r  # older age -> lower risk of kidney failure (competing death)
    assert kfre(**(base | {"sex": "female"})) < r
    assert kfre(**base, years=5) > kfre(**base, years=2)
    assert kfre(**base, region="north_american") > r


def test_kfre_mg_mmol_equals_mg_g():
    assert kfre(60, "female", 40, 30, acr_unit="mg/mmol") == pytest.approx(kfre(60, "female", 40, 30 * 8.84))


def test_kfre_warns_outside_validated_range():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        kfre(50, "male", 75, 100)
    assert any("eGFR < 60" in str(x.message) for x in w)


def test_kfre_rejects_bad_horizon():
    with pytest.raises(ValueError):
        kfre(60, "male", 30, 100, years=3)


# --------------------------------------------------------------------------- CLI


def test_cli_runs(capsys):
    assert main(["--age", "67", "--sex", "male", "--scr", "2.1", "--acr", "450"]) == 0
    out = capsys.readouterr().out
    assert "CKD-EPI 2021 creatinine" in out
    assert "G3bA3" in out
    assert "5-year risk" in out
