"""Command-line interface.

    ckdtools --age 67 --sex male --scr 2.1 --acr 450
    ckdtools --age 54 --sex female --scr 160 --scr-unit umol/L --scys 1.6 --acr 35 --acr-unit mg/mmol
"""

from __future__ import annotations

import argparse
import warnings

from .egfr import egfr_creatinine, egfr_creatinine_cystatin, egfr_cystatin
from .kdigo import classify
from .kfre import kfre


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="ckdtools", description="eGFR, KDIGO staging and KFRE risk.")
    p.add_argument("--age", type=float, required=True, help="years (adults only)")
    p.add_argument("--sex", choices=["female", "male"], required=True)
    p.add_argument("--scr", type=float, help="serum creatinine")
    p.add_argument("--scr-unit", default="mg/dL", choices=["mg/dL", "umol/L"])
    p.add_argument("--scys", type=float, help="serum cystatin C, mg/L")
    p.add_argument("--acr", type=float, help="urine albumin-to-creatinine ratio")
    p.add_argument("--acr-unit", default="mg/g", choices=["mg/g", "mg/mmol"])
    p.add_argument("--region", default="non_north_american", choices=["non_north_american", "north_american"],
                   help="KFRE baseline calibration")
    a = p.parse_args(argv)

    if a.scr is None and a.scys is None:
        p.error("provide --scr and/or --scys")

    results: list[tuple[str, float]] = []
    if a.scr is not None:
        results.append(("CKD-EPI 2021 creatinine", egfr_creatinine(a.scr, a.age, a.sex, a.scr_unit)))
    if a.scys is not None:
        results.append(("CKD-EPI 2012 cystatin C", egfr_cystatin(a.scys, a.age, a.sex)))
    if a.scr is not None and a.scys is not None:
        results.append(
            ("CKD-EPI 2021 creatinine-cystatin C", egfr_creatinine_cystatin(a.scr, a.scys, a.age, a.sex, a.scr_unit))
        )

    print("eGFR (mL/min/1.73 m²)")
    for name, value in results:
        print(f"  {name:<36} {value:6.1f}")
    best_name, best = results[-1]

    if a.acr is not None:
        c = classify(best, a.acr, a.acr_unit)
        print(f"\nKDIGO (using {best_name})")
        print(f"  Stage   {c['stage']}  — {c['gfr_description']}; {c['albuminuria_description']}")
        print(f"  Risk    {c['risk']}")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            r2 = kfre(a.age, a.sex, best, a.acr, a.acr_unit, 2, a.region)
            r5 = kfre(a.age, a.sex, best, a.acr, a.acr_unit, 5, a.region)
        print(f"\nKidney Failure Risk Equation ({a.region.replace('_', '-')} calibration)")
        print(f"  2-year risk  {100 * r2:5.1f}%")
        print(f"  5-year risk  {100 * r5:5.1f}%")
        for msg in w:
            print(f"  Note: {msg.message}")

    print("\nFor research and education only; not a substitute for clinical judgement.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
