"""ckdtools: validated, dependency-free kidney function calculators for clinical research.

* eGFR: race-free CKD-EPI 2021 creatinine, 2012 cystatin C and 2021 creatinine–cystatin C equations
* KDIGO G/A categories and prognosis heat map
* 4-variable Kidney Failure Risk Equation (KFRE), 2- and 5-year
"""

from .egfr import egfr_creatinine, egfr_creatinine_cystatin, egfr_cystatin
from .kdigo import albuminuria_category, classify, gfr_category, kdigo_risk
from .kfre import kfre

__all__ = [
    "egfr_creatinine",
    "egfr_cystatin",
    "egfr_creatinine_cystatin",
    "gfr_category",
    "albuminuria_category",
    "kdigo_risk",
    "classify",
    "kfre",
]
__version__ = "1.0.0"
