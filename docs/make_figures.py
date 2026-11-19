"""Generate the README figures:  python docs/make_figures.py  (needs matplotlib, numpy)."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ckdtools import egfr_creatinine, kdigo_risk, kfre

OUT = Path(__file__).resolve().parent
plt.rcParams.update({"figure.dpi": 150, "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

# ---- 1. KDIGO heat map ---------------------------------------------------------
g_cats = ["G1", "G2", "G3a", "G3b", "G4", "G5"]
g_lab = ["G1  ≥90", "G2  60–89", "G3a 45–59", "G3b 30–44", "G4  15–29", "G5  <15"]
a_cats = ["A1", "A2", "A3"]
a_lab = ["A1\n<30 mg/g", "A2\n30–300 mg/g", "A3\n>300 mg/g"]
colors = {"Low": "#4caf50", "Moderate": "#ffeb3b", "High": "#ff9800", "Very high": "#e53935"}

fig, ax = plt.subplots(figsize=(5.2, 4.6))
for i, g in enumerate(g_cats):
    for j, a in enumerate(a_cats):
        r = kdigo_risk(g, a)
        ax.add_patch(plt.Rectangle((j, len(g_cats) - 1 - i), 1, 1, color=colors[r], ec="white", lw=2))
        ax.text(j + 0.5, len(g_cats) - 0.5 - i, r, ha="center", va="center", fontsize=9)
ax.set_xlim(0, 3)
ax.set_ylim(0, 6)
ax.set_xticks(np.arange(3) + 0.5, a_lab)
ax.set_yticks(np.arange(6) + 0.5, g_lab[::-1])
ax.xaxis.tick_top()
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_title("KDIGO CKD prognosis heat map\n(`ckdtools.kdigo_risk`)", pad=34, fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "kdigo_heatmap.png")
plt.close(fig)

# ---- 2. eGFR vs creatinine and KFRE vs ACR ---------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
scr = np.linspace(0.4, 4, 300)
for age, ls in [(40, "-"), (70, "--")]:
    for sex, c in [("female", "#cf222e"), ("male", "#1f6feb")]:
        axes[0].plot(scr, [egfr_creatinine(s, age, sex) for s in scr], color=c, ls=ls, lw=1.8,
                     label=f"{sex.capitalize()}, {age} y")
for y in (90, 60, 45, 30, 15):
    axes[0].axhline(y, color="grey", lw=0.6, ls=":")
axes[0].set(xlabel="Serum creatinine (mg/dL)", ylabel="eGFR (mL/min/1.73 m²)",
            title="CKD-EPI 2021 (race-free) eGFR", ylim=(0, 140))
axes[0].legend(frameon=False, fontsize=8)

acr = np.logspace(0.5, 3.5, 200)
for egfr, c in [(50, "#2da44e"), (35, "#bf8700"), (20, "#cf222e")]:
    axes[1].plot(acr, [100 * kfre(65, "male", egfr, a, years=5) for a in acr], color=c, lw=1.8,
                 label=f"eGFR {egfr}")
axes[1].axhline(5, color="grey", lw=0.8, ls="--", label="5% (NICE referral threshold)")
axes[1].set_xscale("log")
axes[1].set(xlabel="Urine ACR (mg/g, log scale)", ylabel="5-year kidney failure risk (%)",
            title="KFRE: 65-year-old man (non-North American)")
axes[1].legend(frameon=False, fontsize=8, title="mL/min/1.73 m²", title_fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "egfr_kfre_curves.png")
plt.close(fig)
print("figures written to", OUT)
