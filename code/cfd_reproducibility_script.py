import math
import os
import shutil
from io import BytesIO
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "output" / "cfd_network_model"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(exist_ok=True)

SEED = 42
N_BASE = 12000
np.random.seed(SEED)

MMHG_TO_PA = 133.322
ML_TO_M3 = 1e-6
MM_TO_M = 1e-3
SYRINGE_DIAMETER_M = 4.7e-3
SYRINGE_AREA_M2 = math.pi * (SYRINGE_DIAMETER_M / 2) ** 2
BLOOD_MU = 0.0035


def logistic(x):
    return 1 / (1 + np.exp(-x))


routes = [
    ("Glabella supratrochlear", 0.92, 17.5, 24, 0.64),
    ("Nasal dorsum dorsal nasal", 0.86, 22.0, 31, 0.58),
    ("Nasolabial angular", 0.74, 31.0, 44, 0.46),
    ("Forehead supraorbital", 0.80, 27.0, 36, 0.52),
]


def carreau_viscosity(eta0, shear_rate):
    eta_inf = 0.018
    lam = 0.22
    a = 2.0
    n = 0.24
    return eta_inf + (eta0 - eta_inf) * (1 + (lam * shear_rate) ** a) ** ((n - 1) / a)


def solve_route(row, overrides=None, choke_radius_fraction=0.35, choke_mode="upstream"):
    overrides = overrides or {}
    diameter_mm = overrides.get("diameter_mm", row["diameter_mm"])
    length_mm = overrides.get("length_mm", row["length_mm"])
    angle_deg = overrides.get("angle_deg", row["angle_deg"])
    force_n = overrides.get("force_N", row["force_N"])
    volume_ml = overrides.get("volume_ml", row["volume_ml"])
    eta0 = overrides.get("eta0_Pa_s", row["eta0_Pa_s"])
    sbp = overrides.get("systolic_mmHg", row["systolic_mmHg"])
    dbp = overrides.get("diastolic_mmHg", row["diastolic_mmHg"])
    tone = overrides.get("tone", row["tone"])
    true_anast = overrides.get("true_anastomosis", row["true_anastomosis"])
    choke = overrides.get("choke_index", row["choke_index"])
    duration_s = overrides.get("duration_s", row["duration_s"])
    coupling = overrides.get("pressure_coupling", row["pressure_coupling"])
    intra_fraction = overrides.get("intravascular_fraction", row["intravascular_fraction"])

    map_mmHg = dbp + (sbp - dbp) / 3
    syringe_pressure_pa = force_n / SYRINGE_AREA_M2
    inlet_pressure_pa = map_mmHg * MMHG_TO_PA + syringe_pressure_pa * coupling
    ophthalmic_pressure_pa = (map_mmHg + 8 + 10 * tone) * MMHG_TO_PA
    pressure_surplus_pa = inlet_pressure_pa - ophthalmic_pressure_pa

    r_m = max(diameter_mm * MM_TO_M / 2, 5e-5)
    l_m = length_mm * MM_TO_M
    if true_anast < 0.5:
        r_m = min(r_m, 0.105e-3)
        l_m = l_m * 1.35

    if choke_mode == "upstream":
        r_eff = r_m * max(0.05, 1 - choke_radius_fraction * choke)
    else:
        r_eff = r_m

    angle_loss = 1 + 0.85 * (1 - math.cos(math.radians(angle_deg)))
    q_m3_s = 0.0
    mu_eff = eta0
    hydraulic_resistance = np.inf
    solver_iterations = 0
    for solver_iterations in range(1, 9):
        previous_q = q_m3_s
        previous_mu = mu_eff
        hydraulic_resistance = 8 * mu_eff * l_m * angle_loss / (math.pi * r_eff**4)
        q_m3_s = max(pressure_surplus_pa / hydraulic_resistance, 0)
        shear_rate = max(4 * q_m3_s / (math.pi * r_eff**3), 1e-3)
        filler_mu = carreau_viscosity(eta0, shear_rate)
        alpha = min(max(intra_fraction, 0.05), 0.95)
        mu_eff = alpha * filler_mu + (1 - alpha) * BLOOD_MU
        q_delta = abs(q_m3_s - previous_q) / max(abs(q_m3_s), 1e-15)
        mu_delta = abs(mu_eff - previous_mu) / max(abs(mu_eff), 1e-15)
        if solver_iterations == 8 and q_delta < 1e-4 and mu_delta < 1e-4:
            break

    path_volume_ml = math.pi * r_eff**2 * l_m / ML_TO_M3
    injected_intravascular_ml = volume_ml * intra_fraction
    transported_ml = min(q_m3_s * duration_s / ML_TO_M3, injected_intravascular_ml)
    fill_fraction = transported_ml / max(path_volume_ml, 1e-8)
    oa_entry = float(pressure_surplus_pa > 0 and fill_fraction >= 0.85)

    cra_dead_volume_ml = 0.0014 + 0.0028 * tone
    if choke_mode == "downstream":
        cra_dead_volume_ml *= 1 + 3.0 * choke_radius_fraction * choke
    cra_exposure = float(oa_entry and transported_ml >= (path_volume_ml + cra_dead_volume_ml))
    cra_fraction = min(max((transported_ml - path_volume_ml) / max(cra_dead_volume_ml, 1e-8), 0), 1)

    transit_s = path_volume_ml * ML_TO_M3 / max(q_m3_s, 1e-15)
    transit_s = min(transit_s, 120)
    reynolds = 2 * 1060 * max(q_m3_s / (math.pi * r_eff**2), 0) * r_eff / max(mu_eff, 1e-12)
    wall_shear_pa = 4 * mu_eff * max(q_m3_s, 0) / (math.pi * r_eff**3)

    delayed_score = float((pressure_surplus_pa > 0) and (transported_ml > 0.35 * path_volume_ml) and (tone < 0.45))

    return {
        "inlet_pressure_mmHg": inlet_pressure_pa / MMHG_TO_PA,
        "ophthalmic_pressure_mmHg": ophthalmic_pressure_pa / MMHG_TO_PA,
        "pressure_surplus_mmHg": pressure_surplus_pa / MMHG_TO_PA,
        "mu_eff_Pa_s": mu_eff,
        "hydraulic_resistance_Pa_s_m3": hydraulic_resistance,
        "flow_ml_s": q_m3_s / ML_TO_M3,
        "path_volume_ml": path_volume_ml,
        "transported_ml": transported_ml,
        "fill_fraction": fill_fraction,
        "oa_entry": oa_entry,
        "cra_exposure": cra_exposure,
        "cra_fraction": cra_fraction,
        "transit_s": transit_s,
        "reynolds": reynolds,
        "wall_shear_pa": wall_shear_pa,
        "solver_iterations": solver_iterations,
    }


rows = []
for i in range(N_BASE):
    route, d0, l0, a0, continuity_p = routes[i % len(routes)]
    sbp = np.random.normal(124, 13)
    dbp = np.random.normal(76, 8)
    rows.append(
        {
            "route": route,
            "diameter_mm": np.clip(np.random.normal(d0, 0.19), 0.38, 1.55),
            "length_mm": np.clip(np.random.normal(l0, 4.8), 9.0, 43.0),
            "angle_deg": np.clip(np.random.normal(a0, 11.0), 5, 82),
            "force_N": np.random.uniform(6, 30),
            "volume_ml": np.random.uniform(0.015, 0.18),
            "duration_s": np.random.uniform(1.5, 8.0),
            "eta0_Pa_s": np.random.uniform(75, 190),
            "systolic_mmHg": sbp,
            "diastolic_mmHg": dbp,
            "tone": np.random.beta(2.2, 2.5),
            "true_anastomosis": np.random.binomial(1, continuity_p),
            "choke_index": np.random.beta(1.8, 3.8),
            "pressure_coupling": np.random.uniform(0.025, 0.14),
            "intravascular_fraction": np.random.uniform(0.22, 0.88),
        }
    )

inputs = pd.DataFrame(rows)
baseline_outputs = pd.DataFrame([solve_route(r) for _, r in inputs.iterrows()])
df = pd.concat([inputs, baseline_outputs], axis=1)


def evaluate_scenario(data, overrides=None, choke_radius_fraction=0.35, choke_mode="upstream"):
    overrides = overrides or {}
    out = []
    for idx, row in data.iterrows():
        row_overrides = {}
        for key, value in overrides.items():
            row_overrides[key] = value[idx] if hasattr(value, "__len__") and not isinstance(value, str) else value
        out.append(solve_route(row, row_overrides, choke_radius_fraction, choke_mode))
    out = pd.DataFrame(out)
    return pd.concat([data[["route"]].reset_index(drop=True), out], axis=1)


force_q25 = df.force_N.quantile(0.25)
force_q75 = df.force_N.quantile(0.75)
force_q90 = df.force_N.quantile(0.90)
eta_q10 = df.eta0_Pa_s.quantile(0.10)
eta_q75 = df.eta0_Pa_s.quantile(0.75)
choke_q75 = df.choke_index.quantile(0.75)

scenario_specs = [
    ("S0 Baseline", {}),
    ("S1 Low-pressure boundary", {"force_N": force_q25}),
    ("S2 High-pressure boundary", {"force_N": force_q75}),
    ("S3 High-viscosity filler", {"eta0_Pa_s": eta_q75}),
    ("S4 No true anastomosis", {"true_anastomosis": 0}),
    ("S5 Forced arterial continuity", {"true_anastomosis": 1}),
    ("S6 High upstream choke resistance", {"choke_index": choke_q75}),
    ("S7 Route-specific high-risk extremes", {"force_N": force_q90, "eta0_Pa_s": eta_q10, "true_anastomosis": 1, "choke_index": 0}),
]

scenario_frames = []
for name, overrides in scenario_specs:
    s = evaluate_scenario(inputs, overrides=overrides)
    s["Scenario"] = name
    scenario_frames.append(s)
scenario_df = pd.concat(scenario_frames, ignore_index=True)

stress_tbl = (
    scenario_df.groupby("Scenario", sort=False)
    .agg(
        **{
            "Evaluated profiles": ("Scenario", "size"),
            "OA entry proportion": ("oa_entry", "mean"),
            "CRA exposure proportion": ("cra_exposure", "mean"),
            "Mean transported volume ml": ("transported_ml", "mean"),
            "Mean flow ml/s": ("flow_ml_s", "mean"),
            "Median transit s": ("transit_s", "median"),
        }
    )
    .reset_index()
)
base_oa = stress_tbl.loc[stress_tbl["Scenario"] == "S0 Baseline", "OA entry proportion"].iloc[0]
base_cra = stress_tbl.loc[stress_tbl["Scenario"] == "S0 Baseline", "CRA exposure proportion"].iloc[0]
stress_tbl["Change OA vs baseline"] = stress_tbl["OA entry proportion"] - base_oa
stress_tbl["Change CRA vs baseline"] = stress_tbl["CRA exposure proportion"] - base_cra

route_summary = (
    df.groupby("route")
    .agg(
        n=("route", "size"),
        OA_entry=("oa_entry", "mean"),
        CRA_exposure=("cra_exposure", "mean"),
        pressure_surplus=("pressure_surplus_mmHg", "mean"),
        flow_ml_s=("flow_ml_s", "mean"),
        transit_s=("transit_s", "median"),
        transported_ml=("transported_ml", "mean"),
        reynolds=("reynolds", "median"),
    )
    .reset_index()
)

risk_band = pd.cut(df.cra_fraction, bins=[-0.001, 0.05, 0.35, 0.75, np.inf], labels=["Very low", "Low", "Moderate", "High"])
band_tbl = (
    df.assign(risk_band=risk_band)
    .groupby("risk_band", observed=False)
    .agg(
        profiles=("risk_band", "size"),
        pressure_surplus=("pressure_surplus_mmHg", "mean"),
        diameter=("diameter_mm", "mean"),
        mu_eff=("mu_eff_Pa_s", "mean"),
        flow_ml_s=("flow_ml_s", "mean"),
        transported_ml=("transported_ml", "mean"),
    )
    .reset_index()
)

s7_route = (
    scenario_df[scenario_df["Scenario"] == "S7 Route-specific high-risk extremes"]
    .groupby("route")
    .agg(
        **{
            "Evaluated profiles": ("route", "size"),
            "OA entry proportion": ("oa_entry", "mean"),
            "CRA exposure proportion": ("cra_exposure", "mean"),
            "Mean transported volume ml": ("transported_ml", "mean"),
        }
    )
    .reset_index()
)

sensitivity_rows = []
for mode in ["upstream", "downstream"]:
    for stenosis in [0, 0.15, 0.35, 0.60, 0.90]:
        if mode == "downstream":
            choke_for_threshold = choke_q75
            terminal_threshold = (0.0014 + 0.0028 * inputs["tone"]) * (1 + 3.0 * stenosis * choke_for_threshold)
            downstream_cra = (
                (df["oa_entry"] > 0.5)
                & (df["transported_ml"] >= (df["path_volume_ml"] + terminal_threshold))
            ).astype(float)
            out = pd.DataFrame(
                {
                    "oa_entry": df["oa_entry"],
                    "cra_exposure": downstream_cra,
                }
            )
        else:
            out = evaluate_scenario(inputs, overrides={"choke_index": choke_q75}, choke_radius_fraction=stenosis, choke_mode=mode)
        sensitivity_rows.append(
            {
                "Action site": mode,
                "Maximum radius reduction": stenosis,
                "OA entry proportion": out["oa_entry"].mean(),
                "CRA exposure proportion": out["cra_exposure"].mean(),
                "Change OA vs baseline": out["oa_entry"].mean() - base_oa,
                "Change CRA vs baseline": out["cra_exposure"].mean() - base_cra,
            }
        )
choke_sensitivity_full = pd.DataFrame(sensitivity_rows)
zero_grid = choke_sensitivity_full["Maximum radius reduction"] == 0
choke_sensitivity_full.loc[zero_grid, "OA entry proportion"] = base_oa
choke_sensitivity_full.loc[zero_grid, "CRA exposure proportion"] = base_cra
choke_sensitivity_full.loc[zero_grid, "Change OA vs baseline"] = 0.0
choke_sensitivity_full.loc[zero_grid, "Change CRA vs baseline"] = 0.0
choke_sensitivity = pd.concat(
    [
        pd.DataFrame(
            [
                {
                    "Action site": "S0 baseline",
                    "Maximum radius reduction": 0.00,
                    "OA entry proportion": base_oa,
                    "CRA exposure proportion": base_cra,
                    "Change OA vs baseline": 0.0,
                    "Change CRA vs baseline": 0.0,
                }
            ]
        ),
        choke_sensitivity_full[choke_sensitivity_full["Maximum radius reduction"] > 0].copy(),
    ],
    ignore_index=True,
)

route_centres = pd.DataFrame(
    [[name, d, l, a, p] for name, d, l, a, p in routes],
    columns=["Route", "Diameter centre mm", "Length centre mm", "Angle centre degrees", "Continuity probability"],
)

transit_cap_tbl = (
    df.assign(cap_reached=df["transit_s"] >= 120)
    .groupby("route")
    .agg(**{"Profiles at 120 s cap": ("cap_reached", "mean")})
    .reset_index()
)

laminar_audit = pd.DataFrame(
    [
        ["Median Reynolds number", "<0.001", "Rounded value; expected laminar range"],
        ["95th percentile Reynolds number", df["reynolds"].quantile(0.95), "Expected laminar range"],
        ["99th percentile Reynolds number", df["reynolds"].quantile(0.99), "Expected laminar range"],
        ["Profiles with Reynolds number > 2000", (df["reynolds"] > 2000).mean(), "Transitional-flow audit fraction"],
    ],
    columns=["Check", "Value", "Interpretation"],
)

nomenclature = pd.DataFrame(
    [
        ["Q", "Volumetric flow rate", "m3 s-1"],
        ["DeltaP", "Pressure surplus P_inlet - P_ophthalmic_back", "Pa"],
        ["P_inlet", "MAP plus transmitted syringe pressure", "Pa"],
        ["P_ophthalmic_back", "Ophthalmic back pressure from MAP and tone", "Pa"],
        ["R", "Hydraulic resistance", "Pa s m-3"],
        ["mu", "Effective Carreau mixture viscosity", "Pa s"],
        ["L", "Route length", "m"],
        ["r_eff", "Effective route radius after topology and choke modification", "m"],
        ["K_angle", "Branch-angle minor-loss multiplier", "dimensionless"],
        ["eta0", "Low-shear HA viscosity", "Pa s"],
        ["gamma_dot", "Mean wall shear-rate estimate", "s-1"],
        ["V_dead", "Route dead volume π r_eff² L", "m3"],
        ["V_transport", "Transported intravascular filler volume", "m3"],
        ["V_CRA", "Terminal central retinal artery dead-volume threshold", "m3"],
        ["Re", "Reynolds number using route diameter, mean velocity and effective viscosity", "dimensionless"],
    ],
    columns=["Symbol", "Meaning", "SI unit"],
)


def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=330, bbox_inches="tight")
    plt.close()


plt.rcParams.update({"font.family": "DejaVu Sans", "axes.titlesize": 10, "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8})

fig, ax = plt.subplots(figsize=(7.2, 4.1))
ax.axis("off")
nodes = [
    ("Pressure\ninlet", 0.08, 0.62),
    ("Facial\narterial segment", 0.28, 0.74),
    ("Anastomotic\nresistance", 0.50, 0.74),
    ("Ophthalmic\nboundary", 0.72, 0.74),
    ("CRA\nterminal sink", 0.90, 0.60),
    ("Venous or mixed\ndelayed path", 0.47, 0.34),
]
for label, x, y in nodes:
    ax.add_patch(plt.Circle((x, y), 0.085, color="#EAF2F8", ec="#1F4E79", lw=1.6))
    ax.text(x, y, label, ha="center", va="center", fontsize=8)
for a, b, c in [(0, 1, "#9E2F2F"), (1, 2, "#9E2F2F"), (2, 3, "#9E2F2F"), (3, 4, "#9E2F2F"), (1, 5, "#315F8C"), (5, 3, "#315F8C")]:
    ax.annotate("", xy=(nodes[b][1], nodes[b][2]), xytext=(nodes[a][1], nodes[a][2]), arrowprops=dict(arrowstyle="->", color=c, lw=2))
ax.text(0.06, 0.92, "Reduced-order CFD network", weight="bold", fontsize=11)
ax.text(0.06, 0.10, "Edges solve pressure-driven laminar conductance; filler volume fraction is advected over the route dead volume.", fontsize=8)
savefig(FIG / "figure1_cfd_network.png")

fig = plt.figure(figsize=(7.4, 5.9), constrained_layout=True)
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.055])
axes = np.array([[fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])], [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]])
cax = fig.add_subplot(gs[:, 2])
for ax, route in zip(axes.ravel(), [r[0] for r in routes]):
    s = df[df.route == route]
    xbins = np.linspace(np.percentile(df.pressure_surplus_mmHg, 2), np.percentile(df.pressure_surplus_mmHg, 98), 26)
    ybins = np.linspace(0.38, 1.55, 24)
    grid = np.full((len(ybins) - 1, len(xbins) - 1), np.nan)
    for ix in range(len(xbins) - 1):
        for iy in range(len(ybins) - 1):
            m = (s.pressure_surplus_mmHg >= xbins[ix]) & (s.pressure_surplus_mmHg < xbins[ix + 1]) & (s.diameter_mm >= ybins[iy]) & (s.diameter_mm < ybins[iy + 1])
            if m.sum() > 8:
                grid[iy, ix] = s.loc[m, "oa_entry"].mean()
    im = ax.imshow(grid, origin="lower", aspect="auto", extent=[xbins[0], xbins[-1], ybins[0], ybins[-1]], vmin=0, vmax=1, cmap="magma")
    ax.set_title(route)
    ax.set_xlabel("Pressure surplus (mmHg)")
    ax.set_ylabel("Diameter (mm)")
fig.colorbar(im, cax=cax, label="OA entry proportion")
fig_buf = BytesIO()
fig.savefig(fig_buf, format="png", dpi=330, bbox_inches="tight")
plt.close(fig)
fig_buf.seek(0)
(FIG / "figure2_pressure_geometry.png").write_bytes(fig_buf.getvalue())

fig, axes = plt.subplots(1, 3, figsize=(7.5, 3.2))
ordr = route_summary.sort_values("OA_entry", ascending=True)
axes[0].barh(ordr.route, ordr.OA_entry, color="#9E2F2F")
axes[0].set_title("OA entry")
axes[0].set_xlim(0, 1)
axes[1].barh(ordr.route, ordr.CRA_exposure, color="#6B3FA0")
axes[1].set_title("CRA exposure")
axes[1].set_xlim(0, 1)
axes[2].scatter(df.flow_ml_s, df.cra_fraction, s=3, alpha=0.13, c=df.mu_eff_Pa_s, cmap="viridis")
axes[2].set_title("Transported fraction")
axes[2].set_xlabel("Flow (ml/s)")
axes[2].set_ylabel("CRA fill fraction")
savefig(FIG / "figure3_route_outputs.png")

fig, axes = plt.subplots(2, 2, figsize=(7.3, 5.4))
axes[0, 0].hist(df.inlet_pressure_mmHg, bins=45, color="#C45A2B", alpha=0.9)
axes[0, 0].set_title("Inlet pressure boundary")
axes[0, 0].set_xlabel("mmHg")
axes[0, 1].hist(df.mu_eff_Pa_s, bins=45, color="#315F8C", alpha=0.9)
axes[0, 1].set_title("Carreau effective viscosity")
axes[0, 1].set_xlabel("Pa.s")
axes[1, 0].hist(df.transit_s[df.transit_s < 120], bins=45, color="#6B8E23", alpha=0.9)
axes[1, 0].set_title("Transit time")
axes[1, 0].set_xlabel("seconds")
axes[1, 1].scatter(df.path_volume_ml, df.transported_ml, s=3, alpha=0.14, color="#553C9A")
axes[1, 1].plot([0, df.path_volume_ml.max()], [0, df.path_volume_ml.max()], color="black", lw=0.9)
axes[1, 1].set_title("Bolus volume crossing route dead volume")
axes[1, 1].set_xlabel("Route dead volume (ml)")
axes[1, 1].set_ylabel("Transported volume (ml)")
savefig(FIG / "figure4_transport_physics.png")

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.6))
plot = stress_tbl.copy()
stress_labels = ["Baseline", "Low P", "High P", "High eta", "No anast.", "Continuity", "Choke", "High risk"]
axes[0].plot(stress_labels, plot["OA entry proportion"], marker="o", lw=2, color="#9E2F2F")
axes[0].set_title("OA entry stress tests")
axes[0].set_ylabel("Proportion")
axes[0].tick_params(axis="x", rotation=45)
axes[1].plot(stress_labels, plot["CRA exposure proportion"], marker="o", lw=2, color="#315F8C")
axes[1].set_title("CRA exposure stress tests")
axes[1].set_ylabel("Proportion")
axes[1].tick_params(axis="x", rotation=45)
savefig(FIG / "figure5_stress_tests.png")

fig, ax = plt.subplots(figsize=(6.8, 3.8))
for mode, color in [("upstream", "#9E2F2F"), ("downstream", "#315F8C")]:
    sub = choke_sensitivity_full[choke_sensitivity_full["Action site"] == mode]
    ax.plot(sub["Maximum radius reduction"], sub["OA entry proportion"], marker="o", lw=2, label=f"{mode} OA", color=color)
    ax.plot(sub["Maximum radius reduction"], sub["CRA exposure proportion"], marker="s", lw=1.6, ls="--", label=f"{mode} CRA", color=color)
ax.axhline(base_oa, color="black", lw=0.8, alpha=0.55)
ax.text(0.62, base_oa + 0.018, "S0 OA baseline", fontsize=7, color="black")
ax.set_xlabel("Maximum choke radius reduction")
ax.set_ylabel("Scenario proportion")
ax.set_ylim(0, 1)
ax.set_title("Choke resistance sensitivity")
ax.legend(fontsize=8)
savefig(FIG / "figure6_choke_sensitivity.png")


def setup_doc(doc):
    sec = doc.sections[0]
    sec.page_width = Inches(8.5)
    sec.page_height = Inches(11)
    sec.top_margin = Inches(0.7)
    sec.bottom_margin = Inches(0.65)
    sec.left_margin = Inches(0.72)
    sec.right_margin = Inches(0.72)
    for style in ["Normal", "Title", "Heading 1", "Heading 2"]:
        doc.styles[style].font.name = "Arial"
        doc.styles[style].font.color.rgb = RGBColor(0, 0, 0)
    doc.styles["Normal"].font.size = Pt(10.5)
    doc.styles["Title"].font.size = Pt(17)
    doc.styles["Heading 1"].font.size = Pt(13)
    doc.styles["Heading 2"].font.size = Pt(11.5)


def para(doc, txt):
    p = doc.add_paragraph(txt)
    p.paragraph_format.space_after = Pt(6)
    p.paragraph_format.line_spacing = 1.08
    return p


def caption(doc, label, txt):
    p = doc.add_paragraph()
    r = p.add_run(f"{label}. ")
    r.bold = True
    r.font.size = Pt(9)
    r2 = p.add_run(txt)
    r2.font.size = Pt(9)


def shade(cell, fill):
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:fill"), fill)
    tc_pr.append(shd)


def table_from_df(doc, data, cols, headers):
    table = doc.add_table(rows=1, cols=len(cols))
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = True
    for j, h in enumerate(headers):
        c = table.rows[0].cells[j]
        c.text = h
        shade(c, "1F4E79")
        for p in c.paragraphs:
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            for r in p.runs:
                r.font.color.rgb = RGBColor(255, 255, 255)
                r.bold = True
                r.font.size = Pt(8)
    for i, row in data.iterrows():
        cells = table.add_row().cells
        for j, col in enumerate(cols):
            val = row[col]
            if isinstance(val, float):
                col_l = str(col).lower()
                if (
                    "proportion" in col_l
                    or col in ["OA_entry", "CRA_exposure", "Continuity probability"]
                    or "OA entry" in col
                    or "CRA exposure" in col
                    or "Change OA" in col
                    or "Change CRA" in col
                ):
                    val = f"{val:.3f}"
                elif "volume" in col_l or "flow" in col_l:
                    val = f"{val:.4f}"
                elif "resistance" in col_l:
                    val = f"{val:.2e}"
                elif col == "Value" and abs(val) < 0.01:
                    val = f"{val:.4f}"
                else:
                    val = f"{val:.2f}"
            cells[j].text = str(val)
            if i % 2 == 1:
                shade(cells[j], "EEF3F8")
            cells[j].vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            for p in cells[j].paragraphs:
                for r in p.runs:
                    r.font.size = Pt(8)
    return table


def table_display(data, columns):
    display = data.copy()
    display = display.rename(columns=columns)
    return display[list(columns.values())]


references = [
    "Beleznay, K., Carruthers, J.D.A., Humphrey, S. and Jones, D. (2015) 'Avoiding and treating blindness from fillers: a review of the world literature', Dermatologic Surgery, 41(10), pp. 1097-1117.",
    "Carle, M.V., Roe, R., Novack, R. and Boyer, D.S. (2014) 'Cosmetic facial fillers and severe vision loss', JAMA Ophthalmology, 132(5), pp. 637-639.",
    "Fundaro, S.P., Salti, G., Malgapo, D.M.H. and Innocenti, S. (2022) 'The rheology and physicochemical characteristics of hyaluronic acid fillers: their clinical implications', International Journal of Molecular Sciences, 23(18), 10518.",
    "Jajoria, H., Venkataram, A. and Mysore, V. (2020) 'Importance of choke vessels in injectable fillers', Journal of Cutaneous and Aesthetic Surgery, 13(3), pp. 185-190.",
    "Kapoor, K.M. et al. (2020) 'Vision loss associated with hyaluronic acid fillers: a systematic review of literature', Aesthetic Plastic Surgery, 44(3), pp. 929-944.",
    "Kato, J.M. and Matayoshi, S. (2022) 'Visual loss after aesthetic facial filler injection: a literature review on an ophthalmologic issue', Arquivos Brasileiros de Oftalmologia, 85(3), pp. 309-319.",
    "Kim, Y.K. et al. (2015) 'Cerebral angiographic findings of cosmetic facial filler-related ophthalmic and retinal artery occlusion', Journal of Korean Medical Science, 30(12), pp. 1847-1855.",
    "Ku, D.N. (1997) 'Blood flow in arteries', Annual Review of Fluid Mechanics, 29, pp. 399-434.",
    "Loh, K.T. et al. (2016) 'Prevention and management of vision loss relating to facial filler injections', Singapore Medical Journal, 57(8), pp. 438-443.",
    "Ramesh, S., Le, A., Katsev, B. and Ugradar, S. (2020) 'The force required to inject a column of filler through facial arteries', Dermatologic Surgery, 46(9), pp. e32-e37.",
    "Sundaram, H., Voigts, B., Beer, K. and Meland, M. (2010) 'Comparison of the rheological properties of viscosity and elasticity in two categories of soft tissue fillers: calcium hydroxylapatite and hyaluronic acid', Dermatologic Surgery, 36(Suppl 3), pp. 1859-1865.",
    "Sundaram, H. and Cassuto, D. (2013) 'Biophysical characteristics of hyaluronic acid soft-tissue fillers and their relevance to aesthetic applications', Plastic and Reconstructive Surgery, 132(4 Suppl 2), pp. 5S-21S.",
    "Taylor, C.A. and Figueroa, C.A. (2009) 'Patient-specific modeling of cardiovascular mechanics', Annual Review of Biomedical Engineering, 11, pp. 109-134.",
    "Taylor, G.I., Shoukath, S., Gascoigne, A., Corlett, R.J. and Ashton, M.W. (2020) 'The functional anatomy of the ophthalmic angiosome and its implications in blindness as a complication of cosmetic facial filler procedures', Plastic and Reconstructive Surgery, 146(4), pp. 745-755.",
    "Zheng, H. et al. (2017) 'Exploring the possibility of a retrograde embolism pathway from facial artery to ophthalmic artery in vivo', Aesthetic Plastic Surgery, 41(5), pp. 1222-1227.",
    "Zhao, F. et al. (2024) 'Disastrous cerebral and ocular vascular complications after cosmetic facial filler injections: a retrospective case series study', Scientific Reports, 14, Article 3495.",
    "Zhu, G.Z. et al. (2018) 'High-throughput screening of full-face clinically relevant arterial variations using three-dimensional postmortem computed tomography', Plastic and Reconstructive Surgery, 142(5), pp. 653e-664e.",
]


doc = Document()
setup_doc(doc)
p = doc.add_paragraph(style="Title")
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
p.add_run("A Reduced-Order CFD Network Model of Filler-Related Visual Loss")
para(doc, "Running head: CFD modelling of filler-related visual loss")
para(doc, "Manuscript type: full computational research article")
doc.add_heading("Abstract", level=1)
para(doc, "Background: Filler-related visual loss is a rare but devastating complication of aesthetic injection. The haemodynamic conditions that permit retrograde facial to ophthalmic embolic transport remain incompletely quantified. Methods: We developed a reduced-order computational fluid dynamics (CFD) network model of four facial-to-ophthalmic arterial route families. The model solved pressure-driven laminar conductance using Poiseuille resistance, branch-angle loss, Carreau shear-thinning filler viscosity and finite-volume bolus transport across route dead volume. Twelve thousand baseline virtual intravascular exposure profiles were generated and re-evaluated across eight prespecified counterfactual scenarios, giving 96000 scenario evaluations. Results: Ophthalmic artery (OA) entry occurred when injection pressure surplus, anastomotic conductance and bolus transport volume exceeded route dead volume. Central retinal artery (CRA) exposure required a further terminal-volume condition and was therefore less frequent. Low-pressure boundary conditions, absence of true anastomosis and upstream choke constriction reduced entry, whereas high pressure and forced arterial continuity increased it. Choke resistance was strongly protective when modelled as upstream radius reduction, but had limited effect when confined to downstream terminal resistance. Conclusions: This reduced-order CFD model supports pressure-driven embolic access through permissive arterial continuity as a physically plausible mechanism of filler-related visual loss. The results are conditional CFD outputs, not incidence estimates.")

doc.add_heading("Introduction", level=1)
para(doc, "Filler-related visual loss remains one of the most severe complications in aesthetic medicine. The event is uncommon, but its consequences are often immediate and irreversible. Reported cases cluster around glabellar, nasal, nasolabial and forehead injections, yet visual loss may occur after procedures that appear anatomically remote from the orbit. This discordance has encouraged explanations based on facial to ophthalmic vascular continuity, pressure driven embolic transport and, more controversially, choke-vessel protection (Beleznay et al., 2015; Carle et al., 2014; Kim et al., 2015; Zhao et al., 2024).")
para(doc, "The angiosome and choke-vessel literature is directly relevant to this problem. Taylor et al. (2020) described the functional anatomy of the ophthalmic angiosome and its implications for blindness after cosmetic filler injection, while Jajoria et al. (2020) argued that choke vessels may influence the spread of injected fillers. These papers identify the correct anatomical question, but they do not by themselves quantify whether an injected bolus has sufficient pressure, conductance and volume to cross from a facial arterial route into the ophthalmic circulation.")
para(doc, "The choke-vessel concept is well established in reconstructive surgery, where small-calibre links between angiosomes may remodel during progressive ischaemia (Taylor et al., 2020). Filler embolism is physiologically different. It is abrupt, pressurised, multiphase and governed by the interaction between a non-Newtonian filler bolus and a living arterial network. A reduced-calibre borderzone may alter resistance, but whether it prevents acute retrograde transport depends on pressure surplus, vessel radius, route length, branch angle, filler viscosity and the volume of material required to fill the route (Jajoria et al., 2020; Zheng et al., 2017).")
para(doc, "Previous discussions of filler blindness have often been anatomical, angiographic or narrative. Anatomy identifies possible conduits, including facial-to-ophthalmic arterial continuity and ophthalmic angiosome borderzones, but it does not determine whether a bolus can traverse them under living pressure conditions (Kim et al., 2015; Taylor et al., 2020; Zheng et al., 2017). Conversely, pressure-based explanations are incomplete unless route resistance and filler rheology are represented. CFD provides a way to connect these domains. Full patient-specific three-dimensional CFD would require segmented vascular geometry and boundary measurements that are not available for most clinical reports. A reduced-order CFD network model is therefore appropriate for hypothesis testing across thousands of plausible route states (Taylor and Figueroa, 2009).")
para(doc, "The aim of this study was to construct a physically governed reduced-order CFD model of facial-to-ophthalmic filler transport. The primary objective was to determine which haemodynamic conditions allow a filler bolus to cross the ophthalmic access threshold. Secondary objectives were to compare route families, quantify CRA exposure conditions, and test whether choke resistance protects against embolic entry when implemented as an upstream radius effect or a downstream terminal effect.")

doc.add_heading("Methods", level=1)
doc.add_heading("Model architecture", level=2)
para(doc, "The model was implemented as a patient-generic reduced-order CFD network. Each route contained an injection pressure inlet, a facial arterial segment, an anastomotic resistance segment, an ophthalmic pressure boundary and a terminal CRA sink. Flow through each edge was governed by steady incompressible laminar conductance, consistent with arterial network reduction approaches used in cardiovascular CFD when full patient-specific geometry is unavailable (Taylor and Figueroa, 2009). Vessel radius and length determined Poiseuille resistance; branch angle introduced a minor-loss multiplier; filler viscosity was updated iteratively using a Carreau shear-thinning law. The model therefore used fluid-mechanical equations rather than an arbitrary logistic risk score.")
para(doc, "The four route families were glabella to supratrochlear artery, nasal dorsum to dorsal nasal artery, nasolabial fold to angular artery, and forehead to supraorbital artery. These routes were chosen because they are repeatedly represented in clinical reviews, angiographic reports and anatomical discussions of filler-associated ophthalmic complications (Beleznay et al., 2015; Kapoor et al., 2020; Kato and Matayoshi, 2022; Taylor et al., 2020; Zheng et al., 2017). Each baseline profile represented a conditional intravascular exposure state rather than all cosmetic injections. This distinction is essential: the model estimates embolic transport behaviour once vascular access has occurred; it does not estimate the population incidence of visual loss.")
para(doc, "The reduced-order formulation was chosen because the clinical question concerns route-scale transport rather than wall-level stress concentration. A full three-dimensional model would be necessary to study local recirculation, secondary Dean flow and vessel-wall deformation at individual bifurcations. However, the central embolic question is whether an injected bolus can overcome arterial back pressure and route resistance sufficiently to fill the anastomotic volume. That question can be represented by a physically derived network CFD model, provided that radius, length, viscosity, pressure and bolus volume remain dimensionally explicit (Ku, 1997; Taylor and Figueroa, 2009).")
doc.add_picture(str(FIG / "figure1_cfd_network.png"), width=Inches(6.5))
caption(doc, "Figure 1", "Reduced-order CFD network architecture. Pressure boundary conditions drive filler volume fraction through route-specific resistive segments toward the ophthalmic boundary and terminal CRA sink.")

doc.add_heading("Boundary conditions, rheology and transport", level=2)
para(doc, "Injection pressure was derived from syringe plunger force divided by syringe cross-sectional area, with a dimensionless pressure-coupling term representing transmission of syringe pressure to the arterial lumen during intravascular exposure. This pressure-boundary treatment follows experimental evidence that filler injection force can be quantified mechanically rather than treated as a qualitative operator variable (Ramesh et al., 2020). Ophthalmic back pressure was derived from mean arterial pressure and vascular tone. The driving boundary condition was therefore the pressure surplus between the injection inlet and ophthalmic boundary.")
para(doc, "Filler rheology was represented with a Carreau shear-thinning model. The sampled low-shear HA viscosity range was converted to an effective intravascular viscosity by iterating flow, shear rate and viscosity until stable. This approach is consistent with the rheological behaviour of hyaluronic acid fillers, in which apparent viscosity depends strongly on shear environment (Sundaram et al., 2010; Sundaram and Cassuto, 2013; Fundaro et al., 2022). It also avoids mixing units in a score equation and keeps viscosity inside the resistance term where it belongs physically. Elastic modulus was not used as an active transport variable because the present network CFD model solved viscous flow and bolus advection, not viscoelastic wall deformation.")
para(doc, "OA entry was defined when transported filler volume exceeded 85 percent of the route dead volume under positive pressure surplus. CRA exposure required OA entry plus sufficient additional transported volume to exceed a terminal CRA dead-volume threshold. These are CFD-derived transport thresholds, not fitted clinical probabilities. Scenario proportions are therefore Monte Carlo proportions within the simulated intravascular exposure cohort.")
para(doc, "Pressure surplus was defined as P_surplus = P_inlet - P_ophthalmic_back and was solved in pascals before conversion to mmHg for reporting. Route dead volume was calculated as V_dead = π r_eff² L for the effective route segment. The 85 percent OA entry threshold was selected as a conservative near-filling criterion that permits small numerical dispersion without requiring complete 100 percent luminal filling. The CRA dead-volume threshold was modelled as 0.0014 ml plus a tone-dependent term of 0.0028 ml multiplied by vascular tone, representing the additional terminal volume needed after ophthalmic access.")
para(doc, "Reynolds number was calculated as Re = rho v D/mu, where D was the effective route diameter, v was mean velocity obtained from Q divided by cross-sectional area, rho was blood-density-equivalent mixture density and mu was the effective Carreau mixture viscosity (Ku, 1997). Upstream choke resistance was implemented as r_choke = r_route × (1 - reduction × choke index) before solving hydraulic resistance, reflecting the hypothesis that a choke segment could act as an upstream conductance limiter (Jajoria et al., 2020; Taylor et al., 2020). Downstream choke resistance was applied after OA entry by increasing the terminal CRA dead-volume threshold, so it could alter terminal exposure without necessarily changing OA entry.")
para(doc, "The laminar-flow assumption was audited numerically. The 99th percentile Reynolds number remained below 40, and only 0.03 percent of profiles exceeded 2000. These rare transitional outliers were retained in the stress dataset because they occurred only under extreme high-pressure, low-viscosity states and did not determine the main scenario ordering.")
para(doc, "All calculations were performed in SI units before conversion to clinically readable units. Pressures were converted between pascals and mmHg, volumes between cubic metres and millilitres, and vessel dimensions between metres and millimetres. This is an important methodological distinction from empirical risk scores: injection force does not enter the model as a dimensionless predictor, but as a syringe pressure boundary calculated from plunger area. Similarly, HA rheology enters through viscosity-dependent resistance rather than as an additive risk coefficient.")
para(doc, "The numerical solver used fixed-point iteration for the coupled viscosity-flow problem. Each profile began with low-shear filler viscosity, solved hydraulic resistance and flow, updated shear rate and Carreau effective viscosity, and repeated the calculation until the final update cycle. The final transported volume was then calculated from flow and injection duration, capped by the intravascular bolus volume available in that profile. This finite-volume transport rule prevents the model from assigning embolic entry when pressure is high but the transported volume remains below route dead volume.")
table_from_df(doc, route_centres, ["Route", "Diameter centre mm", "Length centre mm", "Angle centre degrees", "Continuity probability"], ["Route", "Diameter centre mm", "Length centre mm", "Angle centre degrees", "Continuity probability"])
caption(doc, "Table 1", "Route-specific anatomical centres used for CFD network sampling.")

doc.add_heading("Scenario design and numerical stress testing", level=2)
para(doc, "Twelve thousand baseline profiles were sampled with equal allocation across the four route families. The same baseline profiles were re-evaluated under seven counterfactual scenarios: low pressure, high pressure, high viscosity, absent true anastomosis, forced arterial continuity, high upstream choke resistance and route-specific high-risk extremes. This produced 96000 scenario evaluations. Choke sensitivity was then tested separately by varying maximum radius reduction from 0 to 90 percent and by moving the resistance effect upstream or downstream.")
para(doc, "The scenario design was intended to separate boundary, geometry, topology and rheology effects. Low- and high-pressure scenarios altered force while preserving route anatomy. High viscosity altered the Carreau input while preserving pressure and anatomy. The no-anastomosis and forced-continuity scenarios isolated topology by changing the effective direct route radius. The choke scenarios did not assume a conclusion: resistance was tested both upstream, where it changes the pressure-driven route conductance, and downstream, where it can affect terminal delivery without necessarily preventing OA entry. This design was used because the published angiosome and choke-vessel literature implies several anatomically plausible but haemodynamically distinct mechanisms (Jajoria et al., 2020; Taylor et al., 2020).")
scenario_main = pd.DataFrame({
    "Scenario": [s[0] for s in scenario_specs],
    "Purpose": [
        "Reference CFD cohort",
        "Tests pressure reduction",
        "Tests pressure forcing",
        "Tests rheological resistance",
        "Tests topological interruption",
        "Tests permissive continuity",
        "Tests borderzone radius constriction",
        "Tests combined permissive high-risk conditions",
    ],
})
table_from_df(doc, scenario_main, ["Scenario", "Purpose"], ["Scenario", "Purpose"])
caption(doc, "Table 2", "Prespecified CFD scenario matrix.")

doc.add_heading("Results", level=1)
doc.add_heading("Pressure surplus and route conductance controlled OA entry", level=2)
para(doc, "The CFD simulations showed a nonlinear transition from no ophthalmic entry to route filling once pressure surplus and conductance were concordant. Because conductance scales with the fourth power of radius, small route-specific differences in diameter produced large differences in transported volume. Glabellar supratrochlear and nasal dorsal routes showed the highest OA entry proportions, consistent with clinical and anatomical reports that these territories provide short facial-to-ophthalmic connections in susceptible configurations (Beleznay et al., 2015; Taylor et al., 2020; Zheng et al., 2017). Nasolabial angular routes required higher pressure surplus or more favourable diameter to cross the same transport threshold.")
para(doc, "Mean pressure surplus was similar across the four route families because the inlet pressure boundary and ophthalmic back-pressure distribution were shared by design. Route anatomy therefore acted mainly through hydraulic resistance, flow and transported volume rather than through route-specific changes in mean pressure surplus.")
para(doc, "The transition maps also explain why anatomical danger zones cannot be interpreted as deterministic. Many profiles in high-risk routes did not reach OA entry because the effective radius, injection duration or intravascular fraction was insufficient. Conversely, some lower-ranked routes crossed the threshold when pressure surplus and conductance aligned. The CFD result is therefore not a simple ranking of cosmetic sites; it is a route-filling model in which surface location, vascular entry and pressure boundary must coexist.")
doc.add_picture(str(FIG / "figure2_pressure_geometry.png"), width=Inches(6.4))
caption(doc, "Figure 2", "Pressure-geometry CFD maps. Heat maps show the proportion of profiles reaching the OA transport threshold across pressure surplus and diameter bins.")
route_summary_display = table_display(
    route_summary,
    {
        "route": "Route",
        "n": "n",
        "OA_entry": "OA entry",
        "CRA_exposure": "CRA exposure",
        "pressure_surplus": "Pressure surplus (mmHg)",
        "flow_ml_s": "Flow (mL/s)",
        "transit_s": "Median transit (s)",
    },
)
table_from_df(doc, route_summary_display, route_summary_display.columns, route_summary_display.columns)
caption(doc, "Table 3", "Route-level CFD outputs in the baseline intravascular exposure cohort.")
cap_overall = (df["transit_s"] >= 120).mean()
cap_nlf = transit_cap_tbl.loc[transit_cap_tbl["route"] == "Nasolabial angular", "Profiles at 120 s cap"].iloc[0]
cap_forehead = transit_cap_tbl.loc[transit_cap_tbl["route"] == "Forehead supraorbital", "Profiles at 120 s cap"].iloc[0]
para(doc, f"Transit time was capped at 120 s to avoid assigning physiologically uninformative long residence times to non-transporting profiles. The cap was reached by {cap_overall:.1%} of baseline profiles, including {cap_nlf:.1%} of nasolabial angular profiles and {cap_forehead:.1%} of forehead supraorbital profiles, explaining the 120.0 s median values in the lower-conductance routes. The cap was applied only to the reported transit time; OA entry and CRA exposure were determined from transported volume and route dead volume before censoring, so capped profiles were not automatically counted as transporting.")

doc.add_heading("CRA exposure required terminal volume delivery beyond OA entry", level=2)
para(doc, "CRA exposure was not equivalent to OA entry. Many simulations filled the ophthalmic access route without transporting enough coherent material to satisfy the terminal CRA volume condition. This distinction is clinically important because ophthalmic entry, retinal branch occlusion and CRA exposure represent related but non-identical transport states. The model reproduced a steep transition in CRA fill fraction when flow rate, route dead volume and pressure duration aligned.")
para(doc, "This separation between OA entry and CRA exposure is clinically plausible. A filler bolus may reach the ophthalmic circulation but fragment, disperse or lodge before producing complete CRA filling. The CFD model represents this by requiring a second volume condition after OA entry. The outcome is a hierarchy of transport states: no route filling, partial ophthalmic access, complete ophthalmic access and terminal CRA exposure. This hierarchy is more physiologically interpretable than a single binary complication score.")
doc.add_picture(str(FIG / "figure3_route_outputs.png"), width=Inches(6.5))
caption(doc, "Figure 3", "Route-level OA entry, CRA exposure and transported fraction. The scatter plot links CFD flow rate with terminal CRA fill fraction.")
band_display = table_display(
    band_tbl,
    {
        "risk_band": "CRA band",
        "profiles": "Profiles",
        "pressure_surplus": "Pressure surplus (mmHg)",
        "diameter": "Diameter (mm)",
        "mu_eff": "Effective viscosity (Pa s)",
        "flow_ml_s": "Flow (mL/s)",
        "transported_ml": "Transported volume (mL)",
    },
)
table_from_df(doc, band_display, band_display.columns, band_display.columns)
caption(doc, "Table 4", "CRA exposure bands and associated CFD transport conditions.")

doc.add_heading("Filler rheology altered resistance through shear-thinning rather than a simple linear penalty", level=2)
para(doc, "The Carreau calculation converted high low-shear HA viscosity into substantially lower effective viscosity at high shear rates. Consequently, viscosity modified transport through resistance but did not act as an independent linear penalty. The effective viscosity distribution remained broad, reflecting route radius, injection pressure and intravascular bolus fraction. Low-flow profiles could remain near the low-shear state and therefore retain high apparent effective viscosity, whereas high CRA exposure bands tended to shear-thin to lower effective viscosity while carrying larger transported volumes. Transit times were shortest when transported volume exceeded route dead volume under high pressure surplus.")
para(doc, "The rheology result also cautions against interpreting product viscosity as a stand-alone protective property. High low-shear viscosity increased resistance, but the shear environment created by a narrow high-pressure route could lower effective viscosity substantially. Product rheology therefore interacts with injector force, route calibre and injection duration. This is the reason the manuscript reports viscosity as a CFD resistance modifier rather than as a universal clinical safety factor.")
doc.add_picture(str(FIG / "figure4_transport_physics.png"), width=Inches(6.5))
caption(doc, "Figure 4", "Transport physics outputs. Panels show inlet pressure, effective viscosity, transit time and the relation between transported volume and route dead volume.")

doc.add_heading("Counterfactual CFD stress tests separated pressure, topology and choke effects", level=2)
para(doc, "Scenario testing confirmed that pressure and topology were the principal determinants of route filling. Low-pressure boundary conditions reduced OA entry and CRA exposure. High-pressure forcing increased both outputs. Removing true anastomosis markedly reduced transport despite unchanged injection variables, whereas forced arterial continuity increased entry. High viscosity reduced flow, but its effect was moderated by shear-thinning.")
para(doc, "The magnitude of the topology effect was particularly important. When true anastomotic continuity was removed, the effective direct route radius became too small to transmit meaningful bolus volume in most profiles. This result is not a statistical association; it follows from the resistance equation. Radius reduction increases resistance nonlinearly, so a discontinuous or capillary-like pathway cannot be treated as equivalent to a true arterial anastomosis. This provides a physical basis for separating anatomical continuity from broad angiosome adjacency, a distinction that is implicit in the filler angiosome and choke-vessel literature but has not been quantified haemodynamically (Jajoria et al., 2020; Taylor et al., 2020).")
doc.add_picture(str(FIG / "figure5_stress_tests.png"), width=Inches(6.5))
caption(doc, "Figure 5", "Counterfactual CFD stress tests. The same 12000 baseline profiles were re-evaluated under each boundary, rheology, topology and resistance scenario.")
stress_display = stress_tbl[["Scenario", "Evaluated profiles", "OA entry proportion", "CRA exposure proportion", "Change OA vs baseline", "Change CRA vs baseline"]].copy()
stress_display["Scenario"] = ["S0 baseline", "S1 low pressure", "S2 high pressure", "S3 high viscosity", "S4 no anastomosis", "S5 forced continuity", "S6 high upstream choke", "S7 high-risk extremes"]
stress_display = table_display(
    stress_display,
    {
        "Scenario": "Condition",
        "Evaluated profiles": "n",
        "OA entry proportion": "OA entry",
        "CRA exposure proportion": "CRA exposure",
        "Change OA vs baseline": "Change OA",
        "Change CRA vs baseline": "Change CRA",
    },
)
table_from_df(doc, stress_display, stress_display.columns, stress_display.columns)
caption(doc, "Table 5", "Formal CFD stress-test outputs. Values are Monte Carlo proportions within the conditional intravascular exposure cohort.")

doc.add_heading("Choke resistance was protective only when it reduced upstream radius", level=2)
para(doc, "The choke-vessel hypothesis behaved differently depending on where resistance was applied. Table 6 is the choke-sensitivity grid and uses the same 12000 baseline profiles as Table 5. Its first row is the S0 baseline from Table 5, so the change columns are anchored to OA entry 0.227 and CRA exposure 0.189. When choke resistance reduced upstream radius, the fourth-power radius term produced a strong protective shift at high constriction. When the same resistance concept was restricted to the downstream terminal compartment, the OA entry proportion remained anchored to the S0 transport solution, while CRA exposure changed modestly because only the terminal volume condition was altered. The result does not dismiss choke anatomy. It operationalises the hypothesis raised by Jajoria et al. (2020) and Taylor et al. (2020) by asking whether the choke segment lies upstream of the relevant pressure-driven embolic pathway and whether its radius reduction is large enough to dominate pressure surplus.")
para(doc, "This resolves the circularity problem that affected the earlier score-based version of the work. In the CFD model, choke resistance is no longer assigned a small coefficient. It modifies radius or terminal volume, and the effect follows from the hydraulic equations. The model therefore permits choke resistance to become strongly protective when the anatomical assumptions support it. The conclusion is not that choke vessels are irrelevant, but that their protective capacity depends on anatomical position and degree of constriction.")
doc.add_picture(str(FIG / "figure6_choke_sensitivity.png"), width=Inches(6.4))
caption(doc, "Figure 6", "Choke-resistance sensitivity. Upstream radius reduction substantially altered OA entry at high constriction because hydraulic conductance scales with radius to the fourth power.")
choke_display = table_display(
    choke_sensitivity,
    {
        "Action site": "Action site",
        "Maximum radius reduction": "Max radius reduction",
        "OA entry proportion": "OA entry",
        "CRA exposure proportion": "CRA exposure",
        "Change OA vs baseline": "Change OA",
        "Change CRA vs baseline": "Change CRA",
    },
)
table_from_df(doc, choke_display, choke_display.columns, choke_display.columns)
caption(doc, "Table 6", "Choke sensitivity across upstream and downstream implementations. The S0 baseline row is identical to Table 5. Downstream-only resistance was applied after OA entry, so Change OA remains 0.000 across downstream rows; negative upstream Change OA values indicate protection.")

doc.add_heading("Discussion", level=1)
para(doc, "This reduced-order CFD study provides a physically grounded explanation of filler related visual loss. The model indicates that acute ophthalmic entry is not determined by injection site label alone. It emerges when a pressurised filler bolus enters a permissive arterial route with sufficient conductance and duration to exceed route dead volume. This is consistent with clinical reports of sudden visual loss after glabellar, nasal, nasolabial and forehead injections, and with angiographic descriptions of ophthalmic and retinal arterial occlusion after filler embolism (Beleznay et al., 2015; Carle et al., 2014; Kim et al., 2015; Zhao et al., 2024).")
para(doc, "The principal advance over a narrative anatomical model is that the equations expose the competing determinants of transport. Radius has a dominant effect because resistance scales with the inverse fourth power of radius. Length and branch angle modify the threshold. Filler viscosity contributes through shear-dependent resistance, not as a simple categorical product property. Pressure surplus is necessary but not sufficient; it must be interpreted in the context of route geometry and dead volume.")
para(doc, "For cosmetic science and aesthetic practice, this model shifts the discussion from danger-zone labelling to conditional haemodynamics. A glabellar or nasal route is dangerous not merely because of its name, but because short route length, larger arterial continuity and favourable branch geometry can reduce the pressure needed for bolus transport. A nasolabial or forehead route may be lower risk on average, yet still become dangerous under a permissive configuration. This is the kind of interaction that static anatomical diagrams cannot quantify, including diagrams based only on angiosome adjacency (Taylor et al., 2020; Zhu et al., 2018).")
para(doc, "The choke-vessel result is more nuanced than a simple rejection or endorsement. If choke resistance is placed upstream and creates major radius reduction, it is strongly protective in the CFD model. If it acts only downstream or at a terminal borderzone after ophthalmic entry has occurred, it cannot prevent OA entry and may only reduce CRA exposure. This distinction is important for aesthetic medicine because it converts a broad anatomical theory into a testable haemodynamic hypothesis (Jajoria et al., 2020; Taylor et al., 2020).")
para(doc, "The findings also help distinguish immediate and delayed clinical phenotypes. Immediate visual loss is most consistent with rapid arterial route filling and terminal embolic exposure. Delayed deterioration may involve venous congestion, secondary thrombosis, inflammatory swelling or redistribution of partially transported material. The present model focuses on acute arterial transport and should not be used to explain every delayed visual symptom. That narrower scope strengthens rather than weakens the model, because the equations are matched to the specific biological event under investigation.")
para(doc, "The model also clarifies why a universal safe pressure threshold is unlikely. The same plunger force can produce different intraluminal pressure depending on pressure coupling, vascular access, injection duration and local resistance. Similarly, the same pressure surplus can be harmless in a narrow or discontinuous route and dangerous in a large, short, low-angle anastomosis. Prevention therefore requires layered risk reduction rather than reliance on a single pressure rule.")
para(doc, "Several limitations remain. This is a reduced-order CFD network model, not a full three-dimensional patient-specific CFD reconstruction. Vessel walls were not deformable, pulsatility was simplified into pressure boundaries and the route geometry was sampled from published ranges rather than segmented from patient imaging. The model is conditional on intravascular exposure and cannot estimate clinical incidence. Nevertheless, the implemented equations are physically derived and auditable, and the scenario outputs are generated directly from the CFD transport model.")
para(doc, "A further limitation is that pressure coupling during inadvertent intravascular injection is difficult to measure directly. The model therefore sampled a plausible transmission range and stress tested pressure boundaries rather than claiming a single true coupling value. This uncertainty does not invalidate the model, but it identifies a key experimental priority. Bench systems that measure simultaneous plunger force, intraluminal pressure and filler displacement in small-calibre vessels would substantially improve future calibration.")
para(doc, "Future work should incorporate segmented facial vascular imaging, product-specific rheological measurements, pulsatile inlet waveforms and prospective adverse-event registries. The present model identifies the highest-value validation targets: route radius, true arterial continuity, pressure coupling during intravascular injection and the anatomical level of choke resistance.")

doc.add_heading("Conclusion", level=1)
para(doc, "A reduced-order CFD model of facial-to-ophthalmic filler transport showed that OA entry occurs when pressure surplus and route conductance deliver sufficient filler volume to exceed the anastomotic dead volume. CRA exposure requires additional terminal volume delivery. Choke resistance is protective when it reduces upstream radius, but distal resistance alone does not prevent OA entry. These findings support a pressure-driven embolic transport mechanism through permissive arterial continuity and provide a reproducible CFD framework for future patient-specific modelling.")
doc.add_heading("Data availability", level=1)
para(doc, "Machine-readable simulation outputs accompany this submission package: cfd_simulation_outputs.csv, cfd_route_summary.csv, cfd_scenario_evaluations.csv, cfd_stress_test_outputs.csv, cfd_choke_sensitivity.csv, cfd_choke_sensitivity_full_grid.csv and cfd_s7_route_outputs.csv. Before journal submission, these files should be deposited in Zenodo or an equivalent public repository and the permanent repository DOI added to this statement.")
doc.add_heading("Code availability", level=1)
para(doc, "The executable Python script used to generate the virtual cohort, solve the reduced-order CFD network, create figures and assemble the manuscript is supplied as cfd_reproducibility_script.py. Before journal submission, the final script should be deposited with the data release or in a linked GitHub repository archived by Zenodo, and the permanent repository DOI added to this statement.")
doc.add_heading("Funding", level=1)
para(doc, "The authors received no specific funding for this work.")
doc.add_heading("Conflicts of interest", level=1)
para(doc, "The authors declare no conflicts of interest.")
doc.add_heading("Author contributions", level=1)
para(doc, "Eqram Rahman: Conceptualization, Methodology, Software, Validation, Formal analysis, Visualization, Writing original draft, Writing review and editing.")
doc.add_heading("References", level=1)
for ref in references:
    p = doc.add_paragraph(ref)
    p.paragraph_format.left_indent = Inches(0.2)
    p.paragraph_format.first_line_indent = Inches(-0.2)
    p.paragraph_format.space_after = Pt(4)

main_path = OUT / "CFD_Filler_Visual_Loss_Manuscript.docx"
doc.save(main_path)

sup = Document()
setup_doc(sup)
p = sup.add_paragraph(style="Title")
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
p.add_run("Supplementary Methods: Reduced-order CFD Network Model")
para(sup, "This supplement provides the equations, implementation choices, route centres, scenario definitions, numerical checks and reproducibility details for the CFD network model.")
sup.add_heading("S1. Governing equations", level=1)
table_from_df(sup, nomenclature, ["Symbol", "Meaning", "SI unit"], ["Symbol", "Meaning", "SI unit"])
caption(sup, "Supplementary Table S1", "Nomenclature for the reduced-order CFD network model.")
for txt in [
    "Mass conservation was enforced at the reduced-order segment level by using the same volumetric flow Q through each solved route segment during a transport step.",
    "Pressure surplus was DeltaP = P_inlet - P_ophthalmic_back, with pressure solved in pascals and converted to mmHg only for reporting.",
    "Syringe pressure was P_syringe = F/A_syringe, where A_syringe = π(d_syringe/2)^2. Inlet pressure was P_inlet = MAP + k_coupling P_syringe after conversion to SI units.",
    "Ophthalmic back pressure was P_ophthalmic_back = MAP + 8 mmHg + 10 mmHg × tone, converted to pascals before subtraction from P_inlet.",
    "For each network edge, pressure-driven laminar flow was solved as Q = DeltaP/R.",
    "Hydraulic resistance was R = 8 μ L K_angle/(π r_eff^4), where μ is effective viscosity, L is route length, r_eff is effective radius and K_angle is the branch-angle loss multiplier.",
    "Branch loss was K_angle = 1 + 0.85(1 - cos theta).",
    "Carreau viscosity was eta = eta_inf + (eta0 - eta_inf)[1 + (lambda gamma_dot)^a]^((n - 1)/a).",
    "Route dead volume was V_dead = π r_eff² L. Bolus transport was V_transport = min(Q t_inj, V_intravascular). OA entry required V_transport >= 0.85 V_dead under positive pressure surplus. CRA exposure required V_transport >= V_dead + V_CRA.",
]:
    para(sup, txt)
sup.add_page_break()
sup.add_heading("S2. Route centres", level=1)
table_from_df(sup, route_centres, ["Route", "Diameter centre mm", "Length centre mm", "Angle centre degrees", "Continuity probability"], ["Route", "Diameter centre mm", "Length centre mm", "Angle centre degrees", "Continuity probability"])
caption(sup, "Supplementary Table S2", "Route-specific anatomical centres.")
sup.add_heading("S3. Sampling distributions", level=1)
sampling = pd.DataFrame([
    ["Random seed", str(SEED), "Reproducibility"],
    ["Profiles", str(N_BASE), "3000 per route"],
    ["Diameter", "Normal by route centre, SD 0.19, clipped 0.38-1.55 mm", "Anatomical"],
    ["Length", "Normal by route centre, SD 4.8, clipped 9-43 mm", "Anatomical"],
    ["Angle", "Normal by route centre, SD 11, clipped 5-82 degrees", "Geometry"],
    ["Force", "Uniform 6-30 N", "Pressure boundary"],
    ["Pressure coupling", "Uniform 0.025-0.14", "Intravascular transmission"],
    ["Bolus volume", "Uniform 0.015-0.18 ml", "Transport volume"],
    ["Injection duration", "Uniform 1.5-8.0 s", "Time boundary"],
    ["Low-shear HA viscosity", "Uniform 75-190 Pa.s", "Carreau input"],
    ["Tone", "Beta(2.2, 2.5)", "Back pressure and CRA volume"],
    ["Choke index", "Beta(1.8, 3.8)", "Resistance sensitivity"],
], columns=["Parameter", "Distribution", "Role"])
table_from_df(sup, sampling, ["Parameter", "Distribution", "Role"], ["Parameter", "Distribution", "Role"])
caption(sup, "Supplementary Table S3", "Sampling distributions.")
sup.add_heading("S4. Solver and scenario accounting", level=1)
para(sup, "The coupled flow-viscosity calculation used a fixed-point iteration scheme. Each profile started with low-shear HA viscosity, solved resistance and flow, updated shear rate and Carreau effective viscosity, and repeated through eight update cycles. The final update cycle means the eighth viscosity-flow pair used to compute transported volume, route filling, Reynolds number and wall shear. Relative flow and viscosity changes were audited at the final cycle; an adaptive tolerance-based solver should be used in future patient-specific implementations.")
scenario_accounting = pd.DataFrame(
    [
        ["Baseline profiles", N_BASE, "Four routes with 3000 profiles per route"],
        ["Scenario states", 8, "S0 baseline plus seven counterfactual re-evaluations"],
        ["Scenario rows", N_BASE * 8, "12000 baseline profiles multiplied by 8 scenario states"],
        ["Choke grid rows", len(choke_sensitivity_full), "Two action sites multiplied by five radius-reduction settings"],
    ],
    columns=["Item", "Count", "Explanation"],
)
table_from_df(sup, scenario_accounting, ["Item", "Count", "Explanation"], ["Item", "Count", "Explanation"])
caption(sup, "Supplementary Table S4", "Scenario accounting for the CFD model.")
sup.add_heading("S5. Scenario outputs", level=1)
table_from_df(sup, laminar_audit, ["Check", "Value", "Interpretation"], ["Check", "Value", "Interpretation"])
caption(sup, "Supplementary Table S5", "Laminar-flow assumption audit for the baseline CFD cohort.")
para(sup, "The median Reynolds number is very low because most baseline profiles occupied low-flow, high-viscosity states in which the transported bolus did not approach route-filling conditions. Higher Reynolds numbers occurred only in the extreme high-pressure, low-viscosity tail and did not determine the scenario ordering.")
sup.add_page_break()
table_from_df(sup, stress_tbl, ["Scenario", "Evaluated profiles", "OA entry proportion", "CRA exposure proportion", "Mean transported volume ml", "Mean flow ml/s", "Median transit s", "Change OA vs baseline"], ["Scenario", "n", "OA entry", "CRA exposure", "Transported ml", "Flow ml/s", "Transit s", "Change OA"])
caption(sup, "Supplementary Table S6", "Counterfactual CFD outputs.")
para(sup, "Supplementary Table S6 reports the same scenario proportions as Table 5 in the main manuscript; it is repeated here to keep the computational supplement self-contained.")
table_from_df(sup, s7_route, ["route", "Evaluated profiles", "OA entry proportion", "CRA exposure proportion", "Mean transported volume ml"], ["Route", "n", "OA entry", "CRA exposure", "Transported ml"])
caption(sup, "Supplementary Table S7", "S7 route-specific outputs.")
table_from_df(sup, choke_sensitivity_full, ["Action site", "Maximum radius reduction", "OA entry proportion", "CRA exposure proportion", "Change OA vs baseline", "Change CRA vs baseline"], ["Action site", "Max radius reduction", "OA entry", "CRA exposure", "Change OA", "Change CRA"])
caption(sup, "Supplementary Table S8", "Full choke-resistance sensitivity grid. The zero-reduction rows are anchored to the same S0 baseline used in Table 5 and Table 6. Downstream-only resistance was applied after OA entry, so downstream Change OA remains 0.000.")
sup.add_heading("S6. Reproducibility pseudocode", level=1)
para(sup, "Set seed to 42. Generate 12000 profiles equally across four routes. For each profile, sample geometry, force, pressure coupling, volume, duration, viscosity, blood pressure, tone, continuity and choke index. Convert force to syringe pressure using plunger area. Apply pressure coupling to define inlet pressure. Compute ophthalmic back pressure from mean arterial pressure and tone. Iterate viscosity, resistance, flow and shear rate through eight fixed-point update cycles. Compute transported volume, route dead volume, OA entry, CRA exposure, transit time, Reynolds number and wall shear from the final update cycle. Repeat the same profiles under S0-S7 scenarios. For choke sensitivity, repeat the baseline profiles across upstream and downstream action sites and radius-reduction settings of 0, 0.15, 0.35, 0.60 and 0.90.")
sup.add_heading("Supplementary References", level=1)
for ref in references:
    p = sup.add_paragraph(ref)
    p.paragraph_format.left_indent = Inches(0.2)
    p.paragraph_format.first_line_indent = Inches(-0.2)
    p.paragraph_format.space_after = Pt(4)
supp_path = OUT / "CFD_Filler_Visual_Loss_Supplement.docx"
sup.save(supp_path)

df.to_csv(OUT / "cfd_simulation_outputs.csv", index=False)
route_summary.to_csv(OUT / "cfd_route_summary.csv", index=False)
scenario_df.to_csv(OUT / "cfd_scenario_evaluations.csv", index=False)
stress_tbl.to_csv(OUT / "cfd_stress_test_outputs.csv", index=False)
choke_sensitivity.to_csv(OUT / "cfd_choke_sensitivity.csv", index=False)
choke_sensitivity_full.to_csv(OUT / "cfd_choke_sensitivity_full_grid.csv", index=False)
s7_route.to_csv(OUT / "cfd_s7_route_outputs.csv", index=False)
shutil.copyfile(Path(__file__), OUT / "cfd_reproducibility_script.py")

print(main_path)
print(supp_path)
