"""
preprocess.py
-------------
Build a longitudinal panel dataset from Synthea CSV for msm analysis.

Output: data_processed/ckd_panel.csv
Columns:
    patient_id    - unique patient identifier
    age           - patient age in years at observation time
    state         - CKD state (1-5 by eGFR, 6 = death)
    sex           - 1 = male, 0 = female  (time-invariant)
    htn_baseline  - 1 if hypertension present at first CKD observation (time-invariant)

State mapping (eGFR thresholds):
    1: Stage 1  eGFR >= 90
    2: Stage 2  60 <= eGFR < 90
    3: Stage 3  30 <= eGFR < 60
    4: Stage 4  15 <= eGFR < 30
    5: Stage 5  eGFR < 15
    6: Death    (absorbing)
"""

import pandas as pd
import numpy as np
import os

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA = os.path.join(os.path.dirname(__file__), "..", "output", "csv")
OUT  = os.path.join(os.path.dirname(__file__), "..", "data_processed")
os.makedirs(OUT, exist_ok=True)

# ── eGFR → state ──────────────────────────────────────────────────────────────
def egfr_to_state(v):
    if v >= 90: return 1
    if v >= 60: return 2
    if v >= 30: return 3
    if v >= 15: return 4
    return 5

# ─────────────────────────────────────────────────────────────────────────────
# 1. Identify CKD patients
# ─────────────────────────────────────────────────────────────────────────────
print("Loading conditions...")
conditions = pd.read_csv(
    f"{DATA}/conditions.csv",
    usecols=["PATIENT", "START", "DESCRIPTION"],
)
ckd_cond = conditions[
    conditions["DESCRIPTION"].str.contains("Chronic kidney disease", case=False, na=False)
]
ckd_patients = set(ckd_cond["PATIENT"].unique())
print(f"  CKD patients: {len(ckd_patients):,}")

# ─────────────────────────────────────────────────────────────────────────────
# 2. Patient demographics
# ─────────────────────────────────────────────────────────────────────────────
print("Loading patients...")
patients = pd.read_csv(
    f"{DATA}/patients.csv",
    usecols=["Id", "BIRTHDATE", "DEATHDATE", "GENDER"],
)
patients = patients[patients["Id"].isin(ckd_patients)].copy()
patients["BIRTHDATE"] = pd.to_datetime(patients["BIRTHDATE"]).dt.tz_localize(None)
patients["DEATHDATE"] = pd.to_datetime(patients["DEATHDATE"]).dt.tz_localize(None)
patients["sex"] = (patients["GENDER"] == "M").astype(int)
patients = patients.set_index("Id")

# ─────────────────────────────────────────────────────────────────────────────
# 3. eGFR observations → base panel (one row per eGFR measurement)
# ─────────────────────────────────────────────────────────────────────────────
print("Loading eGFR observations (chunked)...")
egfr_parts = []
for chunk in pd.read_csv(
    f"{DATA}/observations.csv",
    usecols=["DATE", "PATIENT", "CODE", "VALUE"],
    chunksize=500_000,
):
    e = chunk[
        chunk["CODE"].isin(["33914-3", "62238-1"])
        & chunk["PATIENT"].isin(ckd_patients)
    ]
    if len(e):
        egfr_parts.append(e)

egfr = pd.concat(egfr_parts, ignore_index=True)
egfr["VALUE"] = pd.to_numeric(egfr["VALUE"], errors="coerce")
egfr = egfr.dropna(subset=["VALUE"])
egfr["DATE"] = pd.to_datetime(egfr["DATE"]).dt.tz_localize(None)

# Keep one measurement per patient per date
egfr = (
    egfr.sort_values(["PATIENT", "DATE", "VALUE"])
    .drop_duplicates(subset=["PATIENT", "DATE"], keep="first")
)
egfr["state"] = egfr["VALUE"].apply(egfr_to_state)

# Merge birthdate and sex
egfr = egfr.merge(
    patients[["BIRTHDATE", "sex"]],
    left_on="PATIENT", right_index=True,
)
egfr["age"] = (egfr["DATE"] - egfr["BIRTHDATE"]).dt.days / 365.25

# ─────────────────────────────────────────────────────────────────────────────
# 4. Add death rows (state = 6) for deceased patients
# ─────────────────────────────────────────────────────────────────────────────
print("Adding death rows...")
dead = patients[patients["DEATHDATE"].notna()].copy()
death_rows = pd.DataFrame({
    "PATIENT":   dead.index,
    "DATE":      dead["DEATHDATE"].values,
    "VALUE":     np.nan,
    "state":     6,
    "BIRTHDATE": dead["BIRTHDATE"].values,
    "sex":       dead["sex"].values,
})
death_rows["age"] = (death_rows["DATE"] - death_rows["BIRTHDATE"]).dt.days / 365.25

panel = pd.concat(
    [egfr[["PATIENT", "DATE", "state", "sex", "age"]],
     death_rows[["PATIENT", "DATE", "state", "sex", "age"]]],
    ignore_index=True,
)
panel["DATE"] = pd.to_datetime(panel["DATE"]).dt.tz_localize(None)
panel = panel.sort_values(["PATIENT", "age"]).reset_index(drop=True)

# ─────────────────────────────────────────────────────────────────────────────
# 5. Hypertension onset → baseline indicator
#    htn_baseline = 1 if diagnosed with hypertension at or before first
#    eGFR observation (i.e., at the time of CKD detection)
# ─────────────────────────────────────────────────────────────────────────────
print("Processing hypertension onset...")

htn_rows = conditions[
    conditions["DESCRIPTION"].str.contains("Hypertension", case=False, na=False)
    & conditions["PATIENT"].isin(ckd_patients)
][["PATIENT", "START"]].copy()
htn_rows["START"] = pd.to_datetime(htn_rows["START"]).dt.tz_localize(None)
htn_onset = htn_rows.groupby("PATIENT")["START"].min()

# First eGFR observation date per patient
first_egfr_date = (
    egfr.sort_values(["PATIENT", "DATE"])
    .drop_duplicates(subset="PATIENT", keep="first")
    .set_index("PATIENT")["DATE"]
)

# htn_baseline = 1 if HTN diagnosed on or before first eGFR observation
htn_baseline_map = {}
for pid in ckd_patients:
    if pid not in first_egfr_date.index:
        continue
    first_date = first_egfr_date[pid]
    onset = htn_onset.get(pid, pd.NaT)
    if pd.isna(onset):
        htn_baseline_map[pid] = 0
    else:
        htn_baseline_map[pid] = int(onset <= first_date)

panel["htn_baseline"] = panel["PATIENT"].map(htn_baseline_map).fillna(0).astype(int)

# ─────────────────────────────────────────────────────────────────────────────
# 6. Final cleanup
# ─────────────────────────────────────────────────────────────────────────────
print("Finalizing dataset...")
panel_final = (
    panel[["PATIENT", "age", "state", "sex", "htn_baseline"]]
    .rename(columns={"PATIENT": "patient_id"})
    .dropna(subset=["age", "state"])
    .sort_values(["patient_id", "age"])
    .reset_index(drop=True)
)

# Deduplicate same (patient, age): keep max state (death row beats eGFR row)
panel_final = (
    panel_final
    .sort_values(["patient_id", "age", "state"], ascending=[True, True, False])
    .drop_duplicates(subset=["patient_id", "age"], keep="first")
    .reset_index(drop=True)
)

# Enforce monotone non-decreasing state per patient (CKD is progressive)
panel_final["state"] = (
    panel_final.groupby("patient_id")["state"].transform("cummax")
)

# Thin within each 6-month window:
#   - No state change in window  → keep only the first record
#   - State changes within window → keep last record before each transition
#                                    AND first record after each transition
panel_final["age_bin"] = (panel_final["age"] // 0.5).astype(int)
panel_final = panel_final.sort_values(["patient_id", "age"]).reset_index(drop=True)

panel_final["_next_state"] = panel_final.groupby(["patient_id", "age_bin"])["state"].shift(-1)
panel_final["_prev_state"] = panel_final.groupby(["patient_id", "age_bin"])["state"].shift(1)
panel_final["_bin_has_change"] = (
    panel_final.groupby(["patient_id", "age_bin"])["state"].transform("nunique") > 1
)

_is_first_in_bin      = ~panel_final.duplicated(subset=["patient_id", "age_bin"], keep="first")
_is_last_before_change = (
    panel_final["_next_state"].notna() &
    (panel_final["state"] != panel_final["_next_state"])
)
_is_first_after_change = (
    panel_final["_prev_state"].notna() &
    (panel_final["state"] != panel_final["_prev_state"])
)

_keep = (
    _is_last_before_change |
    _is_first_after_change |
    (_is_first_in_bin & ~panel_final["_bin_has_change"])
)
panel_final = (
    panel_final[_keep]
    .drop(columns=["age_bin", "_next_state", "_prev_state", "_bin_has_change"])
    .reset_index(drop=True)
)

# Drop patients with only 1 row (no transitions observable)
counts = panel_final.groupby("patient_id").size()
panel_final = panel_final[
    panel_final["patient_id"].isin(counts[counts >= 2].index)
]

# ─────────────────────────────────────────────────────────────────────────────
# 7. Summary and save
# ─────────────────────────────────────────────────────────────────────────────
pts_first = (
    panel_final.sort_values(["patient_id", "age"])
    .drop_duplicates(subset="patient_id", keep="first")
)

n_pts    = panel_final["patient_id"].nunique()
n_rows   = len(panel_final)
n_male   = int(pts_first["sex"].sum())
n_female = n_pts - n_male
n_htn    = int(pts_first["htn_baseline"].sum())
n_htn_m  = int(pts_first.loc[pts_first["sex"]==1, "htn_baseline"].sum())
n_htn_f  = int(pts_first.loc[pts_first["sex"]==0, "htn_baseline"].sum())

print(f"\n{'='*55}")
print(f"Final panel: {n_rows:,} rows, {n_pts:,} patients")
print(f"\n--- Cohort Summary ---")
print(f"  Male:   {n_male:,}  ({100*n_male/n_pts:.1f}%)")
print(f"  Female: {n_female:,}  ({100*n_female/n_pts:.1f}%)")
print(f"\n  Hypertension at CKD diagnosis: {n_htn:,}  ({100*n_htn/n_pts:.1f}%)")
print(f"    Male:   {n_htn_m:,}  ({100*n_htn_m/n_male:.1f}% of males)")
print(f"    Female: {n_htn_f:,}  ({100*n_htn_f/n_female:.1f}% of females)")

print(f"\nState distribution (all rows):")
print(panel_final["state"].value_counts().sort_index().to_string())

out_path = os.path.join(OUT, "ckd_panel.csv")
panel_final.to_csv(out_path, index=False)
print(f"\nPanel saved → {out_path}")
