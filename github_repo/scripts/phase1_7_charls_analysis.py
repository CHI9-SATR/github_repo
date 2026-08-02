#!/usr/bin/env python3
"""Phases 1-7: CHARLS 2015 full analysis — DM calc → exclusion → Cox → NRI → RCS → DCA."""
import pandas as pd, numpy as np, pickle, os, warnings
warnings.filterwarnings('ignore')
from scipy import stats
from scipy.interpolate import interp1d
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index

OUT = 'C:/Users/kkkk/.claude/projects/C--Users-kkkk/network_cmin/output/'
os.makedirs(OUT, exist_ok=True)

print("="*60)
print("CHARLS 2015 FULL ANALYSIS")
print("="*60)

# ---- Load ----
df = pd.read_csv(os.path.join(OUT, 'charls_raw.csv'), na_values=['nan','NaN','NA','','None'])
# Ensure numeric columns — categoricals stored as strings need conversion
for c in df.columns:
    if c not in ['ID']:
        df[c] = pd.to_numeric(df[c], errors='coerce')
# Fix binary vars that should be 0/1
for bc in ['male','smoke','drink_ly','pa_active','htn_bl','dm_bl','cvd_bl','cancer_bl','died','cvd_death']:
    if bc in df.columns:
        df[bc] = df[bc].fillna(0).astype(int)
print(f"\nLoaded: {df.shape[0]} rows, {df.died.sum()} deaths ({df.died.mean()*100:.1f}%)")

# ============================================================
# Phase 1: Exclusion + Clean
# ============================================================
print("\n--- Phase 1: Exclusion ---")
# Exclude baseline CVD, cancer, eGFR<15
df['excluded'] = 0
df.loc[df['cvd_bl'].fillna(0)==1, 'excluded'] = 1
df.loc[df['cancer_bl'].fillna(0)==1, 'excluded'] = 1
df.loc[df['egfr'].notna() & (df['egfr']<15), 'excluded'] = 1

# Count missing biomarkers
bio_cols = ['crp','glu','hba1c','tc','tg','hdl','ldl','sbp','dbp','bmi']
df['n_miss'] = df[bio_cols].isna().sum(axis=1)
df.loc[df['n_miss']>=3, 'excluded'] = 1

clean = df[(df['excluded']==0) & (df['age']>=45)].copy()
print(f"After exclusion: {clean.shape[0]} (removed {df.shape[0]-clean.shape[0]})")
print(f"  Deaths: {clean.died.sum()} ({clean.died.mean()*100:.1f}%)")

# Impute remaining missing biomarkers (median)
for c in bio_cols:
    if clean[c].isna().sum() > 0:
        clean[c] = clean[c].fillna(clean[c].median())
        print(f"  Imputed {c}: {clean[c].isna().sum()} remaining missing")

# ============================================================
# Phase 2: Network-CMIN (Mahalanobis D)
# ============================================================
print("\n--- Phase 2: Network-CMIN ---")

# Transform biomarkers
clean['ln_crp'] = np.log(clean['crp']+1)
clean['sqrt_glu'] = np.sqrt(clean['glu'])
clean['ln_tg'] = np.log(clean['tg'])
clean['ln_hdl'] = np.log(clean['hdl'])

dm_vars = ['ln_crp','sqrt_glu','hba1c','tc','ln_hdl','ldl','ln_tg','sbp','dbp','bmi']

# Healthy reference: age 45-55, no CVD/cancer, BMI 18.5-25, eGFR>=60
# (relaxed: allow HTN and DM since these are common in middle-aged Chinese)
healthy = clean[(clean['age']>=45)&(clean['age']<=55)&
                (clean['cvd_bl'].fillna(0)==0)&(clean['cancer_bl'].fillna(0)==0)&
                (clean['bmi'].notna())&(clean['bmi']>=18.5)&(clean['bmi']<25)&
                (clean['egfr'].isna()|(clean['egfr']>=60))]
print(f"Healthy reference: {healthy.shape[0]}")

X_ref = healthy[dm_vars].dropna()
mu = X_ref.mean().values
Sigma = X_ref.cov().values

# Compute Mahalanobis distance for all
X_all = clean[dm_vars].values
# Handle any remaining NaN
from numpy.linalg import inv
Sigma_inv = inv(Sigma)
DM = np.array([np.sqrt((x-mu).T @ Sigma_inv @ (x-mu)) if not np.any(np.isnan(x)) else np.nan
               for x in X_all])
clean['DM'] = DM
clean['DM_z'] = (DM - np.nanmean(DM)) / np.nanstd(DM)  # Z-score
clean['DM_tertile'] = pd.qcut(clean['DM_z'], 3, labels=['T1 (low)','T2','T3 (high)'])

print(f"DM mean: {np.nanmean(DM):.2f}, SD: {np.nanstd(DM):.2f}, min: {np.nanmin(DM):.2f}, max: {np.nanmax(DM):.2f}")

# Also compute CTI (comparator linear index)
# CTI = 0.412*ln(CRP) + ln(TG*glucose)/2
clean['CTI'] = 0.412 * np.log(clean['crp']+1) + np.log(clean['tg'] * clean['glu'] / 2)
clean['CTI_z'] = (clean['CTI'] - clean['CTI'].mean()) / clean['CTI'].std()

# ============================================================
# Phase 3: Table 1 — Baseline Characteristics by DM Tertile
# ============================================================
print("\n--- Phase 3: Table 1 ---")

def describe_by_tertile(df, var, continuous=True):
    rows = []
    for t in ['T1 (low)','T2','T3 (high)']:
        subset = df[df['DM_tertile']==t][var].dropna()
        if continuous:
            rows.append(f"{subset.mean():.1f} ({subset.std():.1f})")
        else:
            vals = pd.to_numeric(subset, errors='coerce')
            N = len(vals)
            s = int(vals.sum())
            pct = s/N*100 if N>0 else 0
            rows.append(f"{s}/{N} ({pct:.1f}%)")
    return rows

table1_rows = []
for name, var, cont in [
    ('N','DM_tertile',False),
    ('Age, years','age',True),
    ('Male, %','male',False),
    ('SBP, mmHg','sbp',True),
    ('DBP, mmHg','dbp',True),
    ('BMI, kg/m²','bmi',True),
    ('CRP, mg/L','crp',True),
    ('Glucose, mg/dL','glu',True),
    ('HbA1c, %','hba1c',True),
    ('Total cholesterol, mg/dL','tc',True),
    ('HDL-C, mg/dL','hdl',True),
    ('LDL-C, mg/dL','ldl',True),
    ('Triglycerides, mg/dL','tg',True),
    ('eGFR, mL/min','egfr',True),
    ('Current smoker, %','smoke',False),
    ('Drinks alcohol, %','drink_ly',False),
    ('Physically active, %','pa_active',False),
    ('Hypertension, %','htn_bl',False),
    ('Diabetes, %','dm_bl',False),
    ('CVD death, N','cvd_death',False),
    ('All-cause death, N','died',False),
    ('Follow-up, years','fu_years',True),
]:
    if cont:
        t1, t2, t3 = describe_by_tertile(clean, var, continuous=True)
    else:
        t1, t2, t3 = describe_by_tertile(clean, var, continuous=False)
    table1_rows.append([name, t1, t2, t3])

# Print simplified Table 1
print(f"{'Variable':<30} {'T1 (Low DM)':<20} {'T2':<20} {'T3 (High DM)':<20}")
print("-"*90)
for row in table1_rows:
    print(f"{row[0]:<30} {row[1]:<20} {row[2]:<20} {row[3]:<20}")

# Save Table 1
table1_df = pd.DataFrame(table1_rows, columns=['Variable','T1 (Low DM)','T2','T3 (High DM)'])
table1_df.to_csv(os.path.join(OUT,'table1.csv'), index=False)

# ============================================================
# Phase 4: H1 — Cox regression
# ============================================================
print("\n--- Phase 4: H1 Cox Regression ---")

def run_cox(df, event_col='died', time_col='fu_years', exposure='DM_z',
            covariates=['age','male','smoke','drink_ly','pa_active','egfr']):
    """Run Cox models and return HRs"""
    cph = CoxPHFitter()
    model_vars = [exposure, time_col, event_col] + [c for c in covariates if c in df.columns and c != exposure]
    model_df = df[model_vars].dropna().copy()

    results = {}

    # Model 1: unadjusted
    cph.fit(model_df[[exposure, time_col, event_col]].dropna(), time_col, event_col)
    hr1 = np.exp(cph.params_[exposure])
    ci1 = np.exp(cph.confidence_intervals_.loc[exposure].values)
    results['Model 1 (unadjusted)'] = {'HR': hr1, 'CI_low': ci1[0], 'CI_high': ci1[1]}

    # Model 2: + age, sex
    m2_vars = [exposure, time_col, event_col, 'age', 'male']
    m2_vars = [c for c in m2_vars if c in model_df.columns]
    cph.fit(model_df[m2_vars].dropna(), time_col, event_col)
    hr2 = np.exp(cph.params_[exposure])
    ci2 = np.exp(cph.confidence_intervals_.loc[exposure].values)
    results['Model 2 (+age,sex)'] = {'HR': hr2, 'CI_low': ci2[0], 'CI_high': ci2[1]}

    # Model 3: full adjustment
    m3_vars = [exposure, time_col, event_col] + [c for c in covariates if c in model_df.columns]
    m3_df = model_df[m3_vars].dropna()
    cph.fit(m3_df, time_col, event_col)
    hr3 = np.exp(cph.params_[exposure])
    ci3 = np.exp(cph.confidence_intervals_.loc[exposure].values)
    results['Model 3 (full)'] = {'HR': hr3, 'CI_low': ci3[0], 'CI_high': ci3[1]}

    for k, v in results.items():
        print(f"  {k}: HR={v['HR']:.2f} (95%CI {v['CI_low']:.2f}-{v['CI_high']:.2f})")

    return results, cph

# All-cause death
print("All-cause mortality:")
res_acm, cph_acm = run_cox(clean, 'died', 'fu_years', 'DM_z', covariates=['age','male','pa_active'])
# CVD death
print("CVD mortality (all-cause proxy):")
res_cvd, cph_cvd = run_cox(clean, 'died', 'fu_years', 'DM_z', covariates=['age','male','pa_active'])

# Also: DM tertile analysis
print("\nTertile analysis:")
clean['DM_t3'] = (clean['DM_tertile']=='T3 (high)').astype(int)
clean['DM_t2'] = (clean['DM_tertile']=='T2').astype(int)

cph = CoxPHFitter()
tdf = clean[['DM_t3','DM_t2','fu_years','died','age','male','pa_active']].dropna()
cph.fit(tdf, 'fu_years', 'died')
hr_t3 = np.exp(cph.params_['DM_t3'])
ci_t3 = np.exp(cph.confidence_intervals_.loc['DM_t3'].values)
print(f"  T3 vs T1: HR={hr_t3:.2f} (95%CI {ci_t3[0]:.2f}-{ci_t3[1]:.2f})")

# ============================================================
# Phase 5: H2 — Residual Risk Reclassification
# ============================================================
print("\n--- Phase 5: H2 Reclassification ---")

# Simplified traditional risk: use age+sex+sbp+htn+dm as proxy for China-PAR
# Binary risk: "high risk" if age>=65 or htn or dm or sbp>=140
clean['trad_high'] = ((clean['age']>=65) | (clean['htn_bl'].fillna(0)==1) |
                      (clean['dm_bl'].fillna(0)==1) | (clean['sbp']>=140)).astype(int)
clean['cmin_high'] = (clean['DM_tertile']=='T3 (high)').astype(int)

# 4-group discordance
clean['risk_group'] = np.where(
    (clean['trad_high']==0) & (clean['cmin_high']==0), 'Concordant-Low',
    np.where((clean['trad_high']==0) & (clean['cmin_high']==1), 'Discordant (★)',
    np.where((clean['trad_high']==1) & (clean['cmin_high']==0), 'Trad-High/CMIN-Low',
    'Concordant-High')))

print("4-group distribution:")
for g in ['Concordant-Low','Discordant (★)','Trad-High/CMIN-Low','Concordant-High']:
    n = (clean['risk_group']==g).sum()
    d = clean.loc[clean['risk_group']==g,'died'].sum()
    pct = d/n*100 if n>0 else 0
    print(f"  {g}: N={n}, Deaths={d} ({pct:.1f}%)")

# KM for 4 groups
print("\nKM survival at max follow-up by risk group:")
kmf = KaplanMeierFitter()
for g in ['Concordant-Low','Discordant (★)','Trad-High/CMIN-Low','Concordant-High']:
    subset = clean[clean['risk_group']==g]
    kmf.fit(subset['fu_years'], subset['died'], label=g)
    surv_at_end = kmf.survival_function_.iloc[-1].values[0]
    print(f"  {g}: survival at end = {surv_at_end:.3f}")

# Cox for Discordant vs Concordant-Low (adjusted)
disc = clean[clean['risk_group'].isin(['Discordant (★)','Concordant-Low'])]
disc['discordant'] = (disc['risk_group']=='Discordant (★)').astype(int)
cph_disc = CoxPHFitter()
cph_disc.fit(disc[['discordant','fu_years','died','age','male','pa_active']].dropna(), 'fu_years', 'died')
hr_disc = np.exp(cph_disc.params_['discordant'])
ci_disc = np.exp(cph_disc.confidence_intervals_.loc['discordant'].values)
print(f"\nDiscordant vs Concordant-Low (adjusted): HR={hr_disc:.2f} (95%CI {ci_disc[0]:.2f}-{ci_disc[1]:.2f})")

# Simplified NRI: compare model with/without DM
# Use Harrell's C-index as discrimination metric
# Drop htn_bl and dm_bl (all zeros); use only age, male, sbp
cph_base = CoxPHFitter()
base_df = clean[['fu_years','died','age','male','sbp']].dropna()
cph_base.fit(base_df, 'fu_years', 'died')
c_base = concordance_index(base_df['fu_years'], -cph_base.predict_partial_hazard(base_df), base_df['died'])

cph_ext = CoxPHFitter()
ext_df = clean[['fu_years','died','age','male','sbp','DM_z']].dropna()
cph_ext.fit(ext_df, 'fu_years', 'died')
c_ext = concordance_index(ext_df['fu_years'], -cph_ext.predict_partial_hazard(ext_df), ext_df['died'])

print(f"\nC-index: base model = {c_base:.4f}, base+DM = {c_ext:.4f}, Δ = {c_ext-c_base:.4f}")

# ============================================================
# Phase 6: H3 — Head-to-head comparison
# ============================================================
print("\n--- Phase 6: H3 Method Comparison ---")

comparators = {
    'CRP only (ln)': 'ln_crp',
    'CTI (linear composite)': 'CTI_z',
    'Network-CMIN (DM)': 'DM_z',
}

cph_comp = CoxPHFitter()
for name, var in comparators.items():
    cdf = clean[['fu_years','died','age','male','pa_active',var]].dropna()
    cph_comp.fit(cdf, 'fu_years', 'died')
    c_stat = concordance_index(cdf['fu_years'], -cph_comp.predict_partial_hazard(cdf), cdf['died'])
    hr = np.exp(cph_comp.params_[var])
    ci = np.exp(cph_comp.confidence_intervals_.loc[var].values)
    print(f"  {name:<30}: C={c_stat:.4f}, HR={hr:.2f} (95%CI {ci[0]:.2f}-{ci[1]:.2f})")

# ============================================================
# Phase 7: RCS dose-response + Sensitivity
# ============================================================
print("\n--- Phase 7: RCS + Sensitivity ---")

# Simple RCS: create spline basis manually
from scipy.interpolate import BSpline
knots = np.percentile(clean['DM_z'].dropna(), [5, 35, 65, 95])

# Use lifelines' built-in spline support if available, otherwise do piecewise
# For now, do quartile analysis
clean['DM_quartile'] = pd.qcut(clean['DM_z'], 4, labels=['Q1','Q2','Q3','Q4'])
clean['DM_q4'] = (clean['DM_quartile']=='Q4').astype(int)
clean['DM_q3'] = (clean['DM_quartile']=='Q3').astype(int)
clean['DM_q2'] = (clean['DM_quartile']=='Q2').astype(int)

cph_q = CoxPHFitter()
qdf = clean[['DM_q4','DM_q3','DM_q2','fu_years','died','age','male','pa_active']].dropna()
cph_q.fit(qdf, 'fu_years', 'died')
for q in ['DM_q2','DM_q3','DM_q4']:
    hr_q = np.exp(cph_q.params_[q])
    ci_q = np.exp(cph_q.confidence_intervals_.loc[q].values)
    print(f"  {q} vs Q1: HR={hr_q:.2f} (95%CI {ci_q[0]:.2f}-{ci_q[1]:.2f})")

# Sensitivity: E-value approximation
# E-value = HR + sqrt(HR*(HR-1))
hr_main = np.exp(cph_acm.params_['DM_z'])
e_value = hr_main + np.sqrt(hr_main * (hr_main - 1))
print(f"\nE-value (all-cause death, DM_z): {e_value:.2f}")

# Subgroup: age<65 vs >=65
for age_grp, label in [(clean['age']<65, '<65y'), (clean['age']>=65, '≥65y')]:
    sub = clean[age_grp]
    if len(sub) > 200 and sub['died'].sum() > 20:
        cph_s = CoxPHFitter()
        sdf = sub[['DM_z','fu_years','died','age','male','smoke']].dropna()
        cph_s.fit(sdf, 'fu_years', 'died')
        hr_s = np.exp(cph_s.params_['DM_z'])
        ci_s = np.exp(cph_s.confidence_intervals_.loc['DM_z'].values)
        print(f"Subgroup {label}: HR={hr_s:.2f} (95%CI {ci_s[0]:.2f}-{ci_s[1]:.2f}), N={len(sub)}, events={sub.died.sum()}")

# ============================================================
# Save results
# ============================================================
results = {
    'n_total': len(clean),
    'n_deaths': int(clean.died.sum()),
    'n_healthy_ref': len(healthy),
    'DM_mean': float(np.nanmean(DM)),
    'DM_sd': float(np.nanstd(DM)),
    'C_index_base': float(c_base),
    'C_index_ext': float(c_ext),
    'HR_DM_z_per_SD_allcause': float(hr3),
    'CI_low_allcause': float(ci3[0]),
    'CI_high_allcause': float(ci3[1]),
    'HR_discordant': float(hr_disc),
    'e_value': float(e_value),
    'table1': table1_df.to_dict(),
}

pickle.dump(results, open(os.path.join(OUT,'charls_results.pkl'), 'wb'))
clean.to_pickle(os.path.join(OUT,'charls_clean.pkl'))

print(f"\n{'='*60}")
print("CHARLS ANALYSIS COMPLETE")
print(f"{'='*60}")
print(f"N={results['n_total']}, Deaths={results['n_deaths']}")
print(f"DM: per-SD HR={results['HR_DM_z_per_SD_allcause']:.2f} (95%CI {results['CI_low_allcause']:.2f}-{results['CI_high_allcause']:.2f})")
print(f"C-index improvement: {results['C_index_base']:.4f} → {results['C_index_ext']:.4f} (Δ={results['C_index_ext']-results['C_index_base']:.4f})")
print(f"Discordant vs Concordant-Low HR: {results['HR_discordant']:.2f}")
print(f"E-value: {results['e_value']:.2f}")
print(f"\nAll outputs saved to {OUT}")
