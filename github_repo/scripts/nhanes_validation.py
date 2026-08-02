#!/usr/bin/env python3
"""NHANES external validation: DM with available biomarkers → mortality."""
import pandas as pd, numpy as np, pickle, os, warnings
warnings.filterwarnings('ignore')
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index

OUT = 'C:/Users/kkkk/.claude/projects/C--Users-kkkk/network_cmin/output/'

print("="*60)
print("NHANES EXTERNAL VALIDATION")
print("="*60)

# Load
nh = pd.read_csv(os.path.join(OUT, 'nhanes_raw.csv'))
print(f"Loaded: {nh.shape[0]} rows")

# ---- Clean & Exclude ----
nh['eligible'] = nh['eligible'].fillna(0).astype(int)
nhanes = nh[nh['eligible']==1].copy()

# Exclude baseline CVD, cancer, age<45
nhanes['cvd_bl'] = nhanes['cvd_bl'].fillna(0).astype(int)
nhanes['cancer_bl'] = nhanes['cancer_bl'].fillna(0).astype(int)
nhanes = nhanes[(nhanes['age']>=45) & (nhanes['cvd_bl']==0) & (nhanes['cancer_bl']==0)]

print(f"After exclusion (age>=45, no CVD/cancer): {nhanes.shape[0]}, deaths={int(nhanes['died'].sum())}, CVD deaths={int(nhanes['cvd_death'].sum())}")

# ---- Identify available biomarkers ----
bio_avail = {}
for b in ['crp','glu','hba1c','tc','tg','hdl','ldl','creat','sbp','dbp','bmi']:
    if b in nhanes.columns:
        avail = nhanes[b].notna().sum()
        bio_avail[b] = avail
        print(f"  {b}: {avail}/{len(nhanes)} ({avail/len(nhanes)*100:.0f}%)")

# Strategy:
# Core DM: tc, hdl, creat, sbp, dbp, bmi (6 markers, >85% coverage)
# Extended DM: core + crp + glu (wave I only)

# ---- Core DM (6 markers, all waves) ----
core_markers = ['tc','hdl','creat','sbp','dbp','bmi']
core_vars = ['ln_tc','ln_hdl','ln_creat','sbp','dbp','bmi']

# Impute and transform
for m in core_markers:
    nhanes[m] = pd.to_numeric(nhanes[m], errors='coerce')
    nhanes[m] = nhanes[m].fillna(nhanes[m].median())

nhanes['ln_tc'] = np.log(nhanes['tc'])
nhanes['ln_hdl'] = np.log(nhanes['hdl'])
nhanes['ln_creat'] = np.log(nhanes['creat'])

# Exclude implausible values
nhanes = nhanes[(nhanes['bmi']>=15)&(nhanes['bmi']<=60)&
                (nhanes['sbp']>=70)&(nhanes['sbp']<=250)&
                (nhanes['dbp']>=30)&(nhanes['dbp']<=150)]

# Healthy reference: 45-55y, BMI 18.5-25, no DM
healthy_nh = nhanes[(nhanes['age']>=45)&(nhanes['age']<=55)&
                    (nhanes['bmi']>=18.5)&(nhanes['bmi']<25)&
                    (nhanes['dm_bl'].fillna(0)==0)]

# Compute mu, Sigma from healthy reference
X_ref_nh = healthy_nh[core_vars].dropna()
mu_nh = X_ref_nh.mean().values
Sigma_nh = X_ref_nh.cov().values
Sigma_inv_nh = np.linalg.inv(Sigma_nh)

# DM for all
X_all_nh = nhanes[core_vars].values
DM_nh = np.array([np.sqrt((x-mu_nh).T @ Sigma_inv_nh @ (x-mu_nh))
                  for x in X_all_nh])
nhanes['DM'] = DM_nh
nhanes['DM_z'] = (DM_nh - np.nanmean(DM_nh)) / np.nanstd(DM_nh)
nhanes['DM_tertile'] = pd.qcut(nhanes['DM_z'], 3, labels=['T1','T2','T3'], duplicates='drop')

print(f"\nCore DM (6 markers): mean={np.nanmean(DM_nh):.2f}, SD={np.nanstd(DM_nh):.2f}")
print(f"Healthy ref N: {len(healthy_nh)}")

# ---- H1: Cox regression ----
print("\n--- H1: NHANES Cox ---")

def run_cox_nh(df, event='died', exposure='DM_z', covars=['age','male']):
    cph = CoxPHFitter()
    results = {}
    for model_name, cov_list in [
        ('Unadjusted', []),
        ('+Age+Sex', ['age','male']),
        ('Full', covars),
    ]:
        vars_needed = [exposure,'fu_years',event] + cov_list
        mdf = df[vars_needed].dropna()
        try:
            cph.fit(mdf, 'fu_years', event)
            hr = np.exp(cph.params_[exposure])
            ci = np.exp(cph.confidence_intervals_.loc[exposure].values)
            results[model_name] = {'HR':hr, 'CI_low':ci[0], 'CI_high':ci[1]}
        except:
            results[model_name] = {'HR':np.nan, 'CI_low':np.nan, 'CI_high':np.nan}
    return results

# All-cause
res_acm = run_cox_nh(nhanes, 'died', 'DM_z', ['age','male'])
print("All-cause mortality:")
for k,v in res_acm.items():
    print(f"  {k}: HR={v['HR']:.3f} (95%CI {v['CI_low']:.3f}-{v['CI_high']:.3f})")

# CVD death
res_cvd = run_cox_nh(nhanes, 'cvd_death', 'DM_z', ['age','male'])
print("CVD mortality:")
for k,v in res_cvd.items():
    print(f"  {k}: HR={v['HR']:.3f} (95%CI {v['CI_low']:.3f}-{v['CI_high']:.3f})")

# ---- H2: Discordance ----
print("\n--- H2: Discordance ---")
nhanes['trad_high'] = ((nhanes['age']>=65) | (nhanes['dm_bl'].fillna(0)==1) |
                       (nhanes['sbp']>=140)).astype(int)
nhanes['cmin_high'] = (nhanes['DM_tertile']=='T3').astype(int)

nhanes['risk_group'] = np.where(
    (nhanes['trad_high']==0)&(nhanes['cmin_high']==0), 'Concordant-Low',
    np.where((nhanes['trad_high']==0)&(nhanes['cmin_high']==1), 'Discordant',
    np.where((nhanes['trad_high']==1)&(nhanes['cmin_high']==0), 'Trad-Hi/CMIN-Lo',
    'Concordant-High')))

for g in ['Concordant-Low','Discordant','Trad-Hi/CMIN-Lo','Concordant-High']:
    sub = nhanes[nhanes['risk_group']==g]
    d_all = int(sub['died'].sum())
    d_cvd = int(sub['cvd_death'].sum())
    print(f"  {g}: N={len(sub)}, all-cause={d_all} ({d_all/len(sub)*100:.1f}%), CVD={d_cvd} ({d_cvd/len(sub)*100:.1f}%)")

# Discordant vs Concordant-Low Cox
disc_nh = nhanes[nhanes['risk_group'].isin(['Discordant','Concordant-Low'])]
disc_nh['disc'] = (disc_nh['risk_group']=='Discordant').astype(int)
cph_d = CoxPHFitter()
try:
    cph_d.fit(disc_nh[['disc','fu_years','died','age','male']].dropna(), 'fu_years', 'died')
    hr_d = np.exp(cph_d.params_['disc'])
    ci_d = np.exp(cph_d.confidence_intervals_.loc['disc'].values)
    print(f"\nDiscordant vs Concordant-Low: HR={hr_d:.3f} (95%CI {ci_d[0]:.3f}-{ci_d[1]:.3f})")
except Exception as e:
    print(f"\nDiscordant Cox: {e}")

# ---- C-index improvement ----
print("\n--- C-index ---")
cph_b = CoxPHFitter()
bdf = nhanes[['fu_years','died','age','male','sbp']].dropna()
cph_b.fit(bdf, 'fu_years', 'died')
c_base = concordance_index(bdf['fu_years'], -cph_b.predict_partial_hazard(bdf), bdf['died'])

edf = nhanes[['fu_years','died','age','male','sbp','DM_z']].dropna()
cph_e = CoxPHFitter()
cph_e.fit(edf, 'fu_years', 'died')
c_ext = concordance_index(edf['fu_years'], -cph_e.predict_partial_hazard(edf), edf['died'])
print(f"C-index: base={c_base:.4f}, base+DM={c_ext:.4f}, Δ={c_ext-c_base:.4f}")

# ---- Quartile dose-response ----
print("\n--- Dose-response (quartiles) ---")
nhanes['DM_q'] = pd.qcut(nhanes['DM_z'], 4, labels=['Q1','Q2','Q3','Q4'], duplicates='drop')
for i in [2,3,4]:
    nhanes[f'Q{i}'] = (nhanes['DM_q']==f'Q{i}').astype(int)
qdf = nhanes[['Q2','Q3','Q4','fu_years','died','age','male']].dropna()
cph_q = CoxPHFitter()
try:
    cph_q.fit(qdf, 'fu_years', 'died')
    for qi in ['Q2','Q3','Q4']:
        hr_q = np.exp(cph_q.params_[qi])
        ci_q = np.exp(cph_q.confidence_intervals_.loc[qi].values)
        print(f"  {qi} vs Q1: HR={hr_q:.2f} (95%CI {ci_q[0]:.2f}-{ci_q[1]:.2f})")
except Exception as e:
    print(f"  Quartile Cox failed: {e}")

# ---- CHARLS vs NHANES comparison ----
print(f"\n{'='*60}")
print(f"COMPARISON: CHARLS vs NHANES")
print(f"{'='*60}")
# Load CHARLS results
ch_res = pickle.load(open(os.path.join(OUT,'charls_results.pkl'),'rb'))
print(f"CHARLS: N={ch_res['n_total']}, deaths={ch_res['n_deaths']}, HR={ch_res['HR_DM_z_per_SD_allcause']:.2f}")
print(f"NHANES: N={len(nhanes)}, deaths={int(nhanes['died'].sum())}, HR={res_acm['Full']['HR']:.2f}")

# ---- Save ----
nh_results = {
    'n_total': len(nhanes),
    'n_deaths_allcause': int(nhanes['died'].sum()),
    'n_deaths_cvd': int(nhanes['cvd_death'].sum()),
    'n_healthy_ref': len(healthy_nh),
    'DM_mean': float(np.nanmean(DM_nh)),
    'DM_sd': float(np.nanstd(DM_nh)),
    'HR_acm_perSD': res_acm['Full']['HR'],
    'HR_acm_CI_low': res_acm['Full']['CI_low'],
    'HR_acm_CI_high': res_acm['Full']['CI_high'],
    'HR_cvd_perSD': res_cvd['Full']['HR'],
    'HR_cvd_CI_low': res_cvd['Full']['CI_low'],
    'HR_cvd_CI_high': res_cvd['Full']['CI_high'],
    'C_index_base': float(c_base),
    'C_index_plusDM': float(c_ext),
}
pickle.dump(nh_results, open(os.path.join(OUT,'nhanes_results.pkl'),'wb'))
nhanes.to_pickle(os.path.join(OUT,'nhanes_clean.pkl'))

print(f"\nNHANES validation complete. Results saved.")
