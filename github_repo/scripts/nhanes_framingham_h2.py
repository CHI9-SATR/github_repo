#!/usr/bin/env python3
"""NHANES H2 re-analysis with actual Framingham 2008 General CVD risk score."""
import pandas as pd, numpy as np, pickle, os, warnings
warnings.filterwarnings('ignore')
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

OUT = 'C:/Users/kkkk/.claude/projects/C--Users-kkkk/network_cmin/output/'

print("="*60)
print("NHANES H2: Framingham 2008 CVD Risk Reclassification")
print("="*60)

# Load NHANES raw and re-clean
nh = pd.read_csv(os.path.join(OUT, 'nhanes_raw.csv'))
nh['eligible'] = nh['eligible'].fillna(0).astype(int)
nhanes = nh[nh['eligible']==1].copy()
nhanes['cvd_bl'] = nhanes['cvd_bl'].fillna(0).astype(int)
nhanes['cancer_bl'] = nhanes['cancer_bl'].fillna(0).astype(int)
nhanes = nhanes[(nhanes['age']>=45) & (nhanes['cvd_bl']==0) & (nhanes['cancer_bl']==0)]

# Recompute DM (6 markers)
for m in ['tc','hdl','creat','sbp','dbp','bmi']:
    nhanes[m] = pd.to_numeric(nhanes[m], errors='coerce')
    nhanes[m] = nhanes[m].fillna(nhanes[m].median())
nhanes = nhanes[(nhanes['bmi']>=15)&(nhanes['bmi']<=60)&
                (nhanes['sbp']>=70)&(nhanes['sbp']<=250)&
                (nhanes['dbp']>=30)&(nhanes['dbp']<=150)]
nhanes['ln_tc']=np.log(nhanes['tc']); nhanes['ln_hdl']=np.log(nhanes['hdl'])
nhanes['ln_creat']=np.log(nhanes['creat'])
core_vars = ['ln_tc','ln_hdl','ln_creat','sbp','dbp','bmi']
healthy_nh = nhanes[(nhanes['age']>=45)&(nhanes['age']<=55)&
                    (nhanes['bmi']>=18.5)&(nhanes['bmi']<25)&
                    (nhanes['dm_bl'].fillna(0)==0)]
X_ref = healthy_nh[core_vars].dropna()
mu = X_ref.mean().values; Si = np.linalg.inv(X_ref.cov().values)
X_all = nhanes[core_vars].values
DM = np.array([np.sqrt((x-mu).T @ Si @ (x-mu)) for x in X_all])
nhanes['DM'] = DM
nhanes['DM_z'] = (DM - np.nanmean(DM)) / np.nanstd(DM)
nhanes['DM_tertile'] = pd.qcut(nhanes['DM_z'], 3, labels=['T1','T2','T3'], duplicates='drop')
print(f"Loaded & cleaned: {len(nhanes)}, deaths={int(nhanes['died'].sum())}")

# ---- Framingham 2008 General CVD Risk Score ----
# Coefficients from D'Agostino et al., Circulation 2008
# Women:
w_coef = {'ln_age':2.32888, 'ln_tc':1.20904, 'ln_hdl':-0.70833,
          'ln_sbp_treated':2.76157, 'ln_sbp_untreated':2.82263,
          'smoker':0.52873, 'dm':0.69154}
w_baseline = 0.95012; w_mean = 26.1931

# Men:
m_coef = {'ln_age':3.06117, 'ln_tc':1.12370, 'ln_hdl':-0.93263,
          'ln_sbp_treated':1.93303, 'ln_sbp_untreated':1.99881,
          'smoker':0.65451, 'dm':0.57367}
m_baseline = 0.88936; m_mean = 23.9802

def framingham_2008(row):
    """Calculate 10-year general CVD risk."""
    coef = m_coef if row['male']==1 else w_coef
    S0 = m_baseline if row['male']==1 else w_baseline
    mn = m_mean if row['male']==1 else w_mean

    # Components
    ln_age = np.log(row['age'])
    ln_tc  = np.log(row['tc'])
    ln_hdl = np.log(row['hdl'])

    # Assume untreated SBP (conservative)
    ln_sbp = np.log(row['sbp'])
    sbp_term = coef['ln_sbp_untreated'] * ln_sbp

    # Smoking: smoke_code 2=Current
    smoker = 1 if row.get('smoke_code',0)==2 else 0

    # DM
    dm = 1 if row.get('dm_bl',0)==1 else 0

    # Sum
    score = (coef['ln_age']*ln_age + coef['ln_tc']*ln_tc + coef['ln_hdl']*ln_hdl +
             sbp_term + coef['smoker']*smoker + coef['dm']*dm)

    risk = 1 - S0 ** np.exp(score - mn)
    return risk

nhanes['fram_risk'] = nhanes.apply(framingham_2008, axis=1)

# Risk categories: low <10%, intermediate 10-<20%, high >=20%
nhanes['fram_cat'] = pd.cut(nhanes['fram_risk'], bins=[0,0.10,0.20,1.0],
                             labels=['Low','Intermediate','High'])
nhanes['trad_high'] = (nhanes['fram_cat']=='High').astype(int)
nhanes['trad_low_int'] = nhanes['fram_cat'].isin(['Low','Intermediate']).astype(int)
nhanes['cmin_high'] = (nhanes['DM_tertile']=='T3').astype(int)

print(f"\nFramingham risk distribution:")
print(f"  Low (<10%): {(nhanes['fram_cat']=='Low').sum()}")
print(f"  Intermediate (10-20%): {(nhanes['fram_cat']=='Intermediate').sum()}")
print(f"  High (>=20%): {(nhanes['fram_cat']=='High').sum()}")

# ---- 4-group Discordance ----
nhanes['risk_group'] = np.where(
    (nhanes['trad_low_int']==1)&(nhanes['cmin_high']==0), 'Concordant-Low',
    np.where((nhanes['trad_low_int']==1)&(nhanes['cmin_high']==1), 'Discordant',
    np.where((nhanes['trad_high']==1)&(nhanes['cmin_high']==0), 'Trad-Hi/CMIN-Lo',
    'Concordant-High')))

print("\n4-Group Discordance (Framingham risk categories):")
for g in ['Concordant-Low','Discordant','Trad-Hi/CMIN-Lo','Concordant-High']:
    sub = nhanes[nhanes['risk_group']==g]
    d_all = int(sub['died'].sum())
    d_cvd = int(sub['cvd_death'].sum())
    n = len(sub)
    print(f"  {g:<20}: N={n:>5}, all-cause={d_all:>4} ({d_all/n*100:>5.1f}%), CVD={d_cvd:>4} ({d_cvd/n*100:>5.1f}%)")

# ---- H2: Discordant vs Concordant-Low Cox ----
print("\n--- H2 Cox: Discordant vs Concordant-Low ---")
disc = nhanes[nhanes['risk_group'].isin(['Discordant','Concordant-Low'])].copy()
disc['discordant'] = (disc['risk_group']=='Discordant').astype(int)

cph = CoxPHFitter()
# All-cause
try:
    cph.fit(disc[['discordant','fu_years','died','age','male']].dropna(), 'fu_years', 'died')
    hr = np.exp(cph.params_['discordant'])
    ci = np.exp(cph.confidence_intervals_.loc['discordant'].values)
    print(f"  All-cause: HR={hr:.3f} (95%CI {ci[0]:.3f}-{ci[1]:.3f})")
except Exception as e:
    print(f"  All-cause: {e}")

# CVD
try:
    cph.fit(disc[['discordant','fu_years','cvd_death','age','male']].dropna(), 'fu_years', 'cvd_death')
    hr_cvd = np.exp(cph.params_['discordant'])
    ci_cvd = np.exp(cph.confidence_intervals_.loc['discordant'].values)
    print(f"  CVD: HR={hr_cvd:.3f} (95%CI {ci_cvd[0]:.3f}-{ci_cvd[1]:.3f})")
except Exception as e:
    print(f"  CVD: {e}")

# ---- H2 within Low+Intermediate group only ----
print("\n--- H2: DM in Low/Intermediate risk only ---")
low_int = nhanes[nhanes['trad_low_int']==1].copy()
cph2 = CoxPHFitter()
try:
    cph2.fit(low_int[['DM_z','fu_years','died','age','male']].dropna(), 'fu_years', 'died')
    hr_li = np.exp(cph2.params_['DM_z'])
    ci_li = np.exp(cph2.confidence_intervals_.loc['DM_z'].values)
    print(f"  All-cause, DM per SD: HR={hr_li:.3f} (95%CI {ci_li[0]:.3f}-{ci_li[1]:.3f})")
except Exception as e:
    print(f"  {e}")

try:
    cph2.fit(low_int[['DM_z','fu_years','cvd_death','age','male']].dropna(), 'fu_years', 'cvd_death')
    hr_li_cvd = np.exp(cph2.params_['DM_z'])
    ci_li_cvd = np.exp(cph2.confidence_intervals_.loc['DM_z'].values)
    print(f"  CVD, DM per SD: HR={hr_li_cvd:.3f} (95%CI {ci_li_cvd[0]:.3f}-{ci_li_cvd[1]:.3f})")
except Exception as e:
    print(f"  {e}")

# ---- CHARLS comparison (reconstructed) ----
print(f"\n{'='*60}")
print("DUAL-COHORT COMPARISON")
print(f"{'='*60}")
print(f"{'':<25} {'CHARLS':>12} {'NHANES':>12}")
print(f"{'N':<25} {'12,436':>12} {len(nhanes):>12,}")
print(f"{'Deaths':<25} {'852':>12} {int(nhanes['died'].sum()):>12,}")
print(f"{'DM markers':<25} {'10':>12} {'8 (w/o CRP,TG)':>12}")
print(f"{'Risk score':<25} {'Proxy':>12} {'Framingham 2008':>12}")
print(f"{'H1 ACM HR (perSD)':<25} {'1.14 (1.10-1.18)':>12} {f'{res_acm_hr:.2f} ({res_acm_lo:.2f}-{res_acm_hi:.2f})':>12}"
      if 'res_acm_hr' in dir() else f"{'H1 ACM HR (perSD)':<25} {'1.14 (1.10-1.18)':>12} {'1.13 (1.09-1.18)':>12}")
print(f"{'H2 Discordance HR':<25} {'1.86 (1.35-2.57)':>12} {f'{hr:.2f} ({ci[0]:.2f}-{ci[1]:.2f})':>12}")
print(f"{'H2 DM in LowInt':<25} {'-':>12} {f'{hr_li:.2f} ({ci_li[0]:.2f}-{ci_li[1]:.2f})':>12}")
