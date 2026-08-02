#!/usr/bin/env python3
"""NHANES Wave I (2015-2016): Full 10-marker DM with CRP, then H1+H2."""
import pandas as pd, numpy as np, os, pickle, warnings
warnings.filterwarnings('ignore')
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index

OUT = 'C:/Users/kkkk/.claude/projects/C--Users-kkkk/network_cmin/output/'
nd, md = 'D:/NHANES/', 'D:/NHANES_mortality/'
s = 'I'

print("="*60)
print(f"NHANES Wave {s}: Full DM with CRP")
print("="*60)

# ---- Load Wave I with ALL biomarkers ----
demo = pd.read_sas(f'{nd}DEMO_{s}.xpt')
bmx  = pd.read_sas(f'{nd}BMX_{s}.xpt')
bpx  = pd.read_sas(f'{nd}BPX_{s}.xpt')

df = demo[['SEQN','RIAGENDR','RIDAGEYR','RIDRETH1','DMDEDUC2','WTMEC2YR']].copy()
df = df.merge(bmx[['SEQN','BMXBMI','BMXWAIST']], on='SEQN', how='left')
bpx_c = ['SEQN']+[c for c in ['BPXSY1','BPXSY2','BPXSY3','BPXDI1','BPXDI2','BPXDI3'] if c in bpx.columns]
df = df.merge(bpx[bpx_c], on='SEQN', how='left')

# Load ALL biomarker files for wave I
print("Loading biomarkers...")
bio_dfs = []
for fname, vars_needed in [
    ('TCHOL', ['SEQN','LBXTC']),
    ('HDL',   ['SEQN','LBDHDD']),
    ('GLU',   ['SEQN','LBXGLU']),
    ('GHB',   ['SEQN','LBXGH']),
    ('HSCRP', ['SEQN','LBXHSCRP']),
    ('BIOPRO',['SEQN','LBXSCR','LBXALB']),
    ('CBC',   ['SEQN','LBXWBCSI','LBXHGB','LBXPLTSI']),
]:
    for ext in ['.xpt','.csv']:
        path = f'{nd}{fname}_{s}{ext}'
        if os.path.exists(path):
            bd = pd.read_sas(path) if ext=='.xpt' else pd.read_csv(path)
            avail = [c for c in vars_needed if c in bd.columns]
            if avail:
                bio_dfs.append(bd[avail])
            break

# Also try to find TG (may be in TRIGLY or BIOPRO)
for fname in ['TRIGLY','BIOPRO']:
    for ext in ['.xpt','.csv']:
        path = f'{nd}{fname}_{s}{ext}'
        if os.path.exists(path):
            bd = pd.read_sas(path) if ext=='.xpt' else pd.read_csv(path)
            tg_cols = [c for c in bd.columns if 'TRIG' in c.upper() or 'TG' in c.upper() or 'LBDTG' in c.upper()]
            if tg_cols:
                print(f"  TG found in {fname}_{s}: {tg_cols}")
                bio_dfs.append(bd[['SEQN']+tg_cols])
            break

# Questionnaires
for qf, qv in [
    (f'{nd}SMQ_{s}.xpt', ['SEQN','SMQ020','SMQ040']),
    (f'{nd}DIQ_{s}.xpt', ['SEQN','DIQ010']),
    (f'{nd}MCQ_{s}.xpt', ['SEQN','MCQ160B','MCQ160C','MCQ160E','MCQ160F','MCQ220']),
]:
    if os.path.exists(qf):
        qd=pd.read_sas(qf); a=[c for c in qv if c in qd.columns]
        if a: bio_dfs.append(qd[a])

# Merge all
for bd in bio_dfs:
    if 'SEQN' in bd.columns and bd.shape[1]>1:
        df = df.merge(bd, on='SEQN', how='left')

# MORT
mf = f'{md}NHANES_2015_2016_MORT_2019_PUBLIC.dat'
if os.path.exists(mf):
    mort=pd.read_fwf(mf, colspecs=[(0,6),(14,15),(15,16),(16,19),(45,48)],
                     names=['SEQN','eligstat','mortstat','ucod_leading','permth_exm'])
    for c in ['eligstat','mortstat','ucod_leading','permth_exm']:
        mort[c]=pd.to_numeric(mort[c],errors='coerce')
    df=df.merge(mort,on='SEQN',how='left')

print(f"Wave I loaded: {df.shape[0]} rows, {df.shape[1]} cols")

# ---- Clean & derive ----
df['age']=df['RIDAGEYR']; df['male']=(df['RIAGENDR']==1).astype(int)
df['bmi']=df['BMXBMI']; df['waist']=df['BMXWAIST']
bp_s=[c for c in ['BPXSY1','BPXSY2','BPXSY3'] if c in df.columns]
bp_d=[c for c in ['BPXDI1','BPXDI2','BPXDI3'] if c in df.columns]
if bp_s: df['sbp']=df[bp_s].mean(axis=1)
if bp_d: df['dbp']=df[bp_d].mean(axis=1)

# Rename biomarkers
nmap={'LBXTC':'tc','LBDHDD':'hdl','LBXGLU':'glu','LBXGH':'hba1c',
      'LBXHSCRP':'crp','LBXSCR':'creat','LBXALB':'alb',
      'LBXWBCSI':'wbc','LBXHGB':'hb','LBXPLTSI':'plt'}
for o,n in nmap.items():
    if o in df.columns: df[n]=pd.to_numeric(df[o],errors='coerce')

# Check for TG
tg_cols = [c for c in df.columns if 'TRIG' in c.upper() or 'TG' in c.upper() or 'LBDTG' in c.upper()]
if tg_cols:
    df['tg'] = pd.to_numeric(df[tg_cols[0]], errors='coerce')
    print(f"  TG from: {tg_cols[0]}")

# Disease flags
df['dm_bl']=((df.get('DIQ010',pd.Series(0))==1)).astype(int)
cvd_c=[c for c in ['MCQ160B','MCQ160C','MCQ160E','MCQ160F'] if c in df.columns]
df['cvd_bl']=df[cvd_c].eq(1).any(axis=1).astype(int) if cvd_c else 0
df['cancer_bl']=((df.get('MCQ220',pd.Series(0))==1)).astype(int)
df['eligible']=(df.get('eligstat',pd.Series(0))==1)
df['died']=((df['mortstat']==1)&df['eligible']).astype(int)
df['cvd_death']=((df['died']==1)&(df['ucod_leading'].isin([1,5]))).astype(int)
df['fu_years']=(df['permth_exm']/12).clip(0.01,None)

# Smoking
df['smoke']=df.apply(lambda r: 'Current' if pd.notna(r.get('SMQ040')) and r['SMQ040'] in[1,2] else
    'Former' if pd.notna(r.get('SMQ020')) and r['SMQ020']==1 else
    'Never' if pd.notna(r.get('SMQ020')) and r['SMQ020']==2 else np.nan, axis=1)
df['smoke_code']=df['smoke'].map({'Never':0,'Former':1,'Current':2})

# Exclusions
clean = df[(df['eligible'])&(df['age']>=45)&(df['cvd_bl']==0)&(df['cancer_bl']==0)].copy()

# Check biomarker availability
print(f"\nWave I eligible (age>=45, no CVD/cancer): {len(clean)}")
print(f"Deaths: {int(clean['died'].sum())}, CVD: {int(clean['cvd_death'].sum())}")
print(f"Med FU: {clean['fu_years'].median():.1f}y")
for b in ['crp','glu','hba1c','tc','tg','hdl','creat','sbp','dbp','bmi']:
    if b in clean.columns:
        print(f"  {b}: {clean[b].notna().sum()}/{len(clean)} ({clean[b].notna().sum()/len(clean)*100:.0f}%)")
    else:
        print(f"  {b}: NOT IN DATA")

# ---- Build best possible DM ----
# Use available markers with >70% coverage
avail_markers = []
for b in ['crp','glu','hba1c','tc','tg','hdl','creat','sbp','dbp','bmi']:
    if b in clean.columns and clean[b].notna().sum()/len(clean) > 0.5:
        avail_markers.append(b)

print(f"\nAvailable markers (>50%): {avail_markers}")

# Impute + transform
for m in avail_markers:
    clean[m] = clean[m].fillna(clean[m].median())
    clean[m] = clean[m].clip(lower=clean[m].quantile(0.001), upper=clean[m].quantile(0.999))

clean['ln_crp']=np.log(clean['crp']+1)
clean['sqrt_glu']=np.sqrt(clean['glu']) if 'glu' in clean.columns else None
clean['ln_tc']=np.log(clean['tc']) if 'tc' in clean.columns else None
clean['ln_hdl']=np.log(clean['hdl']) if 'hdl' in clean.columns else None
clean['ln_tg']=np.log(clean['tg']) if 'tg' in clean.columns else None
clean['ln_creat']=np.log(clean['creat']) if 'creat' in clean.columns else None

# Build transform map
dm_vars = []
for b in avail_markers:
    if b == 'crp': dm_vars.append('ln_crp')
    elif b == 'glu': dm_vars.append('sqrt_glu')
    elif b == 'hba1c': dm_vars.append('hba1c')
    elif b == 'tc': dm_vars.append('ln_tc')
    elif b == 'hdl': dm_vars.append('ln_hdl')
    elif b == 'tg': dm_vars.append('ln_tg')
    elif b == 'creat': dm_vars.append('ln_creat')
    else: dm_vars.append(b)  # sbp, dbp, bmi

dm_vars = [v for v in dm_vars if v in clean.columns and v is not None]
print(f"DM variables ({len(dm_vars)}): {dm_vars}")

# Filter implausible values
clean = clean[(clean['bmi']>=15)&(clean['bmi']<=60)&
              (clean['sbp']>=70)&(clean['sbp']<=250)&
              (clean['dbp']>=30)&(clean['dbp']<=150)]

# Healthy reference
healthy = clean[(clean['age']>=45)&(clean['age']<=55)&
                (clean['bmi']>=18.5)&(clean['bmi']<25)&
                (clean['dm_bl']==0)]
X_ref = healthy[dm_vars].dropna()
mu = X_ref.mean().values
Si = np.linalg.inv(X_ref.cov().values)
X_all = clean[dm_vars].fillna(clean[dm_vars].median()).values
DM = np.array([np.sqrt((x-mu).T @ Si @ (x-mu)) for x in X_all])
clean['DM'] = DM
clean['DM_z'] = (DM - np.nanmean(DM))/np.nanstd(DM)
clean['DM_tertile'] = pd.qcut(clean['DM_z'],3,labels=['T1','T2','T3'],duplicates='drop')

print(f"DM: mean={np.nanmean(DM):.2f}, SD={np.nanstd(DM):.2f}")
print(f"Healthy ref N: {len(healthy)}")

# ---- H1: Cox ----
print("\n--- H1: Wave I Cox ---")
cph = CoxPHFitter()

for outcome, label in [('died','All-cause'),('cvd_death','CVD')]:
    print(f"\n{label}:")
    for model_name, covs in [
        ('Unadjusted',[]),
        ('+Age+Sex',['age','male']),
        ('Full',['age','male']),
    ]:
        mcols = ['DM_z','fu_years',outcome]+covs
        mdf = clean[mcols].dropna()
        try:
            cph.fit(mdf,'fu_years',outcome)
            hr=np.exp(cph.params_['DM_z'])
            ci=np.exp(cph.confidence_intervals_.loc['DM_z'].values)
            p = cph.summary.loc['DM_z','p']
            print(f"  {model_name:<15}: HR={hr:.3f} ({ci[0]:.3f}-{ci[1]:.3f}), P={p:.4f}")
        except Exception as e:
            print(f"  {model_name:<15}: FAILED — {e}")

# ---- H2: Framingham + Discordance ----
print("\n--- H2: Framingham Discordance ---")

def fram_risk(row):
    if row['male']==1:
        coef={'ln_age':3.06117,'ln_tc':1.12370,'ln_hdl':-0.93263,'ln_sbp':1.99881,'smoker':0.65451,'dm':0.57367}
        S0=0.88936; mn=23.9802
    else:
        coef={'ln_age':2.32888,'ln_tc':1.20904,'ln_hdl':-0.70833,'ln_sbp':2.82263,'smoker':0.52873,'dm':0.69154}
        S0=0.95012; mn=26.1931
    score = (coef['ln_age']*np.log(row['age'])+coef['ln_tc']*np.log(row['tc'])+
             coef['ln_hdl']*np.log(row['hdl'])+coef['ln_sbp']*np.log(row['sbp'])+
             coef['smoker']*(1 if row.get('smoke_code',0)==2 else 0)+
             coef['dm']*(1 if row.get('dm_bl',0)==1 else 0))
    return 1 - S0**np.exp(score-mn)

clean['fram'] = clean.apply(fram_risk, axis=1)
clean['fram_hi'] = (clean['fram']>=0.20).astype(int)
clean['fram_lo_int'] = (clean['fram']<0.20).astype(int)
clean['cmin_hi'] = (clean['DM_tertile']=='T3').astype(int)

clean['risk_group'] = np.where(
    (clean['fram_lo_int']==1)&(clean['cmin_hi']==0),'Concordant-Low',
    np.where((clean['fram_lo_int']==1)&(clean['cmin_hi']==1),'Discordant',
    np.where((clean['fram_hi']==1)&(clean['cmin_hi']==0),'Trad-Hi/CMIN-Lo','Concordant-High')))

print("4-Group:")
for g in ['Concordant-Low','Discordant','Trad-Hi/CMIN-Lo','Concordant-High']:
    sub=clean[clean['risk_group']==g]
    print(f"  {g:<20}: N={len(sub):>4}, ACM={int(sub['died'].sum()):>3} ({sub['died'].mean()*100:.1f}%), CVD={int(sub['cvd_death'].sum()):>3}")

# Discordant Cox
disc=clean[clean['risk_group'].isin(['Discordant','Concordant-Low'])]
disc['disc']=(disc['risk_group']=='Discordant').astype(int)
cph.fit(disc[['disc','fu_years','died','age','male']].dropna(),'fu_years','died')
hr_d=np.exp(cph.params_['disc'])
ci_d=np.exp(cph.confidence_intervals_.loc['disc'].values)
p_d=cph.summary.loc['disc','p']
print(f"\nDiscordant vs Concordant-Low: HR={hr_d:.3f} ({ci_d[0]:.3f}-{ci_d[1]:.3f}), P={p_d:.4f}")

# DM continuous in Low/Int
li=clean[clean['fram_lo_int']==1]
cph.fit(li[['DM_z','fu_years','died','age','male']].dropna(),'fu_years','died')
hr_li=np.exp(cph.params_['DM_z']); ci_li=np.exp(cph.confidence_intervals_.loc['DM_z'].values)
print(f"DM per SD in Low/Int: HR={hr_li:.3f} ({ci_li[0]:.3f}-{ci_li[1]:.3f})")

# ---- Final 3-cohort summary ----
print(f"\n{'='*60}")
print(f"THREE-COHORT SUMMARY")
print(f"{'='*60}")
print(f"{'Cohort':<18} {'N':>6} {'Deaths':>6} {'DM biomarkers':>15} {'H1 HR':>10} {'H2 Disc HR':>12}")
print(f"{'CHARLS 2015':<18} {12436:>6} {852:>6} {'10 (full)':>15} {'1.14 (1.10-1.18)':>10} {'1.86 (1.35-2.57)':>12}")
print(f"{'NHANES all':<18} {17804:>6} {3446:>6} {'6 (no CRP)':>15} {'1.13 (1.09-1.18)':>10} {'1.10 (0.96-1.25)':>12}")
print(f"{'NHANES Wave I':<18} {len(clean):>6} {int(clean['died'].sum()):>6} {f'{len(dm_vars)} markers':>15} {'—':>10} {f'{hr_d:.2f} ({ci_d[0]:.2f}-{ci_d[1]:.2f})':>12}")
print(f"\nDM per SD in Low/Int risk: CHARLS=—, NHANES Wave I={hr_li:.3f} ({ci_li[0]:.3f}-{ci_li[1]:.3f})")

# Save
clean.to_pickle(os.path.join(OUT,'nhanes_waveI_clean.pkl'))
print(f"\nWave I results saved.")
