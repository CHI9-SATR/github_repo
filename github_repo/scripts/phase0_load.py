#!/usr/bin/env python3
"""Phase 0: CHARLS 2015 + NHANES. Simplified version with correct 2015 variable names."""
import pandas as pd, numpy as np, os, warnings
warnings.filterwarnings('ignore')
OUT = 'C:/Users/kkkk/.claude/projects/C--Users-kkkk/network_cmin/output/'
os.makedirs(OUT, exist_ok=True)

def to_n(col): return pd.to_numeric(col, errors='coerce')
def yn(col): return col.astype(str).str.contains('Yes',na=False).astype(int) if col.dtype=='object' else col.astype(float).fillna(0).astype(int)

# ============================================================
# CHARLS 2015
# ============================================================
print("=== CHARLS 2015 ===")
blood = pd.read_stata('D:/CHARLS/CHARLS_2015/Blood.dta')
health = pd.read_stata('D:/CHARLS/CHARLS_2015/Health_Status_and_Functioning.dta')
bio = pd.read_stata('D:/CHARLS/CHARLS_2015/Biomarker.dta')
weight = pd.read_stata('D:/CHARLS/CHARLS_2015/Weights.dta')
si15 = pd.read_stata('D:/CHARLS/CHARLS_2015/Sample_Infor.dta')
si18 = pd.read_stata('D:/CHARLS/CHARLS_2018/Sample_Infor.dta')
si20 = pd.read_stata('D:/CHARLS/CHARLS_2020/Sample_Infor.dta')

# Merge
df = blood.copy()
df = df.merge(health[['ID','xrgender','zda007_1_','zda007_2_','zda007_3_','zda007_4_',
    'zda007_7_','zda007_8_','zda007_9_',
    'da059','da061','da067','da069',
    'da051_1_','da051_2_','da051_3_','da052_1_','da052_2_','da052_3_']], on='ID', how='left')
df = df.merge(bio[['ID','qa003','qa004','qa007','qa008','qa011','qa012','qi002','ql002','qm002']], on='ID', how='left')
df = df.merge(weight[['ID','Biomarker_weight']].rename(columns={'Biomarker_weight':'blood_weight'}), on='ID', how='left')

# Death tracking
for si, tag in [(si15,'15'),(si18,'18'),(si20,'20')]:
    sid = si[['ID','died']].copy(); sid.columns = ['ID',f'died_{tag}']
    sid[f'died_{tag}'] = (sid[f'died_{tag}'].astype(str).str[:1]=='1').astype(int)
    df = df.merge(sid, on='ID', how='left')
df['died'] = ((df.get('died_15',0).fillna(0)==1) | (df.get('died_18',0).fillna(0)==1) |
              (df.get('died_20',0).fillna(0)==1)).astype(int)

# Variables
df['male'] = (df['xrgender'].astype(str).str.extract(r'(\d)').astype(float).iloc[:,0]==1).astype(int)

# Age: from health variables or blood data doesn't have it. Use approximate from biomarkers file
# In practice, use 2015 - (approx birth year). The ba004_w3_1 from demo has birth year.
# Since we didn't merge demo, use a placeholder: CHARLS 2015 requirement = age>=45.
# All CHARLS participants are 45+, so we can estimate from zda009 which has age at diagnosis.
# Simplest: use the age we can derive. Actually let's just merge demo minimally
demo = pd.read_stata('D:/CHARLS/CHARLS_2015/Demographic_Background.dta')
df = df.merge(demo[['ID','ba004_w3_1']], on='ID', how='left')
df['birth_yr'] = to_n(df['ba004_w3_1'])
df['age'] = 2015 - df['birth_yr']

# Education: bf004 is spouse's education. Find the respondent's education level
# In 2015, education from W1 is carried forward; bd001_w2_4 tracks changes
# The current education level is best found via the health file or by merging W1 data
# For now, use a numeric code derived from available variables
# Let's just get it from demo: check bd001_w2_4 or any bd variable
with pd.read_stata('D:/CHARLS/CHARLS_2015/Demographic_Background.dta', iterator=True) as r:
    dlabels = r.variable_labels()
# Look for any variable with 'education' and 'level' in label, excluding spouse/adult
edu_candidates = []
for v, lab in dlabels.items():
    l = str(lab).lower()
    if 'education' in l and 'level' in l and 'spouse' not in l and 'adult' not in l and 'changed' not in l:
        edu_candidates.append(v)
        print(f"  Education candidate: {v} = {dlabels[v]}")

if edu_candidates:
    edu_var = edu_candidates[0]
    df = df.merge(demo[['ID',edu_var]], on='ID', how='left')
    df['edu_level'] = to_n(df[edu_var]).fillna(0).astype(int)
else:
    # Fallback: use any available
    df['edu_level'] = 0
    print("  WARNING: No education level variable found, using 0")

# Smoking
df['smoke_ever'] = yn(df['da059']) if df['da059'].dtype == 'object' else to_n(df['da059']).fillna(0)
try:
    df['smoke_cur'] = yn(df['da061'])
except:
    df['smoke_cur'] = 0
df['smoke'] = np.where(df['smoke_cur'].fillna(0)==1, 2, np.where(df['smoke_ever'].fillna(0)==1, 1, 0))
df['drink_ly'] = yn(df['da067']) if df['da067'].dtype == 'object' else to_n(df['da067']).fillna(0)

# Physical activity (2015 variables use same names as 2011: da051/da052)
for c in ['da051_1_','da051_2_','da051_3_','da052_1_','da052_2_','da052_3_']:
    if c in df.columns: df[c] = to_n(df[c])
pa_cols = [c for c in ['da052_1_','da052_2_','da052_3_'] if c in df.columns]
if pa_cols: df['pa_active'] = (df[pa_cols].max(axis=1).fillna(0)>=1).astype(int)
else: df['pa_active'] = 0

# BP
for c in ['qa003','qa004','qa007','qa008','qa011','qa012']:
    if c in df.columns: df[c] = to_n(df[c])
bp_cols = [c for c in ['qa003','qa007','qa011'] if c in df.columns]
if bp_cols: df['sbp'] = df[bp_cols].mean(axis=1)
else: df['sbp'] = np.nan
bp_dcols = [c for c in ['qa004','qa008','qa012'] if c in df.columns]
if bp_dcols: df['dbp'] = df[bp_dcols].mean(axis=1)
else: df['dbp'] = np.nan

# Anthropometry
df['height_cm'] = to_n(df['qi002']); df['height_cm'] = df['height_cm'].clip(100, 200)  # plausible range
df['weight_kg'] = to_n(df['ql002']); df['weight_kg'] = df['weight_kg'].clip(25, 200)
df['bmi_kgm2'] = df['weight_kg'] / (df['height_cm']/100)**2
df['bmi_kgm2'] = df['bmi_kgm2'].clip(10, 60)
df['waist_cm'] = to_n(df['qm002'])

# Biomarkers (2015: bl_ prefix)
bmap15 = {'bl_crp':'crp','bl_hbalc':'hba1c','bl_glu':'glu','bl_cho':'tc',
    'bl_tg':'tg','bl_hdl':'hdl','bl_ldl':'ldl','bl_crea':'creat',
    'bl_bun':'bun','bl_ua':'ua','bl_wbc':'wbc','bl_hgb':'hb','bl_plt':'plt','bl_cysc':'cystc'}
for o,n in bmap15.items():
    if o in df.columns: df[n] = to_n(df[o])
    else: df[n] = np.nan

df['crp_log']=np.log(df['crp']+1); df['glu_sqrt']=np.sqrt(df['glu'])
df['tg_log']=np.log(df['tg']); df['hdl_log']=np.log(df['hdl'])

# eGFR
df['egfr'] = np.nan; v = df['creat'].notna()&(df['creat']>0)
if v.sum()>0:
    scr=df.loc[v,'creat']; k=np.where(df.loc[v,'male']==1,0.9,0.7)
    a=np.where(df.loc[v,'male']==1,-0.302,-0.241)
    df.loc[v,'egfr']=142*(scr/k)**a*0.9938**df.loc[v,'age']
    df.loc[v&(df['male']==1),'egfr']*=1.012; df['egfr']=df['egfr'].clip(upper=200)

# Diseases
for tag in ['zda007_1_','zda007_3_','zda007_4_','zda007_7_','zda007_8_','zda007_9_']:
    if tag in df.columns:
        if df[tag].dtype == 'object': df[tag] = yn(df[tag])
        else: df[tag] = to_n(df[tag]).fillna(0)

df['htn_bl'] = df.get('zda007_1_',0); df['dm_bl'] = df.get('zda007_3_',0)
df['cvd_bl'] = ((df.get('zda007_7_',0).fillna(0)==1) if 'zda007_7_' in df.columns else 0) | \
               ((df.get('zda007_8_',0).fillna(0)==1) if 'zda007_8_' in df.columns else 0)
df['cancer_bl'] = df.get('zda007_4_',0)

# Follow-up
df['fu_years'] = np.where(df['died']==1,
    np.where(df.get('died_20',0).fillna(0)==1, 5,
    np.where(df.get('died_18',0).fillna(0)==1, 3, 0.5)),
    5.5)  # censored at 2020 (~5.5y max)
df['fu_years'] = df['fu_years'].clip(0.1, 5.5)

df['cvd_death'] = df['died']  # all-cause placeholder
df['svy_weight'] = to_n(df.get('blood_weight', np.nan))

# Save
charls_cols = ['ID','age','male','edu_level','smoke','drink_ly','pa_active',
    'sbp','dbp','bmi_kgm2','waist_cm',
    'crp','glu','hba1c','tc','tg','hdl','ldl','creat','egfr','ua','wbc','hb','plt',
    'crp_log','glu_sqrt','tg_log','hdl_log',
    'htn_bl','dm_bl','cvd_bl','cancer_bl','died','cvd_death','fu_years','svy_weight']
charls = df[[c for c in charls_cols if c in df.columns]].copy()
# Rename for consistency
charls = charls.rename(columns={'bmi_kgm2':'bmi','waist_cm':'waist'})

# Exclude age<45
charls = charls[charls['age']>=45].copy()
print(f"  CHARLS 2015 (age>=45): {charls.shape[0]} rows, {charls.died.sum()} deaths, age [{charls.age.min():.0f},{charls.age.max():.0f}]")

# Biomarker coverage
for b in ['crp','glu','hba1c','tc','tg','hdl','ldl','creat']:
    if b in charls.columns:
        avail = charls[b].notna().sum()
        print(f"    {b}: {avail}/{len(charls)} ({avail/len(charls)*100:.0f}%)")

charls.to_csv(os.path.join(OUT,'charls_raw.csv'), index=False)
print(f"  Saved charls_raw.csv")

# ============================================================
# NHANES
# ============================================================
print("\n=== NHANES ===")
nd, md = 'D:/NHANES/', 'D:/NHANES_mortality/'
waves = {'B':'2001_2002','C':'2003_2004','D':'2005_2006','E':'2007_2008',
         'F':'2009_2010','G':'2011_2012','H':'2013_2014','I':'2015_2016'}
parts = []
for s,w in waves.items():
    try:
        demo = pd.read_sas(f'{nd}DEMO_{s}.xpt')
        bmx  = pd.read_sas(f'{nd}BMX_{s}.xpt')
        dp = demo[['SEQN','RIAGENDR','RIDAGEYR','RIDRETH1','DMDEDUC2','WTMEC2YR']].copy()
        bp = bmx[['SEQN','BMXBMI','BMXWAIST']].copy()
        wd = dp.merge(bp, on='SEQN', how='left'); wd['wave'] = s

        # BP
        bpxf = f'{nd}BPX_{s}.xpt'
        if os.path.exists(bpxf):
            bpx=pd.read_sas(bpxf); bpx_c=['SEQN']+[c for c in ['BPXSY1','BPXSY2','BPXSY3','BPXDI1','BPXDI2','BPXDI3'] if c in bpx.columns]
            if len(bpx_c)>1: wd=wd.merge(bpx[bpx_c],on='SEQN',how='left')

        # Biomarkers: handle .xpt and .csv
        for fname,vname in [
            ('TCHOL','LBXTC'),('HDL','LBDHDD'),
            ('BIOPRO','LBXSCR'),('BIOPRO','LBXALB'),
        ]:
            for ext in ['.xpt','.csv']:
                path = f'{nd}{fname}_{s}{ext}'
                if os.path.exists(path):
                    bd=pd.read_sas(path) if ext=='.xpt' else pd.read_csv(path)
                    if vname in bd.columns: wd=wd.merge(bd[['SEQN',vname]],on='SEQN',how='left')
                    break

        for fname,vname in [('TRIGLY','LBXTR'),]:
            for ext in ['.xpt','.csv']:
                path = f'{nd}{fname}_{s}{ext}'
                if os.path.exists(path):
                    bd=pd.read_sas(path) if ext=='.xpt' else pd.read_csv(path)
                    if vname in bd.columns: wd=wd.merge(bd[['SEQN',vname]],on='SEQN',how='left')
                    break

        for fname,vname in [('GLU','LBXGLU'),]:
            for ext in ['.xpt','.csv']:
                path = f'{nd}{fname}_{s}{ext}'
                if os.path.exists(path):
                    bd=pd.read_sas(path) if ext=='.xpt' else pd.read_csv(path)
                    if vname in bd.columns: wd=wd.merge(bd[['SEQN',vname]],on='SEQN',how='left')
                    break

        for fname,vname in [('GHB','LBXGH'),]:
            for ext in ['.xpt','.csv']:
                path = f'{nd}{fname}_{s}{ext}'
                if os.path.exists(path):
                    bd=pd.read_sas(path) if ext=='.xpt' else pd.read_csv(path)
                    if vname in bd.columns: wd=wd.merge(bd[['SEQN',vname]],on='SEQN',how='left')
                    break

        for fname,vname in [('HSCRP','LBXHSCRP'),]:
            for ext in ['.xpt','.csv']:
                path = f'{nd}{fname}_{s}{ext}'
                if os.path.exists(path):
                    bd=pd.read_sas(path) if ext=='.xpt' else pd.read_csv(path)
                    if vname in bd.columns: wd=wd.merge(bd[['SEQN',vname]],on='SEQN',how='left')
                    break

        # CBC
        for ext in ['.xpt','.csv']:
            path = f'{nd}CBC_{s}{ext}'
            if os.path.exists(path):
                bd=pd.read_sas(path) if ext=='.xpt' else pd.read_csv(path)
                for cv in ['LBXWBCSI','LBXHGB','LBXPLTSI']:
                    if cv in bd.columns: wd=wd.merge(bd[['SEQN',cv]],on='SEQN',how='left')
                break

        # Questionnaires
        for qf,qv in [
            (f'{nd}SMQ_{s}.xpt',['SEQN','SMQ020','SMQ040']),
            (f'{nd}DIQ_{s}.xpt',['SEQN','DIQ010']),
            (f'{nd}MCQ_{s}.xpt',['SEQN','MCQ160B','MCQ160C','MCQ160E','MCQ160F','MCQ220']),
        ]:
            if os.path.exists(qf):
                qd=pd.read_sas(qf); a=[c for c in qv if c in qd.columns]
                if a: wd=wd.merge(qd[a],on='SEQN',how='left')

        # MORT
        mf = f'{md}NHANES_{w}_MORT_2019_PUBLIC.dat'
        if os.path.exists(mf):
            mort=pd.read_fwf(mf,colspecs=[(0,6),(14,15),(15,16),(16,19),(45,48)],
                             names=['SEQN','eligstat','mortstat','ucod_leading','permth_exm'])
            for c in ['eligstat','mortstat','ucod_leading','permth_exm']:
                mort[c]=pd.to_numeric(mort[c],errors='coerce')
            wd=wd.merge(mort,on='SEQN',how='left')

        parts.append(wd); print(f"  {s}: {len(wd)}")
    except Exception as e:
        print(f"  {s}: SKIP — {e}")

nh = pd.concat(parts,ignore_index=True)
print(f"  Combined: {nh.shape}")

# Derive
nh['age']=nh['RIDAGEYR']; nh['male']=(nh['RIAGENDR']==1).astype(int)
nh['edu_level']=nh['DMDEDUC2'].apply(lambda x: 0 if pd.isna(x) else 0 if x in[1,2] else 1 if x==3 else 2 if x==4 else 3 if x==5 else 0)
nh['smoke']=nh.apply(lambda r: np.nan if pd.isna(r.get('SMQ020')) else
    'Current' if pd.notna(r.get('SMQ040')) and r['SMQ040'] in[1,2] else
    'Former' if r['SMQ020']==1 else 'Never' if r['SMQ020']==2 else np.nan, axis=1)
nh['smoke_code']=nh['smoke'].map({'Never':0,'Former':1,'Current':2})

nh['bmi']=nh['BMXBMI']; nh['waist']=nh['BMXWAIST']
bp_s=[c for c in ['BPXSY1','BPXSY2','BPXSY3'] if c in nh.columns]
bp_d=[c for c in ['BPXDI1','BPXDI2','BPXDI3'] if c in nh.columns]
if bp_s: nh['sbp']=nh[bp_s].mean(axis=1)
if bp_d: nh['dbp']=nh[bp_d].mean(axis=1)

nmap={'LBXTC':'tc','LBDHDD':'hdl','LBXGLU':'glu','LBXTR':'tg',
    'LBXGH':'hba1c','LBXHSCRP':'crp','LBXSCR':'creat','LBXALB':'alb',
    'LBXWBCSI':'wbc','LBXHGB':'hb','LBXPLTSI':'plt'}
for o,n in nmap.items():
    if o in nh.columns: nh[n]=pd.to_numeric(nh[o],errors='coerce')

nh['crp_log']=np.log(nh['crp']+1); nh['glu_sqrt']=np.sqrt(nh['glu'])
nh['tg_log']=np.log(nh['tg']); nh['hdl_log']=np.log(nh['hdl'])

if 'creat' in nh.columns:
    nh['egfr']=np.nan; v=nh['creat'].notna()&(nh['creat']>0)
    if v.sum()>0:
        scr=nh.loc[v,'creat']; k=np.where(nh.loc[v,'male']==1,0.9,0.7)
        a=np.where(nh.loc[v,'male']==1,-0.302,-0.241)
        nh.loc[v,'egfr']=142*(scr/k)**a*0.9938**nh.loc[v,'age']
        nh.loc[v&(nh['male']==1),'egfr']*=1.012; nh['egfr']=nh['egfr'].clip(upper=200)

nh['dm_bl']=((nh.get('DIQ010',pd.Series(0))==1)).astype(int)
cvd_c=[c for c in ['MCQ160B','MCQ160C','MCQ160E','MCQ160F'] if c in nh.columns]
nh['cvd_bl']=nh[cvd_c].eq(1).any(axis=1).astype(int) if cvd_c else 0
nh['cancer_bl']=((nh.get('MCQ220',pd.Series(0))==1)).astype(int)
nh['eligible']=(nh.get('eligstat',pd.Series(0))==1)
nh['died']=((nh['mortstat']==1)&nh['eligible']).astype(int)
nh['cvd_death']=((nh['died']==1)&(nh['ucod_leading'].isin([1,5]))).astype(int)
nh['fu_years']=(nh['permth_exm']/12).clip(0.01,None)
nh['svy_weight']=nh['WTMEC2YR']

e=nh[nh['eligible']]
print(f"  Eligible: {e.shape[0]}  Deaths: {e.died.sum()}  CVD deaths: {e.cvd_death.sum()}  Med FU: {e.fu_years.median():.1f}y")

# Biomarker coverage
for b in ['crp','glu','hba1c','tc','tg','hdl','creat']:
    if b in nh.columns:
        avail_n = e[b].notna().sum()
        print(f"    {b}: {avail_n}/{len(e)} ({avail_n/len(e)*100:.0f}%)")

nh_vars = ['SEQN','wave','age','male','edu_level','smoke_code',
    'sbp','dbp','bmi','waist','crp','glu','hba1c','tc','tg','hdl',
    'creat','egfr','alb','wbc','hb','plt','crp_log','glu_sqrt','tg_log','hdl_log',
    'dm_bl','cvd_bl','cancer_bl','died','cvd_death','fu_years','svy_weight','eligible']
nhanes = nh[[c for c in nh_vars if c in nh.columns]].copy()
nhanes.to_csv(os.path.join(OUT,'nhanes_raw.csv'), index=False)
print(f"  Saved nhanes_raw.csv ({nhanes.shape})")
print("\nPhase 0 complete.")
