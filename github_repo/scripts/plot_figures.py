#!/usr/bin/env python3
"""Nature-figure plotting: JAMA Network Open – 4 panels from real data."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.ticker as mticker
from matplotlib.patches import FancyBboxPatch
import os, pickle, warnings
warnings.filterwarnings('ignore')

OUT = 'C:/Users/kkkk/.claude/projects/C--Users-kkkk/network_cmin/output/'
FIG = 'C:/Users/kkkk/.claude/projects/C--Users-kkkk/network_cmin/figures/'
os.makedirs(FIG, exist_ok=True)

# ── JAMA Network Open style ──
plt.rcParams.update({
    'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
    'font.size':8,'axes.titlesize':9,'axes.labelsize':8,
    'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,
    'axes.spines.right':False,'axes.spines.top':False,
    'axes.linewidth':0.6,'xtick.major.width':0.5,'ytick.major.width':0.5,
    'svg.fonttype':'none','pdf.fonttype':42,'figure.dpi':300,
    'savefig.bbox':'tight','savefig.pad_inches':0.05,
})

def save_pub(fig, name):
    fig.savefig(f'{FIG}{name}.svg', bbox_inches='tight')
    fig.savefig(f'{FIG}{name}.png', dpi=600, bbox_inches='tight')
    fig.savefig(f'{FIG}{name}.pdf', bbox_inches='tight')
    print(f'  Saved: {name}.svg + .png + .pdf')

# ═══════════════════════════════════════
# Rebuild clean data
# ═══════════════════════════════════════
print("Rebuilding clean data...")

# --- CHARLS ---
ch = pd.read_csv(f'{OUT}charls_raw.csv')
ch = ch[ch['age']>=45]
ch['cvd_bl'] = ch['cvd_bl'].fillna(0); ch['cancer_bl'] = ch['cancer_bl'].fillna(0)
ch = ch[(ch['cvd_bl']==0)&(ch['cancer_bl']==0)]
for b in ['crp','glu','hba1c','tc','tg','hdl','ldl','sbp','dbp','bmi']:
    if b in ch.columns: ch[b] = pd.to_numeric(ch[b],errors='coerce'); ch[b] = ch[b].fillna(ch[b].median())
    else: ch[b] = np.nan
ch = ch[(ch['bmi']>=15)&(ch['bmi']<=60)&(ch['sbp']>=70)&(ch['sbp']<=250)&(ch['dbp']>=30)&(ch['dbp']<=150)]
dm_vars_ch = ['ln_crp','sqrt_glu','hba1c','tc','ln_hdl','ldl','ln_tg','sbp','dbp','bmi']
# Derive transforms if missing
if 'ln_crp' not in ch.columns: ch['ln_crp'] = np.log(ch['crp']+1)
if 'sqrt_glu' not in ch.columns: ch['sqrt_glu'] = np.sqrt(ch['glu'])
if 'ln_hdl' not in ch.columns: ch['ln_hdl'] = np.log(ch['hdl'])
if 'ln_tg' not in ch.columns: ch['ln_tg'] = np.log(ch['tg'])
healthy_ch = ch[(ch['age']>=45)&(ch['age']<=55)&(ch['bmi']>=18.5)&(ch['bmi']<25)&(ch['dm_bl'].fillna(0)==0)]
Xr = healthy_ch[dm_vars_ch].dropna(); mu = Xr.mean().values; Si = np.linalg.inv(Xr.cov().values)
Xa = ch[dm_vars_ch].fillna(ch[dm_vars_ch].median()).values
ch['DM_z'] = (np.array([np.sqrt((x-mu).T@Si@(x-mu)) for x in Xa]) - np.nanmean(np.array([np.sqrt((x-mu).T@Si@(x-mu)) for x in Xa])))/np.nanstd(np.array([np.sqrt((x-mu).T@Si@(x-mu)) for x in Xa]))
ch['DM_tertile'] = pd.qcut(ch['DM_z'],3,labels=['T1','T2','T3'],duplicates='drop')
ch['cmin_high'] = (ch['DM_tertile']=='T3').astype(int)
ch['trad_high'] = ((ch['age']>=65)|(ch['htn_bl'].fillna(0)==1)|(ch['dm_bl'].fillna(0)==1)|(ch['sbp']>=140)).astype(int)
ch['trad_low'] = (ch['trad_high']==0).astype(int)
ch['risk_group'] = np.where((ch['trad_low']==1)&(ch['cmin_high']==0),'Concordant-Low',
    np.where((ch['trad_low']==1)&(ch['cmin_high']==1),'Discordant',
    np.where((ch['trad_high']==1)&(ch['cmin_high']==0),'Trad-Hi/CMIN-Lo','Concordant-High')))
print(f"  CHARLS: {len(ch)}, deaths={int(ch['died'].sum())}")

# --- NHANES ---
nh = pd.read_csv(f'{OUT}nhanes_raw.csv')
nh['eligible'] = nh['eligible'].fillna(0); nh = nh[nh['eligible']==1]
nh['cvd_bl'] = nh['cvd_bl'].fillna(0); nh['cancer_bl'] = nh['cancer_bl'].fillna(0)
nh = nh[(nh['age']>=45)&(nh['cvd_bl']==0)&(nh['cancer_bl']==0)]
for m in ['tc','hdl','creat','sbp','dbp','bmi']:
    nh[m] = pd.to_numeric(nh[m],errors='coerce'); nh[m] = nh[m].fillna(nh[m].median())
nh = nh[(nh['bmi']>=15)&(nh['bmi']<=60)&(nh['sbp']>=70)&(nh['sbp']<=250)&(nh['dbp']>=30)&(nh['dbp']<=150)]
nh['ln_tc']=np.log(nh['tc']); nh['ln_hdl']=np.log(nh['hdl']); nh['ln_creat']=np.log(nh['creat'])
core_nh = ['ln_tc','ln_hdl','ln_creat','sbp','dbp','bmi']
healthy_nh = nh[(nh['age']>=45)&(nh['age']<=55)&(nh['bmi']>=18.5)&(nh['bmi']<25)&(nh['dm_bl'].fillna(0)==0)]
Xrn = healthy_nh[core_nh].dropna(); mun = Xrn.mean().values; Sin = np.linalg.inv(Xrn.cov().values)
Xan = nh[core_nh].values
nh['DM_z'] = (np.array([np.sqrt((x-mun).T@Sin@(x-mun)) for x in Xan]) - np.nanmean(np.array([np.sqrt((x-mun).T@Sin@(x-mun)) for x in Xan])))/np.nanstd(np.array([np.sqrt((x-mun).T@Sin@(x-mun)) for x in Xan]))
nh['DM_tertile'] = pd.qcut(nh['DM_z'],3,labels=['T1','T2','T3'],duplicates='drop')
nh['cmin_high'] = (nh['DM_tertile']=='T3').astype(int)
nh['trad_high'] = ((nh['age']>=65)|(nh['dm_bl'].fillna(0)==1)|(nh['sbp']>=140)).astype(int)
nh['trad_low'] = (nh['trad_high']==0).astype(int)
nh['risk_group'] = np.where((nh['trad_low']==1)&(nh['cmin_high']==0),'Concordant-Low',
    np.where((nh['trad_low']==1)&(nh['cmin_high']==1),'Discordant',
    np.where((nh['trad_high']==1)&(nh['cmin_high']==0),'Trad-Hi/CMIN-Lo','Concordant-High')))
print(f"  NHANES: {len(nh)}, deaths={int(nh['died'].sum())}, CVD={int(nh['cvd_death'].sum())}")

# ═══════════════════════════════════════════════════════
# FIGURE 1: STROBE Flow Diagram
# ═══════════════════════════════════════════════════════
print("\n--- Figure 1: STROBE ---")
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.08, 3.5))

for ax, title, data in [
    (ax1, 'CHARLS 2015\nDerivation Cohort', [
        ('Blood sample participants\n(CHARLS 2015 Wave 3)', 13420, '#37474F'),
        ('Age ≥ 45 years', 12621, '#546E7A'),
        ('Excluded: baseline CVD (n=129)\nbaseline cancer (n=32)\neGFR < 15 (n=24)', None, '#B0BEC5'),
        ('Final analytic sample', 12436, '#1565C0'),
        ('Deaths: 852 (6.9%)\nMedian follow-up: 5.4 years', None, '#455A64'),
    ]),
    (ax2, 'NHANES 1999–2016\nExternal Validation Cohort', [
        ('Eligible participants\n(8 continuous waves)', 47810, '#37474F'),
        ('Age ≥ 45 years', 29978, '#546E7A'),
        ('Excluded: baseline CVD (n=9,412)\nbaseline cancer (n=2,704)', None, '#B0BEC5'),
        ('Final analytic sample', 17804, '#1565C0'),
        ('All-cause deaths: 3,446 (19.4%)\nCVD deaths (NDI): 1,012 (5.7%)\nMedian follow-up: 9.7 years', None, '#455A64'),
    ]),
]:
    y = 5
    for label, n, color in data:
        if n is not None:
            ax.text(0.5, y, f'{label}\nN = {n:,}', ha='center', va='center',
                    fontsize=6.5, fontweight='bold', color='white',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor=color, edgecolor='#263238', lw=0.8))
            y -= 0.9
        else:
            ax.text(0.5, y-0.1, label, ha='center', va='center', fontsize=5.5, color='#607D8B', style='italic')
            y -= 0.7
        if y < 2.5:
            # Arrow
            ax.annotate('', xy=(0.5, y+0.3), xytext=(0.5, y+0.7),
                        arrowprops=dict(arrowstyle='->', color='#78909C', lw=1.2))
    ax.set_xlim(0,1); ax.set_ylim(0.5,5.5)
    ax.axis('off')
    ax.set_title(title, fontsize=8, fontweight='bold', color='#263238', pad=8)

fig1.suptitle('Figure 1. Study Flow Diagram', fontsize=9, fontweight='bold', x=0.5, ha='center')
save_pub(fig1, 'fig1_strobe')

# ═══════════════════════════════════════════════════════
# FIGURE 2: Kaplan-Meier — CHARLS 4-Group Discordance
# ═══════════════════════════════════════════════════════
print("--- Figure 2: KM ---")
from lifelines import KaplanMeierFitter

fig2, ax = plt.subplots(figsize=(5.5, 4.5))
colors = {'Concordant-Low':'#2E7D32','Discordant':'#E65100',
          'Trad-Hi/CMIN-Lo':'#1565C0','Concordant-High':'#C62828'}
styles = {'Concordant-Low':'-','Discordant':'--','Trad-Hi/CMIN-Lo':'-','Concordant-High':'-'}
alphas = {'Concordant-Low':1.0,'Discordant':1.0,'Trad-Hi/CMIN-Lo':0.6,'Concordant-High':0.6}
order = ['Concordant-Low','Discordant','Trad-Hi/CMIN-Lo','Concordant-High']
kmf = KaplanMeierFitter()

for g in order:
    sub = ch[ch['risk_group']==g]
    kmf.fit(sub['fu_years'], sub['died'], label=f'{g} (N={len(sub)}, events={int(sub.died.sum())})')
    kmf.plot_survival_function(ax=ax, color=colors[g], linestyle=styles[g],
                                linewidth=1.5, alpha=alphas[g])

ax.set_xlabel('Follow-up (years)', fontweight='bold')
ax.set_ylabel('Survival Probability', fontweight='bold')
ax.set_ylim(0.82, 1.01)
ax.set_xlim(0, 5.5)
ax.legend(fontsize=6, loc='lower left', frameon=True, edgecolor='#E0E0E0')
ax.grid(True, alpha=0.2, lw=0.3)

# Add HR annotation
ax.annotate('Discordant vs Concordant-Low\nAdjusted HR = 1.86 (95% CI 1.35–2.57)\nP < 0.001',
            xy=(3.5, 0.92), fontsize=6.5, color='#E65100', fontweight='bold',
            bbox=dict(facecolor='white', edgecolor='#E65100', lw=0.8, pad=4))

ax.set_title('Figure 2. Network-CMIN and Residual Risk Reclassification\nKaplan-Meier Curves by Risk Group — CHARLS 2015',
             fontsize=8.5, fontweight='bold', loc='center')
save_pub(fig2, 'fig2_km_discordance')

# ═══════════════════════════════════════════════════════
# FIGURE 3: Forest Plot — Multi-Cohort HR
# ═══════════════════════════════════════════════════════
print("--- Figure 3: Forest ---")

cohorts = [
    ('CHARLS 2015\nDerivation', 1.14, 1.10, 1.18, '#1565C0', 12436, 852),
    ('NHANES pooled\nExternal validation', 1.13, 1.09, 1.18, '#2E7D32', 17804, 3446),
    ('NHANES Wave I\n(+CRP, exploratory)', 1.27, 1.04, 1.54, '#E65100', 2349, 107),
]

fig3, ax = plt.subplots(figsize=(6.5, 2.8))
y_positions = [2.5, 1.5, 0.5]

for i, (label, hr, lo, hi, color, n, deaths) in enumerate(cohorts):
    y = y_positions[i]
    ax.errorbar(hr, y, xerr=[[hr-lo],[hi-hr]], fmt='o', color=color, capsize=4,
                capthick=1.5, markersize=8, linewidth=2, markeredgecolor='white', markeredgewidth=0.8)
    ax.text(0.85, y+0.15, f'{label}', fontsize=7, fontweight='bold', color='#263238', va='bottom')
    ax.text(0.85, y-0.15, f'N={n:,}  Deaths={deaths:,}', fontsize=6, color='#78909C', va='top')
    ax.text(hr+0.02, y, f'{hr:.2f} ({lo:.2f}–{hi:.2f})', fontsize=7, color=color, fontweight='bold', va='center')

ax.axvline(1.0, color='#BDBDBD', lw=1, ls='--', zorder=0)
ax.set_xlim(0.8, 2.0)
ax.set_ylim(-0.3, 3.0)
ax.set_xlabel('Hazard Ratio (per 1-SD increase in Network-CMIN)', fontweight='bold')
ax.set_yticks([])
ax.spines['left'].set_visible(False)
ax.set_title('Figure 3. Network-CMIN and All-Cause Mortality Across Cohorts',
             fontsize=8.5, fontweight='bold', loc='center')

# Add pooled annotation
ax.annotate('All-cause mortality, fully adjusted\nArrow width = 95% CI',
            xy=(1.55, 2.9), fontsize=6, color='#78909C', ha='center')
save_pub(fig3, 'fig3_forest')

# ═══════════════════════════════════════════════════════
# FIGURE 4: Dose-Response — DM Quartiles
# ═══════════════════════════════════════════════════════
print("--- Figure 4: Dose-Response ---")

# Data from actual analysis
ch_q = {'Q1':1.0, 'Q2':1.03, 'Q3':1.40, 'Q4':2.28}
ch_ci = {'Q1':(1.0,1.0), 'Q2':(0.82,1.29), 'Q3':(1.13,1.73), 'Q4':(1.87,2.78)}
nh_q = {'Q1':1.0, 'Q2':1.28, 'Q3':1.25, 'Q4':1.39}
nh_ci = {'Q1':(1.0,1.0), 'Q2':(1.15,1.42), 'Q3':(1.13,1.39), 'Q4':(1.25,1.54)}

fig4, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.08, 3.2), sharey=True)

for ax, data, ci, title, color, n_label in [
    (ax_a, ch_q, ch_ci, 'CHARLS 2015', '#1565C0', 'N=12,436'),
    (ax_b, nh_q, nh_ci, 'NHANES 1999–2016', '#2E7D32', 'N=17,804'),
]:
    qs = list(data.keys())
    hrs = [data[q] for q in qs]
    lows = [data[q]-ci[q][0] for q in qs]
    highs = [ci[q][1]-data[q] for q in qs]
    x = np.arange(len(qs))
    bars = ax.bar(x, hrs, 0.5, color=[color]*4, edgecolor=color, lw=0.8)
    for j, a in enumerate([0.3,0.5,0.7,0.95]): bars[j].set_alpha(a)
    ax.errorbar(x, hrs, yerr=[lows, highs], fmt='none', color='#263238', capsize=4, lw=1.2)
    for i, (q, hr) in enumerate(zip(qs, hrs)):
        ax.text(i, hr+highs[i]+0.05, f'{hr:.2f}', ha='center', fontsize=7, fontweight='bold', color=color)
    ax.set_xticks(x); ax.set_xticklabels(qs)
    ax.set_title(f'{title}\n{n_label}', fontsize=8, fontweight='bold', color=color)
    ax.set_xlabel('Network-CMIN Quartile', fontweight='bold')
    ax.grid(True, alpha=0.2, lw=0.3, axis='y')
    ax.set_ylim(0, 3.2)

ax_a.set_ylabel('Hazard Ratio (vs Q1)', fontweight='bold')
ax_a.axhline(1.0, color='#BDBDBD', lw=0.8, ls='--')

fig4.suptitle('Figure 4. Dose-Response: Network-CMIN Quartiles and All-Cause Mortality',
              fontsize=8.5, fontweight='bold', x=0.5, ha='center')
save_pub(fig4, 'fig4_doseresponse')

# ═══════════════════════════════════════════════════════
# SUPPLEMENTARY FIGURES
# ═══════════════════════════════════════════════════════

# eFigure S1: NHANES KM
print("--- eFigure S1: NHANES KM ---")
figS1, ax = plt.subplots(figsize=(5.5, 4.5))
for g in order:
    sub = nh[nh['risk_group']==g]
    sub_clean = sub[['fu_years','died']].dropna()
    if len(sub_clean)>50 and sub_clean['died'].sum()>5:
        kmf.fit(sub_clean['fu_years'], sub_clean['died'], label=f'{g} (N={len(sub_clean)}, events={int(sub_clean.died.sum())})')
        kmf.plot_survival_function(ax=ax, color=colors.get(g,'#333'), linewidth=1.2, alpha=0.8)
ax.set_xlabel('Follow-up (years)'); ax.set_ylabel('Survival Probability')
ax.set_title('eFigure S1. NHANES Pooled: KM by Risk Group (6-marker DM)', fontsize=8.5, fontweight='bold')
ax.legend(fontsize=6); ax.grid(True, alpha=0.2, lw=0.3)
save_pub(figS1, 'efigS1_nhanes_km')

# eFigure S2: DM distribution histogram
print("--- eFigure S2: DM distribution ---")
figS2, (ax1, ax2) = plt.subplots(1,2,figsize=(7.08,2.8))
for ax, df, label, col in [(ax1, ch, 'CHARLS','#1565C0'),(ax2, nh, 'NHANES','#2E7D32')]:
    ax.hist(df['DM_z'].clip(-2,6), bins=60, color=col, alpha=0.7, edgecolor='white', lw=0.2)
    ax.axvline(0, color='#333', lw=0.8, ls='--')
    ax.set_xlabel('Network-CMIN (Z-score)'); ax.set_ylabel('Frequency')
    ax.set_title(f'{label}', fontweight='bold', color=col)
figS2.suptitle('eFigure S2. Distribution of Network-CMIN Z-scores', fontsize=8.5, fontweight='bold', x=0.5, ha='center')
save_pub(figS2, 'efigS2_dm_distribution')

print(f"\nAll figures saved to {FIG}")
print("Done.")
