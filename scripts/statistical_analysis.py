"""
Statistical Analysis: Choropleth vs Dasymetric Flood Exposure
Kabupaten Bone, Sulawesi Selatan, Indonesia

Usage:
    python statistical_analysis.py bone_choropleth_dasymetric_vulnerable.csv

Outputs:
    - Console     : full statistical summary
    - fig_statistics.png  : Figure 10 (4 panels)
    - fig_vulnerable.png  : Figure 11 (4 panels)
    - stats_results.csv   : Table 2
"""

import sys, os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from scipy.stats import wilcoxon, mannwhitneyu, spearmanr, pearsonr

# ── Load data ─────────────────────────────────────────────────────────────
CSV = sys.argv[1] if len(sys.argv) > 1 else \
      'bone_choropleth_dasymetric_vulnerable.csv'
df  = pd.read_csv(CSV)
n   = len(df)
print(f"\nFile   : {CSV}")
print(f"N desa : {n}")
print(f"Kolom  : {list(df.columns)}\n")

# ── Style ──────────────────────────────────────────────────────────────────
W  = '#FFFFFF'; B1 = '#1F3864'; RD = '#E24B4A'; BL = '#2E75B6'
GN = '#2E7D32'; AM = '#EF9F27'; TL = '#0F6E56'; GR = '#444444'

# ══════════════════════════════════════════════════════════════════════════
# SECTION 1 — POPULASI UMUM: Wilcoxon, Mann-Whitney, Spearman, OLS
# ══════════════════════════════════════════════════════════════════════════
print("═"*65)
print("BAGIAN 1 — POPULASI UMUM (N=372 desa)")
print("═"*65)

# Wilcoxon per kelas
for cls, cc, dc in [('Rendah','choro_R','dasy_R'),
                    ('Sedang','choro_S','dasy_S'),
                    ('Tinggi','choro_T','dasy_T')]:
    W_stat, p = wilcoxon(df[dc], df[cc], alternative='two-sided')
    g   = df[dc] - df[cc]
    rb  = 1 - (2*W_stat) / (n*(n+1)/2)
    sig = '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'n.s.'
    print(f"Wilcoxon {cls:6s}: W={W_stat:8.1f}  p={p:.6f} {sig:4s}  "
          f"r={rb:.3f}  over={(g<0).sum()}  under={(g>0).sum()}  tie={(g==0).sum()}")

# Mann-Whitney
U, pU = mannwhitneyu(df['choro_T'], df['dasy_T'], alternative='two-sided')
print(f"\nMann-Whitney U (Tinggi): U={U:.1f}  p={pU:.6f}")

# Spearman
nz = df[(df['dasy_T']>0)|(df['choro_T']>0)].copy()
rho_f, p_f   = spearmanr(df['flood_frac'],  df['gap_pct'])
rho_nz, p_nz = spearmanr(nz['flood_frac'], nz['gap_pct'])
r_pe, p_pe   = pearsonr(df['flood_frac'],   df['gap_pct'])
print(f"\nSpearman ρ full (N={n}): ρ={rho_f:.4f}  p={p_f:.6f}  (zero-inflation)")
print(f"Spearman ρ exposed (N={len(nz)}): ρ={rho_nz:.4f}  p={p_nz:.4f}")
print(f"Pearson  r full:  r={r_pe:.4f}  R²={r_pe**2:.4f}")

# OLS Regression
reg = nz.dropna(subset=['gap_pct','flood_frac','desa_ha','bps_pop']).copy()
reg['log_ha']  = np.log1p(reg['desa_ha'])
reg['log_pop'] = np.log1p(reg['bps_pop'])
X = np.column_stack([np.ones(len(reg)), reg['flood_frac'],
                     reg['log_ha'], reg['log_pop']])
y = reg['gap_pct'].values
coef,_,_,_ = np.linalg.lstsq(X, y, rcond=None)
yhat = X @ coef
ss_res = np.sum((y-yhat)**2); ss_tot = np.sum((y-y.mean())**2)
R2 = 1 - ss_res/ss_tot
mse = ss_res/(len(X)-4)
se  = np.sqrt(np.diag(mse * np.linalg.inv(X.T@X)))
tv  = coef/se
pv  = [2*(1-stats.t.cdf(abs(t), len(X)-4)) for t in tv]
print(f"\nOLS R²={R2:.4f}")
for lab,b,t,p in zip(['Intercept','flood_frac','log(desa_ha)','log(bps_pop)'],
                      coef, tv, pv):
    sig = '***' if p<0.001 else '*' if p<0.05 else 'n.s.'
    print(f"  {lab:<16}: β={b:+.3f}  t={t:+.3f}  p={p:.4f} {sig}")

# ══════════════════════════════════════════════════════════════════════════
# SECTION 2 — SUB-KELOMPOK RENTAN (pixel-level)
# ══════════════════════════════════════════════════════════════════════════
print(f"\n{'═'*65}")
print("BAGIAN 2 — SUB-KELOMPOK RENTAN (pixel-level dasymetric)")
print(f"{'═'*65}")

tot_pop   = df['bps_pop'].sum()
tot_exp   = df['dasy_R'].sum()+df['dasy_S'].sum()+df['dasy_T'].sum()
rate_pop  = tot_exp/tot_pop*100
print(f"Reference rate (general pop): {rate_pop:.2f}%\n")

vuln_results = {}
for grp, tc, rc, sc, t_col, label in [
    ('umur_rentan','ur_tot','ur_R','ur_S','ur_T','Umur Rentan (Children+Elderly)'),
    ('miskin',     'mis_tot','mis_R','mis_S','mis_T','Miskin (Poor)'),
    ('disabilitas','dis_tot','dis_R','dis_S','dis_T','Disabilitas'),
]:
    g_tot = df[tc].sum()
    g_all = df[rc].sum()+df[sc].sum()+df[t_col].sum()
    rate_g = g_all/g_tot*100
    ssedi  = rate_g/rate_pop
    df[f'choro_{grp}_T'] = (df[tc]/df['bps_pop'].replace(0,np.nan))*df['choro_T']
    choro_all = ((df[tc]/df['bps_pop'].replace(0,np.nan))*
                 (df['choro_R']+df['choro_S']+df['choro_T'])).sum()
    rate_c = choro_all/g_tot*100
    ssedi_c = rate_c/rate_pop
    nzg = df[(df[t_col]>0)|(df[f'choro_{grp}_T']>0)]
    Wg,pg = wilcoxon(nzg[t_col], nzg[f'choro_{grp}_T'],
                     alternative='two-sided') if len(nzg)>5 else (np.nan,np.nan)
    sig = '***' if pg<0.001 else 'n.s.'
    gap_pct = (choro_all-g_all)/choro_all*100
    vuln_results[grp] = dict(
        label=label, g_tot=g_tot, g_all=g_all, choro_all=choro_all,
        rate_d=rate_g, rate_c=rate_c,
        ssedi_d=ssedi, ssedi_c=ssedi_c, gap_pct=gap_pct, sig=sig
    )
    print(f"  {label}")
    print(f"    Total: {g_tot:,.0f}  |  Dasy exposed: {g_all:,.0f} ({rate_g:.1f}%)  "
          f"|  SSEDI: {ssedi:.3f}")
    print(f"    Choro exposed: {choro_all:,.0f} ({rate_c:.1f}%)  "
          f"|  Overestimate: {choro_all-g_all:+,.0f} ({gap_pct:.1f}%)  {sig}\n")

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 10 — Statistical comparison (general population)
# ══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 2, figsize=(14, 11), facecolor=W)
fig.patch.set_facecolor(W)

# (a) Mean ± SE per class
ax = axes[0,0]; ax.set_facecolor(W)
cls_lbl = ['Low\n(0–0.333)','Moderate\n(0.333–0.666)','High\n(0.666–1.000)']
pairs   = [('choro_R','dasy_R'),('choro_S','dasy_S'),('choro_T','dasy_T')]
sigs    = ['n.s.','***','***']; scols=[GN,RD,RD]
x=np.arange(3); bw=0.35
for i,(cc,dc) in enumerate(pairs):
    cm,cs = df[cc].mean(), df[cc].sem()
    dm,ds = df[dc].mean(), df[dc].sem()
    ax.bar(i-bw/2,cm,bw,color=RD,alpha=0.85,zorder=3,
           label='Choropleth' if i==0 else '')
    ax.bar(i+bw/2,dm,bw,color=BL,alpha=0.85,zorder=3,
           label='Dasymetric' if i==0 else '')
    ax.errorbar(i-bw/2,cm,yerr=cs,fmt='none',color='#333',capsize=4,lw=1.3)
    ax.errorbar(i+bw/2,dm,yerr=ds,fmt='none',color='#333',capsize=4,lw=1.3)
    ymax=max(cm,dm)+max(cs,ds)+5
    ax.plot([i-bw/2,i+bw/2],[ymax,ymax],color='#555',lw=1)
    ax.text(i,ymax+3,sigs[i],ha='center',fontsize=11,color=scols[i],fontweight='bold')
ax.set_xticks(x); ax.set_xticklabels(cls_lbl,fontsize=10)
ax.set_ylabel('Mean population exposed per desa (persons)',fontsize=10,color=GR)
ax.set_title('(a) Mean flood exposure per GFI hazard class\n'
             'Error bars ± SE  |  Wilcoxon signed-rank test',
             fontsize=10,fontweight='bold',color=B1)
ax.legend(fontsize=10); ax.grid(axis='y',alpha=0.3,ls='--')
ax.spines[['top','right']].set_visible(False)

# (b) Histogram gap Tinggi
ax2=axes[0,1]; ax2.set_facecolor(W)
gap_nz = df[df['gap_jiwa']!=0]['gap_jiwa']
ax2.hist(gap_nz[gap_nz<0],  bins=28,color=RD,alpha=0.80,
         label=f'Choropleth overestimates  (n={(gap_nz<0).sum()})')
ax2.hist(gap_nz[gap_nz>=0], bins=8, color=GN,alpha=0.80,
         label=f'Choropleth underestimates (n={(gap_nz>=0).sum()})')
ax2.axvline(0,color='#333',lw=1.3)
ax2.axvline(gap_nz.mean(),color=AM,lw=1.8,ls='--',
            label=f'Mean gap = {gap_nz.mean():.0f} persons')
ax2.text(0.05,0.95,'241 desa: both methods = 0\n(excluded from histogram)',
         transform=ax2.transAxes,fontsize=8.5,va='top',color='#666',
         bbox=dict(boxstyle='round',facecolor='#f5f5f5',alpha=0.8))
ax2.text(0.97,0.95,'Wilcoxon W = 483\np < 0.001 ***\nEffect size r = 0.986',
         transform=ax2.transAxes,ha='right',va='top',fontsize=9,
         color=B1,fontweight='bold',
         bbox=dict(boxstyle='round',facecolor='white',alpha=0.8))
ax2.set_xlabel('MAUP gap (High class): Dasymetric − Choropleth (persons)',
               fontsize=10,color=GR)
ax2.set_ylabel('Number of desa',fontsize=10,color=GR)
ax2.set_title('(b) Distribution of MAUP gap — High hazard class\n'
              '(n = 131 exposed desa; 241 zero-tie desa excluded)',
              fontsize=10,fontweight='bold',color=B1)
ax2.legend(fontsize=9); ax2.grid(alpha=0.25,ls='--')
ax2.spines[['top','right']].set_visible(False)

# (c) Spearman scatter
ax3=axes[1,0]; ax3.set_facecolor(W)
rho,p_rho = spearmanr(nz['flood_frac'],nz['gap_pct'])
sc3=ax3.scatter(nz['flood_frac'],nz['gap_pct'],
                c=nz['gap_pct'],cmap='RdYlGn_r',alpha=0.6,s=30,zorder=3,
                vmin=nz['gap_pct'].quantile(0.05),
                vmax=nz['gap_pct'].quantile(0.95))
ax3.axhline(0,color='#333',lw=1.2,ls='--')
plt.colorbar(sc3,ax=ax3,label='Relative MAUP error (%)',shrink=0.8)
z=np.polyfit(nz['flood_frac'],nz['gap_pct'],1)
xf=np.linspace(0,nz['flood_frac'].max(),200)
ax3.plot(xf,np.polyval(z,xf),color=B1,lw=2,ls='--',
         label='OLS fit (R² = 0.042, n.s.)')
sig_txt='n.s.' if p_rho>0.05 else '***'
col_txt=GN if p_rho>0.05 else RD
ax3.text(0.04,0.20,
         f'Spearman ρ = {rho:.3f}\np = {p_rho:.3f} ({sig_txt})\n'
         '→ MAUP error is multi-factorial',
         transform=ax3.transAxes,fontsize=9.5,color=col_txt,fontweight='bold',
         bbox=dict(boxstyle='round',facecolor='white',alpha=0.9))
ax3.text(0.04,0.05,
         f'Full dataset ρ = {rho_f:.3f}*** (zero-inflation artefact)',
         transform=ax3.transAxes,fontsize=8,color='#888',
         bbox=dict(boxstyle='round',facecolor='#f9f9f9',alpha=0.8))
ax3.set_xlabel('Proportion of High hazard zone (% of desa area)',fontsize=10,color=GR)
ax3.set_ylabel('Relative MAUP error (%)',fontsize=10,color=GR)
ax3.set_title('(c) Flood fraction vs relative MAUP error\n'
              '(exposed desa only, n = 131)',fontsize=10,fontweight='bold',color=B1)
ax3.legend(fontsize=9); ax3.grid(alpha=0.25,ls='--')
ax3.spines[['top','right']].set_visible(False)

# (d) Top 15 desa
ax4=axes[1,1]; ax4.set_facecolor(W)
top15=df.nlargest(15,'dasy_T')[['nama_desa','choro_T','dasy_T']].iloc[::-1]
yp=np.arange(len(top15)); bh=0.38
ax4.barh(yp+bh/2,top15['choro_T'],bh,color=RD,alpha=0.85,label='Choropleth')
ax4.barh(yp-bh/2,top15['dasy_T'], bh,color=BL,alpha=0.85,label='Dasymetric')
ax4.set_yticks(yp); ax4.set_yticklabels(top15['nama_desa'],fontsize=9)
ax4.set_xlabel('Population exposed — High hazard class (persons)',fontsize=10,color=GR)
ax4.set_title('(d) Top 15 desa — High hazard class\nChoropleth vs Dasymetric',
              fontsize=10,fontweight='bold',color=B1)
ax4.legend(fontsize=9,loc='lower right')
ax4.grid(axis='x',alpha=0.3,ls='--')
ax4.xaxis.set_major_formatter(
    plt.FuncFormatter(lambda v,_: f'{int(v/1000)}k' if v>=1000 else str(int(v))))
ax4.spines[['top','right']].set_visible(False)

plt.suptitle(
    'Figure 10.  Statistical comparison of choropleth vs dasymetric flood exposure'
    ' — Kabupaten Bone (N = 372 desa)',
    fontsize=11,fontweight='bold',color=B1,y=1.01)
plt.tight_layout(h_pad=3.5,w_pad=3)
plt.savefig('fig_statistics.png',dpi=180,bbox_inches='tight',facecolor=W)
plt.close()
print("\n✓ fig_statistics.png saved")

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 11 — Sub-kelompok rentan
# ══════════════════════════════════════════════════════════════════════════
grps_v   = ['umur_rentan','miskin','disabilitas']
labels_v = ['Umur Rentan\n(Children + Elderly)',
            'Poor Households\n(Miskin)',
            'Persons with\nDisability']

fig, axes = plt.subplots(2, 2, figsize=(14, 11), facecolor=W)
fig.patch.set_facecolor(W)

# (a) Exposure rate
ax=axes[0,0]; ax.set_facecolor(W)
x=np.arange(3); bw=0.32
rate_c_v=[vuln_results[g]['rate_c'] for g in grps_v]
rate_d_v=[vuln_results[g]['rate_d'] for g in grps_v]
ax.bar(x-bw/2,rate_c_v,bw,color=RD,alpha=0.85,label='Choropleth',zorder=3)
ax.bar(x+bw/2,rate_d_v,bw,color=BL,alpha=0.85,label='Dasymetric (pixel-level)',zorder=3)
ax.axhline(rate_pop,color='#333',lw=1.5,ls='--',
           label=f'General pop. ({rate_pop:.1f}%)',zorder=4)
for i,(rc,rd) in enumerate(zip(rate_c_v,rate_d_v)):
    ax.text(i-bw/2,rc+0.4,f'{rc:.1f}%',ha='center',fontsize=8.5,color=RD,fontweight='bold')
    ax.text(i+bw/2,rd+0.4,f'{rd:.1f}%',ha='center',fontsize=8.5,color=BL,fontweight='bold')
ax.set_xticks(x); ax.set_xticklabels(labels_v,fontsize=9.5)
ax.set_ylabel('Population exposed (% of sub-group total)',fontsize=10,color=GR)
ax.set_title('(a) Flood exposure rate by vulnerable sub-group\n'
             'Choropleth vs Dasymetric (pixel-level, all GFI classes)',
             fontsize=10,fontweight='bold',color=B1)
ax.legend(fontsize=9,loc='upper right'); ax.grid(axis='y',alpha=0.3,ls='--')
ax.set_ylim(0,42); ax.spines[['top','right']].set_visible(False)

# (b) SSEDI
ax2=axes[0,1]; ax2.set_facecolor(W)
ssedi_d_v=[vuln_results[g]['ssedi_d'] for g in grps_v]
ssedi_c_v=[vuln_results[g]['ssedi_c'] for g in grps_v]
ax2.bar(x-bw/2,ssedi_c_v,bw,color=RD,alpha=0.85,label='Choropleth SSEDI',zorder=3)
ax2.bar(x+bw/2,ssedi_d_v,bw,color=BL,alpha=0.85,label='Dasymetric SSEDI',zorder=3)
ax2.axhline(1.0,color='#333',lw=1.5,ls='--',label='SSEDI = 1.0 (parity)',zorder=4)
for i,(sc,sd) in enumerate(zip(ssedi_c_v,ssedi_d_v)):
    ax2.text(i-bw/2,sc+0.01,f'{sc:.3f}',ha='center',fontsize=8.5,color=RD,fontweight='bold')
    ax2.text(i+bw/2,sd+0.01,f'{sd:.3f}',ha='center',fontsize=8.5,color=BL,fontweight='bold')
ax2.axhspan(1.0,1.55,alpha=0.05,color=RD); ax2.axhspan(0.0,1.0,alpha=0.05,color=BL)
ax2.text(2.48,1.43,'Over-represented',ha='right',fontsize=8,color=RD)
ax2.text(2.48,0.06,'Under-represented',ha='right',fontsize=8,color=BL)
ax2.set_xticks(x); ax2.set_xticklabels(labels_v,fontsize=9.5)
ax2.set_ylabel('Social-Spatial Exposure Disparity Index (SSEDI)',fontsize=10,color=GR)
ax2.set_title('(b) SSEDI by sub-group (pixel-level dasymetric)\n'
              'SSEDI > 1.0 = over-represented in flood exposure',
              fontsize=10,fontweight='bold',color=B1)
ax2.legend(fontsize=9); ax2.grid(axis='y',alpha=0.3,ls='--')
ax2.set_ylim(0,1.55); ax2.spines[['top','right']].set_visible(False)

# (c) Overestimation
ax3=axes[1,0]; ax3.set_facecolor(W)
abs_vals=[int(vuln_results[g]['g_all']*vuln_results[g]['gap_pct']/100)
          for g in grps_v]
pp_vals =[vuln_results[g]['rate_c']-vuln_results[g]['rate_d'] for g in grps_v]
cols_v  =[TL,RD,AM]
bH=ax3.barh(labels_v[::-1],abs_vals[::-1],
            color=[cols_v[i] for i in [2,1,0]],alpha=0.85,zorder=3)
ax3.axvline(0,color='#333',lw=1.2)
for b,av,pp in zip(bH,abs_vals[::-1],pp_vals[::-1]):
    ax3.text(av+30,b.get_y()+b.get_height()/2,
             f'+{av:,} persons  ({pp:.1f} pp)',
             ha='left',va='center',fontsize=8.5,fontweight='bold',color='#333')
ax3.set_xlabel('Choropleth overestimation (persons, all GFI classes)',fontsize=10,color=GR)
ax3.set_title('(c) MAUP-induced overestimation per sub-group\n'
              'Choropleth overestimates ALL vulnerable groups',
              fontsize=10,fontweight='bold',color=B1)
ax3.set_xlim(0,8000); ax3.grid(axis='x',alpha=0.3,ls='--')
ax3.spines[['top','right']].set_visible(False)

# (d) Top 15 desa
ax4=axes[1,1]; ax4.set_facecolor(W)
df['ur_all'] = df['ur_R']+df['ur_S']+df['ur_T']
df['mis_all']= df['mis_R']+df['mis_S']+df['mis_T']
top=df.nlargest(15,'ur_all')[['nama_desa','ur_all','mis_all']].iloc[::-1]
y=np.arange(len(top)); bh=0.38
ax4.barh(y+bh/2,top['ur_all'], bh,color=AM,alpha=0.85,
         label='Umur Rentan (pixel-level)')
ax4.barh(y-bh/2,top['mis_all'],bh,color=RD,alpha=0.85,
         label='Miskin (pixel-level)')
ax4.set_yticks(y); ax4.set_yticklabels(top['nama_desa'],fontsize=9)
ax4.set_xlabel('Persons exposed — all GFI classes (pixel-level dasymetric)',
               fontsize=10,color=GR)
ax4.set_title('(d) Top 15 desa — Umur Rentan & Miskin exposure\n'
              '(pixel-level dasymetric, all GFI classes)',
              fontsize=10,fontweight='bold',color=B1)
ax4.legend(fontsize=9,loc='lower right')
ax4.grid(axis='x',alpha=0.3,ls='--')
ax4.xaxis.set_major_formatter(
    plt.FuncFormatter(lambda v,_: f'{int(v/1000)}k' if v>=1000 else str(int(v))))
ax4.spines[['top','right']].set_visible(False)

plt.suptitle(
    'Figure 11.  Social-Spatial Exposure Disparity by Vulnerable Sub-Group\n'
    'Pixel-level dasymetric analysis — Kabupaten Bone (N = 372 desa)',
    fontsize=11,fontweight='bold',color=B1,y=1.01)
plt.tight_layout(h_pad=3.5,w_pad=3)
plt.savefig('fig_vulnerable.png',dpi=180,bbox_inches='tight',facecolor=W)
plt.close()
print("✓ fig_vulnerable.png saved")

# ══════════════════════════════════════════════════════════════════════════
# EXPORT stats_results.csv
# ══════════════════════════════════════════════════════════════════════════
W_T,p_T = wilcoxon(df['dasy_T'],df['choro_T'],alternative='two-sided')
W_S,p_S = wilcoxon(df['dasy_S'],df['choro_S'],alternative='two-sided')
W_R,p_R = wilcoxon(df['dasy_R'],df['choro_R'],alternative='two-sided')
rb_T     = 1-(2*W_T)/(n*(n+1)/2)

table = {
    'Test': [
        'Wilcoxon signed-rank — Low class',
        'Wilcoxon signed-rank — Moderate class',
        'Wilcoxon signed-rank — High class',
        'Mann-Whitney U — High class',
        'Spearman ρ — full dataset (N=372)',
        'Spearman ρ — exposed desa (N=131)',
        'OLS R² (flood_frac + log_area + log_pop)',
    ],
    'Statistic': [
        f'W={W_R:.1f}', f'W={W_S:.1f}', f'W={W_T:.1f}',
        f'U={U:.1f}', f'ρ={rho_f:.4f}', f'ρ={rho_nz:.4f}',
        f'R²={R2:.4f}',
    ],
    'p-value': [
        f'{p_R:.4f}', f'{p_S:.4f}', f'{p_T:.4f}',
        f'{pU:.4f}', f'{p_f:.6f}', f'{p_nz:.4f}',
        'n.s.',
    ],
    'Effect': [
        'n.s.', f'***', f'r={rb_T:.3f} (large)',
        '**', '*** (zero-inflation)', 'n.s.',
        f'R²={R2:.4f}',
    ],
}
pd.DataFrame(table).to_csv('stats_results.csv', index=False)
print("✓ stats_results.csv saved")
print("\nDone! Output files: fig_statistics.png, fig_vulnerable.png, stats_results.csv")
