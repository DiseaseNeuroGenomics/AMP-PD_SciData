#!/usr/bin/env python
# coding: utf-8


import pandas as pd, matplotlib.pyplot as plt, numpy as np, os
import seaborn as sns
from scipy.stats import spearmanr, pearsonr


plt.style.use('tableau-colorblind10')
color = {'donor' + str(k): '' for k in range(1,7)}
for c, k in zip(plt.rcParams['axes.prop_cycle'].by_key()['color'], color.keys()):
    color[k] = c


wet_lab_dir = "wet_lab/"
dry_lab_dir = "dry_lab/"
output_dir = "output_dir/"


# <h3>Plot for Fig 3a</h3>
pre_proc_df = pd.read_csv(os.path.join(dry_lab_dir, 'SciData_fig3a.csv'))

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 5))
sns.barplot(
   data=pre_proc_df.loc[pre_proc_df['status'] == 'kept',], x="Pool_name", y="n_cells",
    ax=ax, color='lightblue'
)
plt.ticklabel_format(style='plain', axis='y')
plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.gca().axes.get_xaxis().set_ticks([]) # Hide x-axis
ml =  pre_proc_df.loc[pre_proc_df['status'] == 'kept', 'n_cells'].mean()
ax.axhline(ml, ls='--', c='r', lw=1)
ax.set_xlabel('Pools')
ax.set_ylabel('Cell Counts')
plt.text(8, ml, 'mean', fontsize=8, va='center', ha='center', backgroundcolor='w')
plt.savefig(os.path.join(output_dir, "SciData_fig3a.pdf"), format='pdf')


# <h1>For Fig 3b</h1>
fig3b_inp = pd.read_csv(os.path.join(dry_lab_dir, "SciData_fig3b.csv"))
fig3b_inp['demux_dons'] = fig3b_inp['demux_dons'].astype('category')
fig3b_inp['final_dons'] = fig3b_inp['final_dons'].astype('category')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 5))
bplot1 = fig3b_inp.boxplot(column=['vir_res'], by='demux_dons', grid=False, rot=0, fontsize=10, notch=True, 
                        ax=ax, patch_artist=True, return_type='dict', 
                        positions = np.arange(len(fig3b_inp['demux_dons'].unique()))-0.25, widths=0.4,)
bplot2 = fig3b_inp.boxplot(column=['cell_counts'], by='final_dons', grid=False, rot=0, fontsize=10, notch=True, 
                        ax=ax, patch_artist=True, return_type='dict', 
                        positions = np.arange(len(fig3b_inp['final_dons'].unique()))+0.25, widths=0.4,)
fig.suptitle('Cell counts distribution')
for b in [bplot1, bplot2]:
    for k in b.keys():
        for i, (p, c) in enumerate(zip(b[k]['boxes'], color.values())):
            p.set_facecolor(c)
            p.set_alpha(0.5)
            if k == 'vir_res':
                p.set_linestyle('--')
                p.set_linewidth(3)
                if i == 0:
                    plt.plot([], color='black', linestyle='--', linewidth=3, label='Before QC')
                    plt.legend()
            else:
                p.set_linestyle('-')
                p.set_linewidth(1)
                if i == 1:
                    plt.plot([], color='black', linestyle='-', linewidth=1, label='After QC')
                    plt.legend()
                    

ax.set_xlabel(None)
ax.set_ylabel("cell counts")
ax.yaxis.grid(True)
ax.text(0.5,-0.1, "Highest --------------------------------------------------------------------> Lowest", size=12, ha="center", 
         transform=ax.transAxes, fontweight='bold')
ax.set_title('')
ax.set_xticks([i for i in np.arange(len(fig3b_inp['demux_dons'].unique()))])
ax.set_xticklabels(fig3b_inp['demux_dons'].cat.categories)
plt.savefig(os.path.join(output_dir, 'SciData_fig3b_1.pdf'), format='pdf')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 5))
bplot1 = fig3b_inp.boxplot(column=['init_ratio'], by='demux_dons', grid=False, rot=0, fontsize=10, notch=True, 
                        ax=ax, patch_artist=True, return_type='dict', positions = np.arange(len(fig3b_inp['demux_dons'].unique()))-0.25, widths=0.4,)
bplot2 = fig3b_inp.boxplot(column=['final_ratio'], by='final_dons', grid=False, rot=0, fontsize=10, notch=True, 
                        ax=ax, patch_artist=True, return_type='dict', positions = np.arange(len(fig3b_inp['final_dons'].unique()))+0.25, widths=0.4,)
fig.suptitle('Cell counts distribution')
for b in [bplot1, bplot2]:
    for k in b.keys():
        for i, (p, c) in enumerate(zip(b[k]['boxes'], color.values())):
            p.set_facecolor(c)
            p.set_alpha(0.5)
            if k == 'init_ratio':
                p.set_linestyle('--')
                p.set_linewidth(3)
                if i == 0:
                    plt.plot([], color='black', linestyle='--', linewidth=3, label='Before QC')
                    plt.legend()
            else:
                p.set_linestyle('-')
                p.set_linewidth(1)
                if i == 1:
                    plt.plot([], color='black', linestyle='-', linewidth=1, label='After QC')
                    plt.legend()

            
ax.set_ylabel("Ratio")
ax.yaxis.grid(True)
ax.text(0.5,-0.1, "Highest --------------------------------------------------------------------> Lowest", size=12, ha="center", 
         transform=ax.transAxes, fontweight='bold')
ax.set_xlabel(None)
ax.set_title('')
ax.set_xticks([i for i in np.arange(len(fig3b_inp['demux_dons'].unique()))])
ax.set_xticklabels(fig3b_inp['demux_dons'].cat.categories)

plt.savefig(os.path.join(output_dir, 'SciData_fig3b_2.pdf'), format='pdf')


# <h2>Plot for Fig 3c</h2>
fig3c_inp = pd.read_csv(os.path.join(dry_lab_dir, "SciData_fig3c.csv"))

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 5))
sns.scatterplot(
   data=fig3c_inp,
    x="n_cells_high_rep", y="n_cells_low_rep",
    ax=ax,
)
spear_corr = spearmanr(fig3c_inp['n_cells_high_rep'], fig3c_inp['n_cells_low_rep'])
pear_corr = pearsonr(fig3c_inp['n_cells_high_rep'], fig3c_inp['n_cells_low_rep'])
min_val = max(fig3c_inp['n_cells_high_rep'].min(), fig3c_inp['n_cells_low_rep'].min())-1000
max_val = max(fig3c_inp['n_cells_high_rep'].max(), fig3c_inp['n_cells_low_rep'].max())+1000
x_vals = np.linspace(min_val, max_val, 50)
sns.lineplot(
    x=x_vals,
    y=x_vals,
    ax=ax,
    lw=1, ls='--',
    color='red'
)
plt.ticklabel_format(style='plain', axis='both')
ax.text(15000, 17000, 'y=x', rotation=45, fontsize=8, va='center', ha='center', backgroundcolor='w', rotation_mode='anchor',
              transform_rotates_text=True, color='red')
sns.regplot(data=fig3c_inp, x="n_cells_high_rep", y="n_cells_low_rep", ax=ax, 
            line_kws=dict(color="green", lw=1, ls='-'),
           )

textstr = '\n\n'.join([
    'spearman: '+str(round(spear_corr.correlation,2)) + '\np=' + str(np.format_float_scientific(spear_corr.pvalue, precision=2)),
    'pearson: ' + str(round(pear_corr.statistic,2)) + '\np=' + str(np.format_float_scientific(pear_corr.pvalue, precision=2)),
])

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

# place a text box in upper left in axes coords
ax.text(15000, 30000, textstr, fontsize=9,
        verticalalignment='top', bbox=props, va='center', ha='center', backgroundcolor='w')
ax.set_xlabel('Rep (Higher cell count)')
ax.set_ylabel('Rep (Lower cell count)')
plt.gca().ticklabel_format(axis='both', style='sci', scilimits=(0, 0))
plt.savefig(os.path.join(output_dir, 'SciData_fig3c.pdf'), format='pdf')
