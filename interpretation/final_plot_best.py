import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.ticker import AutoMinorLocator
import matplotlib.gridspec as gridspec
import matplotlib as mpl


# read obs
file_obs = '/home0/cxd170430/NatureGeo/interpretation/plot_data/obs.txt'
obs = np.loadtxt(file_obs)
obs_offset = obs[:,0]/1000.0
obs_pspp = obs[:,1]
obs_error = obs[:,2]

# read syn fluid
file_syn_f1 = '/home0/cxd170430/NatureGeo/interpretation/plot_data/f0.08.txt'
syn_f1 = np.loadtxt(file_syn_f1)
syn_f1_offset = syn_f1[:,0]/1000.0
syn_f1_pspp = syn_f1[:,1]

file_syn_f2 = '/home0/cxd170430/NatureGeo/interpretation/plot_data/f0.12.txt'
syn_f2 = np.loadtxt(file_syn_f2)
syn_f2_offset = syn_f2[:,0]/1000.0
syn_f2_pspp = syn_f2[:,1]

file_syn_f3 = '/home0/cxd170430/NatureGeo/interpretation/plot_data/f0.16.txt'
syn_f3 = np.loadtxt(file_syn_f3)
syn_f3_offset = syn_f3[:,0]/1000.0
syn_f3_pspp = syn_f3[:,1]

# read syn melt
file_syn_m1 = '/home0/cxd170430/NatureGeo/interpretation/plot_data/m0.01.txt'
syn_m1 = np.loadtxt(file_syn_m1)
syn_m1_offset = syn_m1[:,0]/1000.0
syn_m1_pspp = syn_m1[:,1]

file_syn_m2 = '/home0/cxd170430/NatureGeo/interpretation/plot_data/m0.08.txt'
syn_m2 = np.loadtxt(file_syn_m2)
syn_m2_offset = syn_m2[:,0]/1000.0
syn_m2_pspp = syn_m2[:,1]

file_syn_m3 = '/home0/cxd170430/NatureGeo/interpretation/plot_data/m0.16.txt'
syn_m3 = np.loadtxt(file_syn_m3)
syn_m3_offset = syn_m3[:,0]/1000.0
syn_m3_pspp = syn_m3[:,1]

# read fluid misfit
file_misfit_f = '/home0/cxd170430/NatureGeo/interpretation/plot_data/misfit_fluid.txt'
misfit_f = np.loadtxt(file_misfit_f)
pf = misfit_f[:-1,0]
chif = misfit_f[:-1,5]

# read melt misfit
file_misfit_m = '/home0/cxd170430/NatureGeo/interpretation/plot_data/misfit_melt.txt'
misfit_m = np.loadtxt(file_misfit_m)
pm = misfit_m[:-1,0]
chim = misfit_m[:-1,5]

#---------------------------------
fig = plt.figure(figsize=(13,6))
gs = fig.add_gridspec(nrows=2, ncols=2, width_ratios=[1,1], height_ratios=[1,1], wspace=0.21, hspace=0.07)

# ps/pp ratio vs. offset
ax1 = fig.add_subplot(gs[:,0])

ax1.errorbar(obs_offset, obs_pspp, yerr=obs_error, capsize=3.5, marker='o', markersize=3.7, color='black', linestyle='none')

ax1.plot(syn_f1_offset, syn_f1_pspp, color='tab:red', lw=1.25, label='0.08')
ax1.plot(syn_f2_offset, syn_f2_pspp, color='tab:blue', lw=1.25, label='0.12')
ax1.plot(syn_f3_offset, syn_f3_pspp, color='tab:orange', lw=1.25, label='0.16')
ax1.fill_between(syn_f1_offset, syn_f1_pspp, syn_f3_pspp, color='mistyrose', alpha=0.3)

ax1.plot(syn_m1_offset, syn_m1_pspp, color='tab:red', linestyle='dashed', lw=1.25, label='0.01')
ax1.plot(syn_m2_offset, syn_m2_pspp, color='tab:blue', linestyle='dashed', lw=1.25, label='0.08')
ax1.plot(syn_m3_offset, syn_m3_pspp, color='tab:orange', linestyle='dashed', lw=1.25, label='0.16')
#ax1.fill_between(syn_m1_offset, syn_m1_pspp, syn_m3_pspp, color='linen', alpha=0.3)
ax1.set_xlabel('Offset (km)',fontsize=14)
ax1.set_ylabel('PzS / PzP ratio',fontsize=14)

ax1.set_xlim(2.13, 5.07)
ax1.set_ylim(0.56, 1.91)
ax1.set_xticks([2.2, 2.6, 3, 3.4, 3.8, 4.2, 4.6, 5])
ax1.set_yticks([0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8])
ax1.xaxis.set_minor_locator(tck.AutoMinorLocator(4))
ax1.yaxis.set_minor_locator(tck.AutoMinorLocator(4))

ax1.tick_params(axis='x', labelsize=14, length=7)
ax1.tick_params(axis='y', labelsize=14, length=7)
# ticks look
ax1.tick_params(which='minor', length=3.8, width=1.2)
ax1.tick_params(width=1.2)
# frame look
plt.setp(ax1.spines.values(), linewidth=1.8)

ax1.legend(fontsize=11, title='Porosity',title_fontsize=11)
ax1.text(4.2,0.78,'scH$_{2}$O', fontsize=14, rotation=5, rotation_mode='anchor')
ax1.text(3.8,1.162,'Rhyolite melt', fontsize=14, rotation=15, rotation_mode='anchor')
ax1.text(-0.15, 1.0, 'A', transform=ax1.transAxes, fontname='DejaVu Sans', size=18, weight='bold')

# misfit vs. porosity
ax2 = fig.add_subplot(gs[0,1])
ax2.plot(pm, chim, color='tab:red', linestyle='dashed', linewidth=1.25, marker='s', markersize=6, label='Rhyolite melt')

ax2.set_xlim(0, 0.26)
ax2.set_ylim(85, 820)
ax2.set_yticks([100, 275, 450, 625, 800])
ax2.yaxis.set_minor_locator(tck.AutoMinorLocator(3))
ax2.tick_params(axis='x', labelsize=14, length=7)
ax2.tick_params(axis='y', labelsize=14, length=7)
# ticks look
ax2.tick_params(which='minor', length=3.8, width=1.2)
ax2.tick_params(width=1.2)
# frame look
plt.setp(ax2.spines.values(), linewidth=1.8)

ax2.text(-0.155, 1.0, 'B', transform=ax2.transAxes, fontname='DejaVu Sans', size=18, weight='bold')



ax3 = fig.add_subplot(gs[1,1])
ax3.plot(pf, chif, color='tab:blue', linewidth=1.25, marker='o', markersize=6, label='scH$_{2}$O') # markeredgecolor=None, markerfacecolor=None

ax3.set_xlabel('Porosity',fontsize=14)
ax3.set_ylabel('$\chi^{2}$ misfit',fontsize=14)
ax3.yaxis.set_label_coords(-0.116, 1.05)
ax3.set_xlim(0, 0.26)
ax3.set_ylim(-2.5, 65)
ax3.set_yticks([0, 15, 30, 45, 60])
ax3.xaxis.set_minor_locator(tck.AutoMinorLocator(5))
ax3.yaxis.set_minor_locator(tck.AutoMinorLocator(3))
ax3.tick_params(axis='x', labelsize=14, length=7)
ax3.tick_params(axis='y', labelsize=14, length=7)
# ticks look
ax3.tick_params(which='minor', length=3.8, width=1.2)
ax3.tick_params(width=1.2)
# frame look
plt.setp(ax3.spines.values(), linewidth=1.8)

ax2.spines.bottom.set_visible(False)
ax3.spines.top.set_visible(False)
ax2.xaxis.tick_top()
ax2.tick_params(labeltop=False)
ax3.xaxis.tick_bottom()
ax2.tick_params(axis='x',top=False)

h1, l1 = ax2.get_legend_handles_labels()
h2, l2 = ax3.get_legend_handles_labels()
ax2.legend(h1+h2, l1+l2, fontsize=13)

d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax2.plot([0, 1], [0, 0], transform=ax2.transAxes, **kwargs)
ax3.plot([0, 1], [1, 1], transform=ax3.transAxes, **kwargs)

fig.savefig('/home0/cxd170430/NatureGeo/interpretation/test.png', dpi=300, facecolor='w', bbox_inches='tight')
