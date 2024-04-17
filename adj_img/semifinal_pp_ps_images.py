from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.ticker import AutoMinorLocator
import matplotlib.gridspec as gridspec


# input parameters-----------------------------------------------------------------------------------------------------
file_pp = '/home0/cxd170430/NatureGeo/adj_img/data/vz_pp_dt0_final.rsf@'
file_ps = '/home0/cxd170430/NatureGeo/adj_img/data/R_plus_T_final.rsf@'
file_dvs = '/home0/cxd170430/NatureGeo/vs_model/dVs_2d_model/dVs_ross.bin'

nz = 251
dz = 0.04
nx = 751
dx = 0.04
oz = 0.0
ox = 0.0

show = True
#----------------------------------------------------------------------------------------------------------------------


# load data
z = np.linspace(oz, (nz-1)*dz, nz, endpoint=True)
x = np.linspace(ox, (nx-1)*dx, nx, endpoint=True)

ppimg = np.fromfile(file_pp, dtype=np.float32)
ppimg = ppimg.reshape(nz,nx,order='F')
ppimg_norm = ppimg/np.max(abs(ppimg))

psimg = np.fromfile(file_ps, dtype=np.float32)
psimg = psimg.reshape(nz,nx,order='F')
psimg_norm = psimg/np.max(abs(psimg))

interp_x, interp_z = np.mgrid[0:30:151j, 0:10:51j]
reflec_perc = np.fromfile(file_dvs, dtype=np.float32)
reflec_perc = reflec_perc.reshape(51,151,order='F')

#---------------------------------------------

fig = plt.figure(figsize=(10,7.5))
gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1,1], hspace=0.5)


# PP image
ax2 = fig.add_subplot(gs[0,0])
pic1 = ax2.imshow(ppimg_norm, cmap="RdBu", vmin=-0.52, vmax=0.52, extent=(x[0],x[-1],z[-1],z[0]), interpolation='none', aspect='equal')
ax2.set_xlabel('Distance (km)',fontsize=15)
ax2.set_ylabel('Depth (km)',fontsize=15)

# plot contours
CS2 = ax2.contour(interp_x.T, interp_z.T, reflec_perc, [-30, -25, -20], colors='black', linewidths=0.7)
position3 = [(16.8, 8.6), (27.2, 8.2), (1.2, 4.5)]
ax2.clabel(CS2, inline=True, fontsize=13, manual=position3)

# colorbars look
cbar2=fig.colorbar(pic1, ax=ax2, location='right', orientation='vertical', shrink=0.7, anchor=(-0.1,0.5))
cbar2.set_label('PP reflectivity',fontsize=13)
cbar2.set_ticks([-0.5,0,0.5])
cbar2.ax.tick_params(labelsize=13, length=5, width=1.0)
cbar2.ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
cbar2.ax.tick_params(which='minor', length=5, width=1.0)
cbar2.outline.set_linewidth(1.8)
# axis-labels to mark
ax2.set_yticks([0,2,4,6,8,10])
ax2.xaxis.set_minor_locator(tck.AutoMinorLocator(5))
ax2.yaxis.set_minor_locator(tck.AutoMinorLocator(4))
# axis-labels look
ax2.tick_params(axis='x', labelsize=15, length=7, top=True)
ax2.tick_params(axis='y', labelsize=15, length=7, right=True)
# ticks look
ax2.tick_params(which='minor', length=3.8, width=1.0, top=True, right=True)
ax2.tick_params(width=1.0)
# frame look
plt.setp(ax2.spines.values(), linewidth=2.0)
# annotations
ax2.annotate('', xy=(0,0), xycoords='data', xytext=(0,-1.5), textcoords='data', size=13,
             arrowprops=dict(arrowstyle="-",facecolor='black',linewidth=1.4))
ax2.annotate('', xy=(30,0), xycoords='data', xytext=(30,-1.5), textcoords='data', size=13,
             arrowprops=dict(arrowstyle="-",facecolor='black',linewidth=1.4))
ax2.annotate('', xy=(12.3,0.05), xycoords='data', xytext=(12.3,-1.5), textcoords='data', size=13,
             arrowprops=dict(arrowstyle="-",facecolor='black',linewidth=1.4))
ax2.text(9.8, -1.7, 'Mud Volcano', fontsize=15)
ax2.text(-0.35, -1.7, 'X', fontsize=15)
ax2.text(29.65, -1.7, 'X$^{\prime}$', fontsize=15)
ax2.text(-0.13, 1.2, 'A', transform=ax2.transAxes, fontname='DejaVu Sans', size=19, weight='bold')


# PS image
ax3 = fig.add_subplot(gs[1,0], sharex=ax2)
pic2 = ax3.imshow(psimg_norm, cmap="RdBu", vmin=-0.42, vmax=0.42, extent=(x[0],x[-1],z[-1],z[0]), interpolation='none', aspect='equal')
ax3.set_xlabel('Distance (km)',fontsize=15)
ax3.set_ylabel('Depth (km)',fontsize=15)

# plot contours
CS3 = ax3.contour(interp_x.T, interp_z.T, reflec_perc, [-30, -25, -20], colors='black', linewidths=0.7)
position3 = [(16.8, 8.6), (27.2, 8.2), (1.2, 4.5)]
ax3.clabel(CS3, inline=True, fontsize=13, manual=position3)

# colorbars look
cbar3=fig.colorbar(pic2, ax=ax3, location='right', orientation='vertical', shrink=0.7, anchor=(-0.1,0.5))
cbar3.set_label('PS reflectivity',fontsize=13)
cbar3.set_ticks([-0.4,0,0.4])
cbar3.ax.tick_params(labelsize=13, length=5, width=1.0)
cbar3.ax.yaxis.set_minor_locator(AutoMinorLocator(n=4))
cbar3.ax.tick_params(which='minor', length=5, width=1.0)
cbar3.outline.set_linewidth(1.8)
# axis-labels to mark
ax3.set_yticks([0,2,4,6,8,10])
ax3.xaxis.set_minor_locator(tck.AutoMinorLocator(5))
ax3.yaxis.set_minor_locator(tck.AutoMinorLocator(4))
# axis-labels look
ax3.tick_params(axis='x', labelsize=15, length=7, top=True)
ax3.tick_params(axis='y', labelsize=15, length=7, right=True)
# ticks look
ax3.tick_params(which='minor', length=3.8, width=1.0, top=True, right=True)
ax3.tick_params(width=1.2)
# frame look
plt.setp(ax3.spines.values(), linewidth=2.0)
# annotations
ax3.annotate('', xy=(0,0), xycoords='data', xytext=(0,-1.5), textcoords='data', size=13,
             arrowprops=dict(arrowstyle="-",facecolor='black',linewidth=1.4))
ax3.annotate('', xy=(30,0), xycoords='data', xytext=(30,-1.5), textcoords='data', size=13,
             arrowprops=dict(arrowstyle="-",facecolor='black',linewidth=1.4))
ax3.annotate('', xy=(12.3,0.05), xycoords='data', xytext=(12.3,-1.5), textcoords='data', size=13,
             arrowprops=dict(arrowstyle="-",facecolor='black',linewidth=1.4))
ax3.text(9.8, -1.7, 'Mud Volcano', fontsize=15)
ax3.text(-0.35, -1.7, 'X', fontsize=15)
ax3.text(29.65, -1.7, 'X$^{\prime}$', fontsize=15)
ax3.text(-0.13, 1.2, 'B', transform=ax3.transAxes, fontname='DejaVu Sans', size=19, weight='bold')



if show:
    plt.show()
else:
    fig.savefig('./figs/pp_ps_images.png', dpi=300, facecolor='w', bbox_inches='tight')

