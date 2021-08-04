import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


# Plot style

plt.rcParams["font.family"] = "Times New Roman"
plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)

fsize = 18
tsize = 15
tdir = 'in'
major = 5.0
minor = 3.0
lwidth = 1.8
lhandle = 2.0

plt.style.use('default')
plt.rcParams['font.size'] = fsize
plt.rcParams['legend.fontsize'] = tsize
plt.rcParams['text.usetex']=True
plt.rcParams['xtick.minor.visible'], plt.rcParams['xtick.top'] = True,True
plt.rcParams['ytick.minor.visible'], plt.rcParams['ytick.right'] = True,True
plt.rcParams['xtick.direction'], plt.rcParams['ytick.direction'] = 'in','in' 

left = b1 = 0.15
width = 0.75
h1 = 0.75
height = 0.75

n, mb, mbe, ms, mse, ssfr, sfc, eta, etae, vc = np.loadtxt('MassLoading.dat',unpack=True,dtype='U24,f,f,f,f,f,f,f,f,f')

eta_HI = eta[0:11]
etae_HI = etae[0:11]
lms_HI = ms[0:11]
lmse_HI = mse[0:11]
vc_HI = vc[0:11]
ssfr_HI= ssfr[0:11]

# McQuinn eta errors from wind assumptions, giving

etae_HI = etae_HI * 0.0 + (12.5 / 37.5) * eta_HI

eta_UV = eta[11:]
etae_UV = etae[11:]
lms_UV = ms[11:]
lmse_UV = mse[11:]
vc_UV = vc[11:]
ssfr_UV = ssfr[11:]

# Ciampa et al. 2021 - from range for IVC and HVC
eta_LMC = [4.545,999]
etae_LMC =[1.011,999]
lms_LMC = [np.log10(3e9),999]
ssfr_LMC = [0.32/3e9,999]

# Comparison with FIRE 
lms_sim=np.arange(6,11,0.1)
eta_sim = 3.6 * (10 ** lms_sim / 1.0e10) ** (-0.35)

fig=plt.figure(figsize=(8,6))
ax=fig.add_axes([left,b1,width,height])
ax.set_xlim(6.5,10.99)
ax.set_ylim(0.01,780.9)
ax.set_ylabel('$\dot{M_o}$/SFR', labelpad = 10)
ax.set_xlabel('log($M_*/M_\odot)$', labelpad = 10,)
ax.set_yscale('log')

ax.errorbar(lms_UV, eta_UV, yerr = etae_UV, xerr = lmse_UV , fmt = 'None',  ecolor = 'k',elinewidth = 0.5,zorder=0)
ax.errorbar(lms_HI, eta_HI, yerr = etae_HI, xerr = lmse_HI , fmt = 'None',  ecolor = 'k',elinewidth = 0.5,zorder=0)
ax.errorbar(lms_LMC, eta_LMC, yerr = etae_LMC, fmt = 'None',  ecolor = 'k',elinewidth = 0.5,zorder=0)

hi = ax.scatter(lms_HI,eta_HI,marker='o',s=200,c=ssfr_HI,cmap='PuRd',vmax=20,vmin=0.05, label = 'Ha, McQuinn et al. 2019')
uv = ax.scatter(lms_UV,eta_UV,marker='s',s=120,c=ssfr_UV,cmap='PuRd',vmax=20,vmin=0.05, label = 'UV, Chisholm et al. 2017')
lmc = ax.scatter(lms_LMC,eta_LMC,marker='d',s=200,c=ssfr_LMC,
                 cmap='PuRd',vmax=20,vmin=0.05,edgecolor='k', label = 'LMC Ha, Ciampa et al. 2021')
ax.plot(lms_sim,eta_sim, c='k', alpha=0.5,lw=2,zorder=0,label='FIRE, Muratov et al. 2015')

ax.legend(loc='upper right')
leg = ax.get_legend()
leg.legendHandles[1].set_color('None')
leg.legendHandles[1].set_edgecolor('k')
leg.legendHandles[2].set_color('None')
leg.legendHandles[2].set_edgecolor('k')
leg.legendHandles[3].set_color('None')
leg.legendHandles[3].set_edgecolor('k')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(hi, cax=cax, label='sSFR ($10^{-10}$ yr$^{-1}$)')

plt.savefig('MassLoading.pdf')
