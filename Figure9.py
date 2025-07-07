from matplotlib import pyplot as plt, ticker as mticker
import numpy as np
from GMI_Finesse_Funcs import *
finesse.configure(plotting=True)
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.style.use('bmh')
plt.rcParams.update({'font.family': 'times'})

gmipower = 125
ligopower =  gmipower
finesse.configure(plotting=True)
pend = 1
nsr = 1
rad_pressure = 1
phi0= 0
phi3 = 0
dphi = 0
losses = [0,1]
rp = 'with' if rad_pressure else 'without'
fig, ax = plt.subplot_mosaic([[0],[1]],figsize=(3,3))

for ind,i in enumerate(losses):
    lossy = i
    phi1 = np.array([.1,.03,.01])
    phi2 = 90-phi1-.00
    ligo_ = ligo(power=ligopower,lossy=lossy,dphi=dphi,pendula=pend,nsr=nsr)
    gmi0 = gmi(power=gmipower,phi2=phi2[0],phi1=phi1[0],dphi=dphi,phi0=phi0,phi3=phi3,lossy=lossy,pendula=pend,nsr=nsr)
    gmi1 = gmi(power=gmipower,phi2=phi2[1],phi1=phi1[1],dphi=dphi,phi0=phi0,phi3=phi3,lossy=lossy,pendula=pend,nsr=nsr)
    gmi2 = gmi(power=gmipower,phi2=phi2[2],phi1=phi1[2],dphi=dphi,phi0=phi0,phi3=phi3,lossy=lossy,pendula=pend,nsr=nsr)

    ligo_model = ligo_.run("xaxis(fsig.f, log, 1, 1M, 1000)")
    # mi_ = mi(power=gmipower,lossy=lossy,nsr=nsr)
    # mi_model = mi_.run("xaxis(fsig.f, log, 1, 1M, 1000)")

    out0 = gmi0.run("parallel(xaxis(fsig.f, log, 1, 1M, 1000,name='default'),xaxis(phi1,lin,0,1,500,name='trans'))")
    out1 = gmi1.run("parallel(xaxis(fsig.f, log, 1, 1M, 1000,name='default'),xaxis(phi1,lin,0,1,500,name='trans'))")
    out2 = gmi2.run("parallel(xaxis(fsig.f, log, 1, 1M, 1000,name='default'),xaxis(phi1,lin,0,1,500,name='trans'))")

    phi = np.logspace(0,6,1001)
    label0 = 'aLIGO' if ind == 1 else None
    label1 = f'GMI, $\phi_N={2 * phi1[0]:.2f}$°' if ind == 1 else None
    label2 = f'GMI, $\phi_N={2 * phi1[1]:.2f}$°' if ind == 1 else None
    label3 = f'GMI, $\phi_N={2 * phi1[2]:.2f}$°' if ind == 1 else None
    ax[ind].plot(phi,ligo_model[f'NSR_{rp}_RP'],label=label0,color='k')
    ax[ind].plot(phi,out0['default'][f'NSR_{rp}_RP'],label=label1,color='r')
    ax[ind].plot(phi,out1['default'][f'NSR_{rp}_RP'],label=label2,color='g')
    ax[ind].plot(phi,out2['default'][f'NSR_{rp}_RP'],label=label3,color='b')
    # ax.plot(phi,mi_model[f'NSR_{rp}_RP'],label='Michelson',color='purple')
    ax[ind].set_ylim(1e-25,1e-15)
    ax[ind].set_xlim(1, 1e4)
    ax[ind].loglog()
    ax[ind].xaxis.set_major_locator(mticker.LogLocator(numticks=999))
    ax[ind].xaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
    ax[ind].yaxis.set_major_locator(mticker.LogLocator(numticks=6))
    ax[ind].yaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
ax[1].set_xlabel('Signal Frequency [Hz]',fontsize=10,font='Times New Roman')
# ax[1].set_ylabel('Strain Sensitivity [/$\sqrt{Hz}$]',fontsize=11)
ax[0].set_title('Common-Mode Detection\na.) Lossless Mirrors',fontsize=9,font='Times New Roman')
ax[1].set_title('b.) T=15x10$^{-6}$, L=5x10$^{-6}$',fontsize=9,font='Times New Roman')

# plt.subplots_adjust(wspace=0.4)
# ax[1].legend(ncols=1,bbox_to_anchor=(1.05,.5),loc='center left')
# fig.legend(ncols=2,loc='outside lower center')
fig.supylabel(r'Strain Sensitivity $\left[\text{h}/\sqrt{\text{Hz}}\right]$',fontsize=10,font="Times New Roman",x=0)
fig.legend(ncols=2,bbox_to_anchor=(.5, 0.02),loc='upper center',
          bbox_transform=fig.transFigure,prop=FontProperties(family='times',size=8))
plt.tight_layout(pad=.1)
plt.savefig('Figure9.png',bbox_inches='tight',dpi=400)
plt.show()
