import matplotlib.pyplot as plt
import numpy as np
from GMI_Finesse_Funcs import *
finesse.configure(plotting=True)
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.style.use('bmh')
plt.rcParams.update({'font.family': 'times'})
plt.rcParams['mathtext.fontset'] = 'cm'  # Use Computer Modern font
plt.rcParams['font.size'] = 8  # Set the base font size
from matplotlib import pyplot as plt, ticker as mticker

gmipower = 125
ligopower =  gmipower
finesse.configure(plotting=True)
pend = 1
nsr = 1
rad_pressure = 1
phi0= 0
phi3 = 0
losses = [0,1]
rp = 'with' if rad_pressure else 'without'
fig, ax = plt.subplots(2,1,figsize=(3,3))

for ind,i in enumerate(losses):
    lossy = i
    phi2 = np.array([89.9,89.97,89.99])
    phi1 = 90-phi2
    # ligo_ = ligo(power=ligopower,lossy=lossy,pendula=pend,nsr=nsr)
    gmi0 = gmi(power=gmipower,phi2=phi2[0],phi1=phi1[0],phi0=phi0,phi3=phi3,lossy=lossy,pendula=pend,nsr=nsr)
    gmi1 = gmi(power=gmipower,phi2=phi2[1],phi1=phi1[1],phi0=phi0,phi3=phi3,lossy=lossy,pendula=pend,nsr=nsr)
    gmi2 = gmi(power=gmipower,phi2=phi2[2],phi1=phi1[2],phi0=phi0,phi3=phi3,lossy=lossy,pendula=pend,nsr=nsr)

    # ligo_model = ligo_.run("parallel(xaxis(fsig.f, log, 1, 1M, 1000,name='default'),xaxis(phi1,lin,0,.2,1000,name='trans'))")
    # mi_ = mi(power=gmipower,lossy=lossy,nsr=nsr)
    # mi_model = mi_.run("xaxis(fsig.f, log, 1, 1M, 1000)")

    out0 = gmi0.run("parallel(xaxis(fsig.f, log, 1, 1M, 1000,name='default'),xaxis(phi1,lin,0,.2,1000,name='trans'))")
    out1 = gmi1.run("parallel(xaxis(fsig.f, log, 1, 1M, 1000,name='default'),xaxis(phi1,lin,0,.2,1000,name='trans'))")
    out2 = gmi2.run("parallel(xaxis(fsig.f, log, 1, 1M, 1000,name='default'),xaxis(phi1,lin,0,.2,1000,name='trans'))")

    phi = np.linspace(0,.4,1001)
    label0 = 'aLIGO' if ind == 1 else None
    label1 = f'$\phi_N={2*phi1[0]:.2f}$°' if ind == 1 else None
    label2 = f'$\phi_N={2*phi1[1]:.2f}$°' if ind==1 else None
    label3 = f'$\phi_N={2*phi1[2]:.2f}$°' if ind==1 else None
    # ax[ind].plot(phi,ligo_model['trans']['output']/gmipower,label=label0,color='k')
    ax[ind].plot(phi,out0['trans']['output']/gmipower,label=label1,color='r')
    ax[ind].plot(phi,out1['trans']['output']/gmipower,label=label2,color='g')
    ax[ind].plot(phi,out2['trans']['output']/gmipower,label=label3,color='b')
    # ax.plot(phi,mi_model[f'NSR_{rp}_RP'],label='Michelson',color='purple')
    # ax[ind].set_ylim(1e-25,1e-17)
    ax[ind].semilogy()
    # ax[ind].set_ylim(0,1)
    ax[ind].set_ylabel('Transmission', fontsize=11, font="Times New Roman")
    ax[ind].yaxis.set_major_locator(mticker.LogLocator(numticks=5))
    ax[ind].yaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
    ax[ind].set_ylim(1e-8, 1)
ax[1].set_xlabel('$\phi_E$ $(^\circ)$',fontsize=11,font='Times New Roman')
# ax[1].set_ylabel('Strain Sensitivity [/$\sqrt{Hz}$]',fontsize=11)
ax[0].set_title('a.) Lossless Mirrors',fontsize=9,font='Times New Roman')
ax[1].set_title(r'b.) T=15x10$^{-6}$, L=5x10$^{-6}$',fontsize=9,font='Times New Roman')

# ax[1].text(s="x10",x=0.062,y=0.008,fontsize='7')
# ax[1].text(s=r"x10$^3$",x=0.015,y=0.011,fontsize='7')
# plt.subplots_adjust(wspace=0.4)
# ax[1].legend(ncols=1,bbox_to_anchor=(1.05,.5),loc='center left')
# fig.legend(ncols=2,loc='outside lower center')
fig.legend(ncols=3,bbox_to_anchor=(.5, 0.05),loc='upper center',
          bbox_transform=fig.transFigure,prop=FontProperties(family='times',size=7))
plt.tight_layout()
plt.savefig('Figure4.png',bbox_inches='tight',dpi=400)
plt.show()
