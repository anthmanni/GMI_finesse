from matplotlib import pyplot as plt, ticker as mticker
import numpy as np
from GMI_Finesse_Funcs import *
finesse.configure(plotting=True)
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.style.use('bmh')
plt.rcParams.update({'font.family': 'times'})
# plt.rcParams['font.size']=8

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
fig, ax = plt.subplots(2,2,figsize=(5,4))
phi = np.logspace(0,6,1001)

for ind,i in enumerate(losses):
    lossy = i
    phi1 = np.array([.1,.03,.01])
    ligo_ = ligo(power=ligopower,lossy=lossy,pendula=pend,nsr=nsr)
    ligo_model = ligo_.run("xaxis(fsig.f, log, 1, 1M, 1000)")
    for ind2,j in enumerate([0,0.0005]):
        phi2 = 90 - phi1 - j
        gmi0 = gmi(power=gmipower,phi2=phi2[0],phi1=phi1[0],phi0=phi0,phi3=phi3,lossy=lossy,pendula=pend,nsr=nsr)
        gmi1 = gmi(power=gmipower,phi2=phi2[1],phi1=phi1[1],phi0=phi0,phi3=phi3,lossy=lossy,pendula=pend,nsr=nsr)
        gmi2 = gmi(power=gmipower,phi2=phi2[2],phi1=phi1[2],phi0=phi0,phi3=phi3,lossy=lossy,pendula=pend,nsr=nsr)

        # mi_ = mi(power=gmipower,lossy=lossy,nsr=nsr)
        # mi_model = mi_.run("xaxis(fsig.f, log, 1, 1M, 1000)")

        out0 = gmi0.run("parallel(xaxis(fsig.f, log, 1, 1M, 1000,name='default'),xaxis(phi1,lin,0,1,500,name='trans'))")
        out1 = gmi1.run("parallel(xaxis(fsig.f, log, 1, 1M, 1000,name='default'),xaxis(phi1,lin,0,1,500,name='trans'))")
        out2 = gmi2.run("parallel(xaxis(fsig.f, log, 1, 1M, 1000,name='default'),xaxis(phi1,lin,0,1,500,name='trans'))")

        label0 = 'aLIGO' if ind == 1 & ind2 == 1 else None
        label1 = f'GMI, $\phi_N={2*phi1[0]:.2f}$°' if ind == 1 & ind2 == 1 else None
        label2 = f'GMI, $\phi_N={2*phi1[1]:.2f}$°' if ind==1  & ind2 == 1 else None
        label3 = f'GMI, $\phi_N={2*phi1[2]:.2f}$°' if ind==1  & ind2 == 1 else None
        print(label0,label1,label2,label3)
        ax[ind,ind2].plot(phi,ligo_model[f'NSR_{rp}_RP'],label=label0,color='k')
        ax[ind,ind2].plot(phi,out0['default'][f'NSR_{rp}_RP'],label=label1,color='r')
        ax[ind,ind2].plot(phi,out1['default'][f'NSR_{rp}_RP'],label=label2,color='g')
        ax[ind,ind2].plot(phi,out2['default'][f'NSR_{rp}_RP'],label=label3,color='b')
        # ax.plot(phi,mi_model[f'NSR_{rp}_RP'],label='Michelson',color='purple')
        ax[ind,ind2].set_ylim(1e-25,1e-19)
        ax[ind,ind2].loglog()
        ax[ind,ind2].xaxis.set_major_locator(mticker.LogLocator(numticks=999))
        ax[ind,ind2].xaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
        ax[ind,ind2].yaxis.set_major_locator(mticker.LogLocator(numticks=6))
        ax[ind,ind2].yaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
fig.supxlabel('Signal Frequency [Hz]',fontsize=10,font='Times New Roman')
# ax[1].set_ylabel('Strain Sensitivity [/$\sqrt{Hz}$]',fontsize=11)
ax[0,0].set_title(r'$\delta_{\text{off}}=0$'+'\na.) Lossless Mirrors',fontsize=10,font='Times New Roman',pad=10)
ax[0,1].set_title(r'$\delta_{\text{off}}=0.001°$'+'\nb.) Lossless Mirrors',fontsize=10,font='Times New Roman',pad=10)
ax[1,0].set_title('c.) T=15x10$^{-6}$, L=5x10$^{-6}$',fontsize=10,font='Times New Roman',pad=10)
ax[1,1].set_title('d.) T=15x10$^{-6}$, L=5x10$^{-6}$',fontsize=10,font='Times New Roman',pad=10)

# plt.subplots_adjust(hspace=1)
# ax[1].legend(ncols=1,bbox_to_anchor=(1.05,.5),loc='center left')
# fig.legend(ncols=2,loc='outside lower center')
fig.supylabel(r'Strain Sensitivity $\left[\text{h}/\sqrt{\text{Hz}}\right]$',fontsize=10,font="Times New Roman",x=0)
fig.legend(ncols=2,bbox_to_anchor=(.5, 0.01),loc='upper center',
          bbox_transform=fig.transFigure,prop=FontProperties(family='times',size=8))
plt.tight_layout(pad=.2)
plt.savefig('Figure7.png',bbox_inches='tight',dpi=400)
plt.show()
