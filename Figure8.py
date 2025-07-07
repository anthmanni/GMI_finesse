from matplotlib import pyplot as plt, ticker as mticker
import numpy as np
from GMI_Finesse_Funcs import *
finesse.configure(plotting=True)
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.style.use('bmh')
plt.rcParams.update({'font.family': 'times','text.latex.preamble':r"\usepackage{siunitx}"})

gmipower = 125
ligopower =  gmipower
phi_2 = 89.9
phi_1 = 90-phi_2
lossy = 1
nsr = 1
phi0 = [0,90]
phi3 = 0
rp = 'with'
prmt = 0.001

ligo_ = ligo(power=ligopower,lossy=lossy)
ligo_model = ligo_.run("xaxis(fsig.f, log, 1, 1M, 1000)")
mi_ = mi(power=gmipower,lossy=lossy,nsr=nsr)
mi_model = mi_.run("xaxis(fsig.f, log, 1, 1M, 1000)")

gmi0 = [gmi(power=gmipower,phi2=phi_2,phi1=phi_1,phi0=i,phi3=phi3,lossy=lossy) for i in phi0]
gmi1 = [gmi(power=gmipower,phi2=phi_2,phi1=phi_1,phi0=i,prmT=prmt,phi3=phi3,lossy=lossy,pr=1) for i in phi0]
gmi2 = [gmi(power=gmipower,phi2=phi_2,phi1=phi_1,phi0=i,phi3=phi3,lossy=lossy,sr=1) for i in phi0]
gmi3 = [gmi(power=gmipower,phi2=phi_2,phi1=phi_1,phi0=i,phi3=phi3,prmT=.03,lossy=lossy,sr=1,pr=1) for i in phi0]

out0=[i.run("parallel(xaxis(fsig.f, log, 1, 1M, 1000,name='default'),xaxis(phi1,lin,0,1,500,name='trans'))") for i in gmi0]
out1=[i.run("parallel(xaxis(fsig.f, log, 1, 1M, 1000,name='default'),xaxis(phi1,lin,0,1,500,name='trans'))") for i in gmi1]
out2=[i.run("parallel(xaxis(fsig.f, log, 1, 1M, 1000,name='default'),xaxis(phi1,lin,0,1,500,name='trans'))") for i in gmi2]
out3=[i.run("parallel(xaxis(fsig.f, log, 1, 1M, 1000,name='default'),xaxis(phi1,lin,0,1,500,name='trans'))") for i in gmi3]

phi = np.logspace(0,6,1001)

fig, ax = plt.subplots(3,1,figsize=(3,6),constrained_layout=True)
namedict = {0:'PR-GMI',1:'SR-GMI',2:'DR-GMI'}

for ind,i in enumerate([out1,out2,out3]):
    ax[ind].plot(phi,ligo_model['NSR_with_RP'],label='aLIGO',color='k')
    ax[ind].plot(phi,out0[0]['default']['NSR_with_RP'],label='GMI',color='r')
    # ax[ind].plot(phi, out0[1]['default']['NSR_with_RP'], label='Standard GMI, $\phi_1=\pi$', color='white',linestyle='dotted')
    ax[ind].plot(phi,i[0]['default']['NSR_with_RP'],label=namedict[ind]+'\n$\phi_1=0$',color='b')
    ax[ind].plot(phi, i[1]['default']['NSR_with_RP'], label=namedict[ind]+'\n$\phi_1=\pi$', color='g')
    ax[ind].loglog()
    ax[ind].legend(ncols=2,facecolor='None',frameon=False,fontsize=8)
    ax[ind].set_ylim(1e-24,1e-19)
    ax[ind].set_xlim(1, 5e4)
    ax[ind].xaxis.set_major_locator(mticker.LogLocator(numticks=999))
    ax[ind].xaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
ax[2].set_xlabel('Signal Frequency [Hz]',fontsize=10)
fig.supylabel(r"Strain Sensitivity $\left[\text{h}/\sqrt{\text{Hz}}\right]$",fontsize=10)
plt.savefig('Figure8.png',bbox_inches='tight',dpi=500)
# plt.tight_layout()
plt.show()



