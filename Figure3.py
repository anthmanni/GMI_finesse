import matplotlib.style
import numpy as np
from matplotlib import pyplot as plt, ticker as mticker
from numpy import pi
import matplotlib as mp
mp.style.use('bmh')
plt.rcParams.update({'font.size': 12, 'font.family': 'times'})

r1 = np.sqrt(1)
r2 = np.sqrt(1)
r3 = np.sqrt(1)
ec1 = lambda phi1,phi2,phi3: 0.5 * (1 - r2*r3*np.exp(1j*(phi2+phi3)))/(1 - 0.5 * r3*np.exp(1j*phi3)*(r1*np.exp(1j*phi1) + r2*np.exp(1j*phi2)))
ec2 = lambda phi1,phi2,phi3: 0.5 * (1 - r1*r3*np.exp(1j*(phi1+phi3)))/(1 - 0.5 * r3*np.exp(1j*phi3)*(r1*np.exp(1j*phi1) + r2*np.exp(1j*phi2)))
pc = lambda phi1,phi2: (r1**2+r2**2-2*r1*r2*np.cos(phi1-phi2))/(16*(1-np.cos(phi1)-np.cos(phi2))+4*(r1**2+r2**2)+8*r1*r2*np.cos(phi1-phi2))
eigamma = lambda phi1,phi2,phi3: (r1*np.exp(1j*phi1) + r2*np.exp(1j*phi2) - 2*r1*r2*r3*np.exp(1j*(phi1+phi2+phi3))) / (r3*np.exp(1j*phi3)*(r1*np.exp(1j*phi1) + r2*np.exp(1j*phi2)) - 2)
et = lambda phi1,phi2,phi3: 0.5 * (eigamma(phi1,phi2,phi3) + 1)
et_ = lambda phi1,phi2,phi3: np.cos(np.angle(eigamma(phi1,phi2,phi3))/2)**2
er = lambda phi1,phi2,phi3: 0.5 * (eigamma(phi1,phi2,phi3) - 1)
er_ = lambda phi1,phi2,phi3: np.sin(np.angle(eigamma(phi1,phi2,phi3))/2)**2

phi2 = [pi/16,pi/8,pi/4,pi/2,pi] # fano at phi_2,phi_3 --> pi
phi3 = 0
phi_1 = np.linspace(1e-5,2*pi,10000)
clrs = ['b','g','orange','r','k']
gamma = lambda phi1,phi2,phi3: np.angle(eigamma(phi1,phi2,phi3))

fig,ax = plt.subplots(2,1,figsize=(3,5))
for ind,j in enumerate(phi2):
    _eig = [eigamma(i,j,phi3) for i in phi_1]
    _gam = np.angle(_eig)
    r = np.sin(_gam/2)**2
    t =  np.cos(_gam/2)**2
    ec = [abs((ec1(i, j, phi3) - ec2(i, j, phi3)) / 2)**2 for i in phi_1]
    label = '0' if j<.001 else f'π/{np.round(pi/j,0):.0f}'
    label = 'π' if j==pi else label
    ax[0].plot(phi_1/pi,_gam/pi,label=f'$\phi_N$={label}',color=clrs[ind])
    # ax[1].plot(phi_1 / pi, r,label=f'$\phi_N$={label}')
    ax[1].plot(phi_1 / pi, t,label=f'$\phi_N$={label}',color=clrs[ind])
    # ax[2].plot(phi_1 / pi, ec, label=f'$\phi_N$={label}')
ax[0].set_xlabel('$\phi_E$ ($\pi$ units)')
ax[1].set_xlabel('$\phi_E$ ($\pi$ units)')
ax[0].set_ylabel('Phase $\gamma(\phi_N,\phi_E)$')
ax[1].set_ylabel('Transmittance')
# ax[2].set_xlabel('$\phi_E$ ($\pi$ units)')
# ax[2].set_xlabel('$\phi_E$ ($\pi$ units)')
# ax[2].set_xlabel('$\phi_E$ ($\pi$ units)')
# ax[1].set_ylabel('Reflectance')
# ax[2].set_ylabel(r'$|E_c|$')
# ax[2].set_xlim(0,2)
# ax[2].set_ylim(0,30)
[ax[i].set_xlim(0,2) for i in [0,1]]
[ax[i].set_ylim(0,1) for i in [1]]
ax[0].set_ylim(-1,1)
ax[0].set_yticks([-1,-0.5,0,0.5,1],['-π','-π/2',0,'π/2','π'])
ax[0].legend(bbox_to_anchor=(1.04, .5), loc="center left")
ax[1].legend(bbox_to_anchor=(1.04, .5), loc="center left")
for ind in [0,1]:
    ax[ind].xaxis.set_major_locator(mticker.LinearLocator(5))
    ax[ind].xaxis.set_minor_locator(mticker.LinearLocator(21))
    ax[ind].yaxis.set_major_locator(mticker.LinearLocator(5))
    ax[ind].yaxis.set_minor_locator(mticker.LinearLocator(21))
plt.subplots_adjust(hspace=0.5)

plt.savefig('Figure3.png',bbox_inches='tight',dpi=400)
plt.show()