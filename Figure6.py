from matplotlib import pyplot as plt, ticker as mticker
from numpy import pi, sin,cos
import numpy as np
from GMI_Finesse_Funcs import *
from scipy.constants import c,epsilon_0,hbar
finesse.configure(plotting=True)
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.style.use('bmh')
plt.rcParams.update({'font.family': 'times'})

power = 125
e0 = (2*power/(c*epsilon_0))**.5
wlen = 1.064e-6
f = np.logspace(0,6,1001)
k = 2*pi*f/c
f0 = 1/wlen
k0 = f0*2*pi
omega0 = k0*c
Larm = 4e3
k0deltaL = 1.282*pi/180
k0deltaL2 = 0.0006708*pi/180
noise = (power*hbar*omega0*np.sin(3.4*k0deltaL)**2)**0.5
print(noise)

# agw = e0*abs(c*f0*sin(k*Larm)*sin(2*k0deltaL)/f)
agw = power*abs(c*f0*sin(k*Larm)*sin(2*k0deltaL)/f)

mi_ = mi(power=power,lossy=0,nsr=1)
mi_model = mi_.run("xaxis(fsig.f, log, 1, 1M, 1000)")

gmi1 = gmi(power=power,phi2=89.985,phi1=0.015,lossy=0,pendula=0,nsr=1)
out1 = gmi1.run("xaxis(fsig.f, log, 1, 1M, 1000)")

r1 = np.sqrt(1)
r2 = np.sqrt(1)
r3 = np.sqrt(1)
ec1 = lambda phi1,phi2,phi3: 0.5 * (1 - r2*r3*np.exp(1j*(phi2+phi3)))/(1 - 0.5 * r3*np.exp(1j*phi3)*(r1*np.exp(1j*phi1) + r2*np.exp(1j*phi2)))
ec2 = lambda phi1,phi2,phi3: 0.5 * (1 - r1*r3*np.exp(1j*(phi1+phi3)))/(1 - 0.5 * r3*np.exp(1j*phi3)*(r1*np.exp(1j*phi1) + r2*np.exp(1j*phi2)))
eigamma = lambda phi1,phi2,phi3: (r1*np.exp(1j*phi1) + r2*np.exp(1j*phi2) - 2*r1*r2*r3*np.exp(1j*(phi1+phi2+phi3))) / (r3*np.exp(1j*phi3)*(r1*np.exp(1j*phi1) + r2*np.exp(1j*phi2)) - 2)
et = lambda phi1,phi2,phi3: 0.5 * (eigamma(phi1,phi2,phi3) + 1)
et_ = lambda phi1,phi2,phi3: cos(np.angle(eigamma(phi1,phi2,phi3))/2)**2
er = lambda phi1,phi2,phi3: 0.5 * (eigamma(phi1,phi2,phi3) - 1)
er_ = lambda phi1,phi2,phi3: sin(np.angle(eigamma(phi1,phi2,phi3))/2)**2

phi2 = .03*pi/180
phi3 = 0
phi1 = 2*pi-phi2

ec = abs((ec1(phi1,phi2,phi3)-ec2(phi1,phi2,phi3))/2)
agw_gmi = ec**2*power*abs(c*f0*sin(k*Larm)*sin(2*k0deltaL2)/f)

mi_nsr = abs((2*pi*f/np.sin(k*Larm))*(hbar/(2*power*omega0))**0.5)
gmi_nsr = abs((2*pi*f/np.sin(k*Larm))*(hbar/(2*ec**2*power*omega0))**0.5)

fig, ax = plt.subplot_mosaic([['left'], ['right']],figsize=(3,3))

ax['left'].plot(f,agw,color='b',label='MI (Analytical)')
ax['left'].plot(f,abs(mi_model['signal']),linestyle='dashed',color='grey',label='MI (Finesse)')
ax['left'].plot(f,agw_gmi,color='r',label='GMI (Analytical)')
ax['left'].plot(f,abs(out1['signal']),color='k',label='GMI (Finesse)',linestyle='dashed')
ax['left'].loglog()
# ax['left'].set_xlabel('Signal Frequency [Hz]',fontsize=12,font='Times New Roman')
ax['left'].set_ylabel(r'Magnitude [W/h]',fontsize=10,font='Times New Roman')
# plt.plot(f,noise/agw)
ax['right'].plot(f,mi_nsr,color='b')#,label='Michelson (Analytical)')
ax['right'].plot(f,mi_model['NSR_without_RP'],linestyle='dashed',color='grey')#,label='Michelson (Finesse)')
# plt.plot(f,noise/agw_gmi)
ax['right'].plot(f,gmi_nsr,color='r')#,label='Grover-Michelson (Analytical)')
ax['right'].plot(f,out1['NSR_without_RP'],color='k',linestyle='dashed')#,label='Grover-Michelson (Finesse)')
ax['right'].set_xlabel('Signal Frequency [Hz]',fontsize=10,font='Times New Roman')
ax['right'].set_ylabel(r'NSR $\left[\text{h}/\sqrt{\text{Hz}}\right]$',fontsize=10,font='Times New Roman')
ax['right'].loglog()

for ind in ['left','right']:
    ax[ind].xaxis.set_major_locator(mticker.LogLocator(numticks=999))
    ax[ind].xaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
    ax[ind].yaxis.set_major_locator(mticker.LogLocator(numticks=6))
    ax[ind].yaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
# fig.legend(loc='upper center',ncols=4,borderaxespad=0)
plt.tight_layout()
fig.legend(ncols=2,bbox_to_anchor=(.5, 0.05),loc='upper center',
          bbox_transform=fig.transFigure)
# plt.subplots_adjust(wspace=0.3)
plt.savefig('Figure6.png',bbox_inches='tight',dpi=400)
plt.show()