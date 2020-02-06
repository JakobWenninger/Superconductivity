import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad
from scipy.constants import k,hbar, elementary_charge
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid",
              {'axes.edgecolor': '.2',
               'axes.facecolor': 'white',
               'axes.grid': True,
               'axes.linewidth': 0.5,
               'figure.facecolor': 'white',
               'grid.color': '.8',
               'grid.linestyle': u'-',
               'legend.frameon': True,
               'xtick.color': '.15',
               'xtick.direction': u'in',
               'xtick.major.size': 3.0,
               'xtick.minor.size': 1.0,
               'ytick.color': '.15',
               'ytick.direction': u'in',
               'ytick.major.size': 3.0,
               'ytick.minor.size': 1.0,
               })
sns.set_context("poster")

tDebye =  276 # K
omegaDebye = k/hbar*tDebye
tC = 9.2 # K
delta0 = 1.76 * k * tC

#approximated equation Tinkham eq. 3.54
apprx = lambda tau : 1.74*np.sqrt(1 - tau)

#the function to be integrated / BCS temperature dependent energy gap; Tinkham eq. 3.53 
#https://physics.stackexchange.com/questions/54200/superconducting-gap-temperature-dependence-how-to-calculate-this-integral
#Equation to integrate
funcToIntegrate = lambda z , sig , tau: np.divide(np.tanh(np.multiply(.882 , np.multiply(np.divide(sig , tau),np.sqrt(1 + np.multiply(z , z))))),np.sqrt(1 + np.multiply(z , z)))
#Equation to Solve
funcToSolve = lambda sig , *tau: + quad(funcToIntegrate,0,np.divide(hbar*omegaDebye,np.multiply(delta0,sig)),args=(sig,tau))[0] -np.arcsinh(hbar*omegaDebye/delta0) 

tTC =np.arange(0,1.001,.001) # T/TC values
deltaDelta0 = []
deltaDelta0apprx = []
for i in tTC: # probe for T/TC values
    deltaDelta0.append(fsolve(funcToSolve,.7, args=(i))) #.5 is a first guess

deltaDelta0apprx =(apprx(tTC))
deltaDelta0apprx = np.hstack([deltaDelta0apprx,0])
#add a little in normal state
tTC=np.hstack([tTC,1.2])
deltaDelta0.append(0)
deltaDelta0 = np.array(deltaDelta0)
deltaDelta0apprx = np.array(deltaDelta0apprx)
#plot result
plt.plot(tTC,deltaDelta0,label = 'Self Consistent Equation') 
plt.plot(tTC,deltaDelta0apprx,label = r'$T \approx T_\mathrm{C}$ Approximation')
plt.hlines(y=0.95,xmin=0,xmax=1.2, linestyle='--',label = r'${\Delta (T)}/{\Delta_0(T=0\,\mathrm{K})}=0.95$')
plt.xlabel(r'${T}/{T_\mathrm{C}}$')
plt.ylabel(r'${\Delta}/{\Delta_0(T=0\,\mathrm{K})}$')
plt.grid(True)
plt.xlim(0,1.2)
plt.ylim(0,1.2)
plt.rcParams.update({'font.size': 12})
legend = plt.legend(loc='best', shadow=False,ncol=1)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()  # all the text.Text instance in the legend
llines = leg.get_lines()  # all the lines.Line2D instance in the legend
plt.setp(ltext, fontsize='small')
plt.setp(llines, linewidth=1.5)      # the legend linewidth
plt.tight_layout()
plt.savefig('BCS_Temperature_Dependent_Gap_Energy.pdf', bbox_inches="tight")
plt.show()
