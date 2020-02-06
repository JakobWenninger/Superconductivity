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

#Delta is 1
nSnN = lambda eOverDelta : np.divide(eOverDelta,np.sqrt(np.subtract(np.multiply(eOverDelta,eOverDelta),1)))


energies = np.arange(1.001,5,0.001)# avoid 1 since it causes division by 0
#plot result
plt.figure()
plt.plot(energies,nSnN(energies),label=r'''Superconductor $N_\mathrm{S} = N_\mathrm{N}(0\,\mathrm{eV})\cdot\frac{|E|}{\sqrt{E^2-\Delta^2}}$''') 
#plt.plot(tTC,deltaDelta0apprx,label = r'$T \approx T_\mathrm{C}$ Approximation')
plt.hlines(y=1,xmin=0,xmax=5,linestyle='--', label = r'Normal Metal $N_\mathrm{N}(\xi = 0\,\mathrm{eV})$')
rect = plt.Rectangle((0,0),1,5, color='g',hatch='//',fill=False, label='Band Gap Region')
plt.gca().add_patch(rect)
plt.xlabel(r'${E}/{\Delta}$')
plt.ylabel(r'$N_\mathrm{S}(E)/N_\mathrm{N}(\xi =0\,\mathrm{eV})$')
plt.grid(True)
plt.xlim(0,4)
plt.ylim(0,5)
plt.rcParams.update({'font.size': 12})
handles,labels = plt.gca().get_legend_handles_labels()
order = [0,2,1]
legend = plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='best', shadow=False,ncol=1)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()  # all the text.Text instance in the legend
llines = leg.get_lines()  # all the lines.Line2D instance in the legend
plt.setp(ltext, fontsize='small')
plt.setp(llines, linewidth=1.5)      # the legend linewidth
plt.tight_layout()
plt.savefig('Superconducting_States_over_Normal_States.pdf', bbox_inches="tight")
plt.show()
