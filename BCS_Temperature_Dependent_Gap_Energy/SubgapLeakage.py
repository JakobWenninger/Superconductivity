import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad
from scipy.constants import k,hbar, elementary_charge
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.special import kn
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

rn=13 # Ohm
temperatures=[3,4,5,6,7]
delta =  np.multiply(10./2. * k * 3.5 , [.997,.985,.958,.907,.828])# values form implicit BCS theory
#delta =  np.multiply(10./2. * k * 3.5 , np.sqrt(1-np.divide(temperatures,10)))
#delta =  379e-6 *elementary_charge#Al

eV = lambda v: np.multiply(elementary_charge,v)
eV2kT = lambda v,t : np.divide(eV(v),2*k*t)
temdecay = lambda t,d : np.exp(-np.divide(d,np.multiply(t,k))) # constant against voltage
sqr = lambda v,d: np.sqrt(np.divide(2*d,np.add(eV(v),2*d)))
sinhK0 = lambda v,t: np.multiply(np.sinh(eV2kT(v,t)),kn(0,eV2kT(v,t)))#approaches expectation value
current = lambda v,t,d : np.multiply(np.multiply(np.multiply(np.multiply(
                                            #np.divide(2,np.multiply(1,rn)),
                                            np.divide(2,
                                                  np.multiply(elementary_charge,rn)),
                                            temdecay(t,d)),
                                            sqr(v,d)),
                                            np.add(eV(v),d)),
                                            sinhK0(v,t))
#unit is ok too
voltages=np.arange(0,3e-3,0.001e-3)
plt.figure()
for t in range(len(temperatures)): 
    plt.plot(voltages*1000,current(voltages,temperatures[t],delta[t])*1e6,label='$T = %r\,\mathrm{K}$'%temperatures[t])
#plt.plot(voltages,sqr(voltages),label='sqr')
#plt.axvline(np.divide(delta,elementary_charge))#ok
plt.xlabel('$V_0$ [mV]')
plt.ylabel('$I_\mathrm{SG}\,\,[\mathrm{\mu A}]$')
plt.grid(True)
plt.xlim(0,2.0)
plt.ylim(0,25)
plt.rcParams.update({'font.size': 12})
legend = plt.legend(loc='best', shadow=False,ncol=1)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()  # all the text.Text instance in the legend
llines = leg.get_lines()  # all the lines.Line2D instance in the legend
plt.setp(ltext, fontsize='small')
plt.setp(llines, linewidth=1.5)      # the legend linewidth
plt.tight_layout()
plt.savefig('SubgapLeakage.pdf', bbox_inches="tight")
plt.show()
