#____________________________________________________________________________#
# Nice CV plots from uncorrected raw data. Employs capacitive corrections and#
# iR correction. add line                                                   #
#                                                                            #
# Author: Kevin Krempl                                            10/06/18   #
#____________________________________________________________________________#

# Import modules
import numpy as np
from scipy.interpolate import interp1d
import pickle
import matplotlib.pyplot as plt

class cyclic_voltamogram:

    def __init__(self, cvfilepath, **kwargs):
        self.CVraw = np.genfromtxt(cvfilepath)
        self.solres=kwargs.get('solres', 0)
        if kwargs.get('capcorr', 'None')!='None':
            self.capcorr=np.genfromtxt(kwargs.get('capcorr'))
        else:
            self.capcorr=np.zeros(self.CVraw.shape)
        self.shift=kwargs.get('shift', 0)
        self.surfarea=kwargs.get('surfarea', 1)
        self.name=kwargs.get('name', 'empty')

    @property
    def CVcorr(self):
        solres = self.solres
        shift = self.shift
        sa = self.surfarea
        capcorr = self.capcorr

        CVcorr1= np.column_stack((self.CVraw[:,0]+shift+abs(self.CVraw[:,1])*solres/1000, self.CVraw[:,1]/sa))
        CVcorr2 = np.column_stack((self.capcorr[:,0]+shift+abs(self.capcorr[:,1])*solres/1000, self.capcorr[:,1]/sa))

        func = interp1d(CVcorr2[:,0], CVcorr2[:,1], fill_value = 'extrapolate')
        CVcorr3 = np.column_stack((CVcorr1[:,0],CVcorr1[:,1]-func(CVcorr1[:,0])))

        return [CVcorr1, CVcorr2, CVcorr3]

    def mkplot(self, figname, **kwargs):

        CV=self.CVcorr[2]
        plt.style.use('mystandart')

        plt.figure(figsize = (6.5, 5))
        plt.plot(CV[:,0], CV[:,1], 'b', label=self.name)
        if kwargs.get('range', 0)!= 0:
            plt.axis(kwargs.get('range'))
        plt.xlabel(u'$\it E$ ${[V_{RHE}]}$')
        #plt.xlabel(u'$\it E_{iR \u2010 corr}$ ${[V_{RHE}]}$')
        plt.ylabel(r'$\it i_{geo}$ $[\frac{mA}{cm^2}]$')
        plt.legend(loc = 4)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(figname + '.png')



class CV_plot:

    def __init__(self, CVs, **kwargs):
        self.CVs = CVs
        self.electrolyte=kwargs.get('electrolyte')
        self.gas=kwargs.get('gas')
        self.scanrate=kwargs.get('scanrate')

    def mkplot(self, figname, **kwargs):

        plt.style.use('mystandart')
        plt.figure(figsize = (6.5, 5))
        for CV in self.CVs:
            currentCV = CV.CVcorr[2]
            plt.plot(currentCV[:,0], currentCV[:,1], label=CV.name)

        if kwargs.get('range', 0)!= 0:
            plt.axis(kwargs.get('range'))
        plt.xlabel(u'$\it E$ ${[V_{RHE}]}$')
        plt.ylabel(r'$\it i_{geo}$ $[\frac{mA}{cm^2}]$')
        strConditions = r'$10\,\frac{mV}{s}$ in $0.1 \, M$ $HClO_4$'
        plt.figtext(0.6,0.5, strConditions, ha = 'left', va = 'center', bbox=dict(facecolor='white', edgecolor='black', linestyle='--'))
        plt.legend(loc = 4)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(figname + '.png')
