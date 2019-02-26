#____________________________________________________________________________#
# Nice CV plots from uncorrected raw data. Employs capacitive corrections and#
# iR correction. add line                                                   #
#                                                                            #
# Author: Kevin Krempl                                            10/06/18   #
#____________________________________________________________________________#

# Import modules
import numpy as np
from scipy.interpolate import interp1d
import brewer2mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

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
        if self.capcorr[0,0]==0:
            CVcorr3=CVcorr1
        else:
            func = interp1d(CVcorr2[:,0], CVcorr2[:,1], fill_value = 'extrapolate')
            CVcorr3 = np.column_stack((CVcorr1[:,0],CVcorr1[:,1]-func(CVcorr1[:,0])))

        return [CVcorr1, CVcorr2, CVcorr3]

    def mkplot(self, figname, **kwargs):

        CV=self.CVcorr[kwargs.get('scan')]
        plt.style.use('mystandart')
        ax=plt.gca()

        plt.figure(figsize = (6.5, 5))
        plt.plot(CV[:,0], CV[:,1], 'g', label=self.name)
        if kwargs.get('range', 0)!= 0:
            plt.axis(kwargs.get('range'))
        plt.xlabel(u'$\it E_{iR\u2010free}$ ${[V_{RHE}]}$')
        #plt.xlabel(u'$\it E_{iR \u2010 corr}$ ${[V_{RHE}]}$')
        plt.ylabel(r'$\it i_{geo}$ $[\frac{mA}{cm^2}]$')
        plt.legend(loc = 4)
        strConditions = kwargs.get('scanrate') + r'$\,\frac{mV}{s}$ in ' + kwargs.get('electrolyte')
        plt.text(0.6, 0.4, strConditions, transform=ax.transAxes, ha = 'left', va = 'top', bbox=dict(facecolor='white', edgecolor='black', linestyle='--'))
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(figname + '.png')



class CV_plot:

    def __init__(self, CVs, **kwargs):
        self.CVs = CVs
        self.electrolyte=kwargs.get('electrolyte', '')
        self.gas=kwargs.get('gas', '')
        self.scanrate=kwargs.get('scanrate', '')

    def mkplot(self, figname, **kwargs):

        plt.style.use('mystandart')
        plt.figure(figsize = (6.5, 5))
        for CV in self.CVs:
            currentCV = CV.CVcorr[2]
            plt.plot(currentCV[:,0], currentCV[:,1], label=CV.name)

        ax=plt.gca()
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        if kwargs.get('range', 0)!= 0:
            plt.axis(kwargs.get('range'))
        plt.xlabel(u'$\it E_{iR\u2010corr}$ ${[V_{RHE}]}$')
        plt.ylabel(r'$\it i_{geo}$ $[\frac{mA}{cm^2}]$')
        strConditions = self.scanrate + r'$\,\frac{mV}{s}$ in ' + self.electrolyte
        plt.text(0.05, 0.05, strConditions, transform=ax.transAxes, ha = 'left', va = 'bottom', bbox=dict(facecolor='white', edgecolor='black', linestyle='--'))
        plt.legend(loc = 4)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(figname + '.png')


class Raman_spectra:

    def __init__(self, filepath, **kwargs):
        self.xydata = np.genfromtxt(filepath)
        self.name = kwargs.get('name', 'undefined')
        self.range = kwargs.get('range', [500, 3000])
        self.D = kwargs.get('D', False)

    def mkplot(self, figname, **kwargs):
        xydata=self.xydata
        plt.style.use('mystandart')

        plt.figure(figsize = (6.5, 5))
        plt.plot(xydata[:,0], xydata[:,1], 'red', label=self.name)
        plt.xlabel(u'Raman shift $[cm^{-1}]$')
        #plt.xlabel(u'$\it E_{iR \u2010 corr}$ ${[V_{RHE}]}$')
        plt.ylabel(r'Raman Intensity [a.u.]')
        ax=plt.gca()
        ax.set_yticklabels('')
        ax.set_yticks([])
        ax.set_xlim(self.range)
        plt.text(0.42, 0.9, 'G', ha = 'right', va = 'top', transform=ax.transAxes)
        plt.text(0.86, 0.8, '2D', ha = 'right', va = 'top', transform=ax.transAxes)
        if self.D==True:
            plt.text(0.3, 0.8, 'D', ha = 'right', va = 'top', transform=ax.transAxes)
        plt.legend(loc = 1)
        #plt.grid(True)
        plt.tight_layout()
        plt.savefig(figname + '.png')

class Raman_plot:

    def __init__(self, spectras, **kwargs):
        self.spectras = spectras
        self.range = kwargs.get('range', [500, 3000])
        self.D = kwargs.get('D', False)

    def mkplot(self, figname, **kwargs):

        plt.style.use('mystandart')
        plt.figure(figsize = (6.5, 5))

        x=0
        for key, value in self.spectras.iteritems():
            xydata=value.xydata
            plt.plot(xydata[:,0], normalize_vec(xydata[:,1])+x, label=value.name)
            x=x+1

        plt.xlabel(u'Raman shift $[cm^{-1}]$')
        #plt.xlabel(u'$\it E_{iR \u2010 corr}$ ${[V_{RHE}]}$')
        plt.ylabel(r'Raman Intensity [a.u.]')
        ax=plt.gca()
        ax.set_yticklabels('')
        ax.set_yticks([])
        ax.set_xlim(self.range)
        plt.text(0.42, 0.9, 'G', ha = 'right', va = 'top', transform=ax.transAxes)
        plt.text(0.86, 0.8, '2D', ha = 'right', va = 'top', transform=ax.transAxes)
        if self.D==True:
            plt.text(0.3, 0.8, 'D', ha = 'right', va = 'top', transform=ax.transAxes)
        plt.legend(loc = 2)
        #plt.grid(True)
        plt.tight_layout()
        plt.savefig(figname + '.png')






class XRD_diffrac:

    def __init__(self, filepath, **kwargs):
        self.xydata = np.genfromtxt(filepath)
        self.name = kwargs.get('name', 'undefined')

    def mkplot(self, figname, **kwargs):
        xydata=self.xydata
        plt.style.use('mystandart')

        plt.figure(figsize = (6.5, 5))
        plt.plot(xydata[:,0], xydata[:,1], 'b', label=self.name)
        #ax.set_yticklabels('')
        #ax.set_yticks([])
        if kwargs.get('range', 0)!= 0:
            plt.axis(kwargs.get('range'))
        plt.xlabel(r'$\it 2\Theta$ Mo K$\alpha$ [$\degree$]')
        #plt.xlabel(u'$\it E_{iR \u2010 corr}$ ${[V_{RHE}]}$')
        plt.ylabel(r'$\it Intensity$ [a.u.]')
        plt.legend(loc = 1)
        #plt.grid(True)
        plt.tight_layout()
        plt.savefig(figname + '.png')

class XRD_plot:

    def __init__(self, diffracs, **kwargs):
        self.diffracs = diffracs
        self.range = kwargs.get('range', [15, 45])
        self.zoom = kwargs.get('zoom', 1)


    def mkplot(self, figname, **kwargs):

        n_diffracs = len(self.diffracs)
        plt.ion
        plt.style.use('mystandart')

        fig, axarr=plt.subplots(n_diffracs, 1, sharex=True, figsize = (6.5, 1.7*n_diffracs))

        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        plt.grid(False)
        plt.xlabel(r'$\it 2\Theta$ Mo K$\alpha$ [$\degree$]')
        plt.ylabel(r'$\it Intensity$ [a.u.]', labelpad=-20)
        colors = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666']
        for i in range(0,n_diffracs):
            current_diffrac = self.diffracs[i]
            current_ax = axarr[i]
            xydata = current_diffrac.xydata
            current_ax.plot(xydata[:,0], xydata[:,1], label=current_diffrac.name, color=colors[i])
            ymin, ymax = current_ax.get_ylim()
            current_ax.set_ylim([ymin*self.zoom, ymax*self.zoom])
            if i == 0:
                current_ax.spines['bottom'].set_visible(False)
                current_ax.tick_params(bottom=False, top=True)


            elif i == n_diffracs-1:
                current_ax.spines['top'].set_visible(False)
                current_ax.tick_params(bottom=True, top=False)

            else:
                current_ax.spines['bottom'].set_visible(False)
                current_ax.spines['top'].set_visible(False)
                current_ax.tick_params(bottom=False, top=False)

            current_ax.set_yticklabels('')
            current_ax.set_yticks([])
            current_ax.set_xlim(self.range)
            current_ax.vlines([19.657, 22.752, 32.469, 38.275, 40.044],0,1, transform=current_ax.get_xaxis_transform(), colors='red', linewidths=1, linestyles='dashed')
            current_ax.vlines([20.178, 28.687, 35.325, 41.018, 46.122],0,1, transform=current_ax.get_xaxis_transform(), colors='red', linewidths=1, linestyles='dashed')
            current_ax.vlines([13.732, 16.135, 19.468, 27.761, 25.48],0,1, transform=current_ax.get_xaxis_transform(), colors='blue', linewidths=1, linestyles='dashed')
            current_ax.legend(loc=1)

        txt=u'dotted red: $\gamma\u2010$Fe, dotted blue: $FeO_x$'
        fig.text(0.5, 0.01, txt, wrap=True, horizontalalignment='center')
        fig.tight_layout()
        fig.subplots_adjust(hspace=0)
        fig.savefig(figname + '.png', bbox_inches='tight')


def normalize_vec(vector):
    x = (vector-vector.min())/(vector.max()-vector.min())
    return x
