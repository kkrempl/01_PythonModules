#____________________________________________________________________________#
# Nice CV plots from uncorrected raw data. Employs capacitive corrections and#
# iR correction. Test test   changes                                                      #
#                                                                            #
# Author: Kevin Krempl                                            10/06/18   #
#____________________________________________________________________________#

# Import modules
import matplotlib
import numpy as np
from scipy.interpolate import interp1d
import pickle

class cyclic_voltamogram:

    def __init__(self, cvfilepath, **kwargs):
        self.CVraw = np.genfromtxt(cvfilepath)
        self.solres=kwargs.get('solres')
        if kwargs.get('capcorr', 'None')!='None':
            self.capcorr=np.genfromtxt(kwargs.get('capcorr'))
        else:
            self.capcorr=np.zeros(self.CVraw.shape)
        self.shift=kwargs.get('shift')
        self.surfarea=kwargs.get('surfarea')
        self.options=kwargs.get('options',[True, True, True, True])

    @property
    def CVcorr(self):
        solres = self.solres
        shift = self.shift
        sa = self.surfarea
        capcorr = self.capcorr

        if self.options[0]==False:
            solres = 0
        if self.options[1]==False:
            shift = 0
        if self.options[2]==False:
            sa = 1
        if self.options[3]==False:
            capcorr = np.zeros(self.CVraw.shape)

        CVcorr1= np.column_stack((self.CVraw[:,0]+shift+abs(self.CVraw[:,1])*solres/1000, self.CVraw[:,1]/sa))
        CVcorr2 = np.column_stack((self.capcorr[:,0]+shift+abs(self.capcorr[:,1])*solres/1000, self.capcorr[:,1]/sa))

        func = interp1d(CVcorr2[:,0], CVcorr2[:,1], fill_value = 'extrapolate')
        CVcorr3 = np.column_stack((CVcorr1[:,0],CVcorr1[:,1]-func(CVcorr1[:,1])))

        return [CVcorr1, CVcorr2, CVcorr3]
