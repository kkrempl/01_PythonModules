"""Finding dft-calculated lattice constants with golden section search
"""

#| - Import Modules
from ase.build import bulk
from ase.io import write
from espresso import espresso

#__|

#| - Inputs

"""
name: str
        Chemical symbol or symbols as in 'MgO' or 'NaCl'.
crystalstructure: str
        Must be one of sc, fcc, bcc, hcp, diamond, zincblende,
        rocksalt, cesiumchloride, fluorite or wurtzite.
intervall: [x,y]
        Lower (x) and upper (y) bound of the lattice constant.

"""

name = 'Ni'
crystalstructure = 'fcc'
intervall = [0.5,0.5]


#__|

#| - Function calculating potential energy for given lattice constants

def energy_fcc(a,c=None):

    atoms=bulk(name, crystalstructure, a=a, c=c)

    params= {
        "output": {
            "avoidio": True,
            "removesave": True,
            "removewf": True,
            "wf_collect": False,
            },

        "outdir": 'calcdir_opt',
        "kpts": (8,8,8),
        "parflags": None,
        "xc": 'BEEF-vdW',
        "nbands": -50,
        "convergence": {
            "nmix": 20,
            "diag": 'david',
            "energy": 2e-06,
            "mixing_mode": 'local-TF',
            "maxsteps": 500,
            "mixing": 0.2
            },
        "dw": 8000.0,
        "outdir": 'calcdir',
        "pw": 800,
        "noncollinear": False,
        "dipole": {
            "status": True
            },
        "beefensemble": True,
        "printensemble": True,
        "output": {
            "removesave": True
            },
        "sigma": 0.005,
        "spinpol": True
        }

    atoms.set_calculator(params)
    clean_up_dft()
    return atoms.get_potential_energy()


print(energy_fcc(3.524))
