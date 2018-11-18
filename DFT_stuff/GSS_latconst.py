#/usr/bin/env python
#| - SLURM HEADER
#above line selects special python interpreter which knows all the paths
#SBATCH -p iric,owners
#################
#set a job name
#SBATCH --job-name=recon1
#################
#a file for job output, you can check job progress
#SBATCH --output=myjob.out
#################
# a file for errors from the job
#SBATCH --error=myjob.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH --time=00:30:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#memory per node; default is 4000 MB per CPU
#SBATCH --mem-per-cpu=4000
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=ALL
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=kkrempl@stanford.edu
#################
#task to run per node; each node has 16 cores
#SBATCH --ntasks-per-node=16
#################
#__|
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

        "outdir": "calcdir_opt",
        "kpts": ["8","8","8"],
        "parflags": None,
        "xc": "BEEF-vdW",
        "nbands": -50,
        "convergence": {
            "nmix": 20,
            "diag": "david",
            "energy": 2e-06,
            "mixing_mode": "local-TF",
            "maxsteps": 500,
            "mixing": 0.2
            },
        "dw": 8000.0,
        "outdir": "calcdir",
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

    atoms.set_calculator(**params)
    return atoms.get_potential_energy()


print(energy_fcc(3.524))
