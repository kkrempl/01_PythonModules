#!/usr/bin/env python
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
#SBATCH --mail-type=FAIL
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

import scipy.optimize as opt

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
init_lat_const= [3.5, 3.5]
calc= espresso(
    output = {
        'avoidio': True,
        'removesave': True,
        'removewf': True,
        'wf_collect': False,
        },

    kpts = (6,6,6),
    parflags = None,
    xc = 'BEEF-vdW',
    nbands = -50,
    convergence = {
        'nmix': 20,
        'diag': 'david',
        'energy': 2e-06,
        'mixing_mode': 'local-TF',
        'maxsteps': 500,
        'mixing': 0.2
        },
    dw = 6000.0,
    outdir = 'calcdir',
    pw = 600,
    noncollinear = False,
    dipole = {
        'status': False
        },
    beefensemble = True,
    printensemble = True,
    sigma = 0.005,
    spinpol = True
    )

#__|

#| - Functions calculating potential energy for given lattice constants

def energy_1D(lat_const):

    atoms=bulk(name, crystalstructure, a=lat_const[0])
    atoms.set_calculator(calc)
    atoms.set_initial_magnetic_moments([6])
    sol = atoms.get_potential_energy()
    print(str(sol) + '; ' + str(lat_const))
    return sol

def energy_2D(lat_const):

    atoms=bulk(name, crystalstructure, a=lat_const[0], c=lat_const[1])
    atoms.set_calculator(calc)
    atoms.set_initial_magnetic_moments([6])
    sol = atoms.get_potential_energy()
    print(str(sol) + '; ' + str(lat_const[0]) + '; ' + str(lat_const[1]))
    return sol
#__|

#| - Minimize potential energy
if crystalstructure == 'hcp':
    res = opt.minimize(energy_2D,
                            init_lat_const,
                            method='BFGS',
                            tol=1e-6
                            )
    atoms_out = bulk(name, crystalstructure, a=res.x[0], c=res.x[1])

else:
    res = opt.minimize_scalar(energy_1D,
                            bracket=(init_lat_const[0]-0.2, init_lat_const[0]-0.2),
                            method='Brent',
                            tol=1e-6
                            )
    atoms_out = bulk(name, crystalstructure, a=res.x)

print(res.x)

write(name+'_'+crystalstructure+'_bulk_optimized.traj', atoms_out)

#__|
