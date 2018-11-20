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
"""findes matching unit cells for graphene and a previously optimized bulk crystal
"""

#| - Import Modules
from mpinterfaces.interface import Interface
from mpinterfaces.transformations import *
from mpinterfaces.utils import *

from pymatgen.io.ase import AseAtomsAdaptor

from ase.io import write


#__|


#| - Inputs
bulk_filename = 'Cobulk.POSCAR'
graphene_filename = 'graph.POSCAR'

surface_cut = [0,0,1]

separation = 3
nlayers_2d = 2
nlayers_substrate = 2
# Lattice matching algorithm parameters
max_area = 400
max_mismatch = 4
max_angle_diff = 1
r1r2_tol = 0.01

#__|

#| - Generate heterstructures

#impletment loop over a list of defined surface cuts

substrate_bulk = Structure.from_file(bulk_filename)
substrate_slab = Interface(substrate_bulk, hkl = surface_cut, min_thick = 10, min_vac = 25, primitive = False, from_ase = True)
mat2d_slab = slab_from_file([0,0,1], graphene_filename)

# get aligned lattices
substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(substrate_slab, mat2d_slab, max_area = max_area, max_mismatch = max_mismatch, max_angle_diff = max_angle_diff, r1r2_tol = r1r2_tol)

# merge substrate and mat2d in all possible
#ways
hetero_interfaces = generate_all_configs(mat2d_slab_aligned, substrate_slab_aligned, nlayers_2d, nlayers_substrate, separation)

hetero_interfaces_ase = AseAtomsAdaptor.get_atoms(hetero_interfaces)

write('heterostructure.traj', hetero_interfaces_ase)

#__|

#| - Main loop
