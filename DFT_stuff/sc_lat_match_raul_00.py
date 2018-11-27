#!/usr/bin/env python

"""Create heterointerfaces between graphene and slab surface.

TEST 181121 - RF
I'm working on this in my personal branch - Raul Flores

Todo:
    # Impletment loop over a list of defined surface cuts

Author(s): Kevin Krempl, Raul Flores
"""

#| - Import Modules
import pickle
import os

from mpinterfaces.transformations import (
    Structure,
    get_aligned_lattices,
    generate_all_configs,
    )

from mpinterfaces.utils import slab_from_file
from mpinterfaces.interface import Interface

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from ase import io
#__|

#| - Script Inputs
# strain_sys = "support"  # 'support' or 'overlayer'
strain_sys = "overlayer"  # 'support' or 'overlayer'

bulk_filename = 'Cobulk.cif'
graphene_filename = 'init_graphene.cif'

surface_cuts = [[0, 0, 1],[1,0,0],[1,0,1],[1,1,0],[2,1,0]]

separation = 3
nlayers_2d = 1
nlayers_substrate = 1

# Lattice matching algorithm parameters
max_area = 1000
max_mismatch = 0.02
max_angle_diff = 0.5
r1r2_tol = 0.01
#__|

#| -

out_dir = "01_heterostructures"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

mat2d_slab = slab_from_file([0, 0, 1], graphene_filename)
substrate_bulk = Structure.from_file(bulk_filename)

analyzer=SpacegroupAnalyzer(substrate_bulk)
substrate_bulk_conv=analyzer.get_conventional_standard_structure()

writer=CifWriter(substrate_bulk_conv)
writer.write_file("conv_" + bulk_filename)


#| - Loop over all surface surface_cuts

for surface_cut in surface_cuts:

    substrate_slab = Interface(
        substrate_bulk_conv,
        hkl=surface_cut,
        min_thick=10,
        min_vac=25,
        primitive=False,
        from_ase=True,
        )

    if strain_sys == "support":
        lower_mat = mat2d_slab
        upper_mat = substrate_slab
    elif strain_sys == "overlayer":
        lower_mat = substrate_slab
        upper_mat = mat2d_slab

    # mat2d_slab_aligned, substrate_slab_aligned = get_aligned_lattices(
    lower_mat_aligned, upper_mat_aligned = get_aligned_lattices(
        lower_mat,
        upper_mat,
        max_area=max_area,
        max_mismatch=max_mismatch,
        max_angle_diff=max_angle_diff,
        r1r2_tol=r1r2_tol,
        )

    # merge substrate and mat2d in all possible ways
    hetero_interfaces = generate_all_configs(
        lower_mat_aligned,
        upper_mat_aligned,
        nlayers_2d,
        nlayers_substrate,
        separation,
        )

    out_dir = os.path.join(
        '01_heterostructures',
        str(surface_cut[0]) + str(surface_cut[1]) + str(surface_cut[2])
        )

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for i, iface in enumerate(hetero_interfaces):
        atoms = AseAtomsAdaptor.get_atoms(iface)
        io.write(
            os.path.join(
                out_dir,
                'structure_' + str(i + 1).zfill(2) + '.traj',
                ),
            atoms,
            )
