from ase.spacegroup import crystal
from ase.spacegroup import Spacegroup as ASESpacegroup
from ase.atoms import Atoms
import numpy as np

def getSpeciesListCIF(inputFile):
    ''' inputFile is the name of an input file 
        use ASE at this point to quickly parse file 
        DOES NOTHING ABOUT FRACTIONAL OCCUPANCY '''

    from ase.io.cif import parse_cif
    aseParser = parse_cif(inputFile)

    for name, c in aseParser:
        scaled_positions = np.array([c['_atom_site_fract_x'],
                                     c['_atom_site_fract_y'],
                                     c['_atom_site_fract_z']]).T
        crys = crystal(c['_atom_site_type_symbol'],
                        basis=scaled_positions,
                        cellpar=[c['_cell_length_a'], c['_cell_length_b'], c['_cell_length_c'],
                                 c['_cell_angle_alpha'], c['_cell_angle_beta'], c['_cell_angle_gamma']],
                        spacegroup=c['_symmetry_int_tables_number'])
        atoms = Atoms(symbols = c['_atom_site_type_symbol'],
                      scaled_positions=scaled_positions,
                      cell = [c['_cell_length_a'], c['_cell_length_b'], c['_cell_length_c'],
                              c['_cell_angle_alpha'], c['_cell_angle_beta'], c['_cell_angle_gamma']],
                      info = {'spacegroup': ASESpacegroup(c['_symmetry_int_tables_number'])})

        print 'pbc=False', len(crys), len(atoms)
        yield atoms
#getSpeciesListCIF('testingScripts/CSPPCM_example.cif')
