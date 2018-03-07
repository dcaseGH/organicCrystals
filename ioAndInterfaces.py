def ccdcCrystalToASE(ccdcCrystal):
    ''' from CCDC crystal class, make ASE Atoms class '''
    from ase import Atoms
    from ase.spacegroup import Spacegroup
    aseCrystal = Atoms([x.atomic_symbol     for x in ccdcCrystal.molecule.atoms],
                       [list(x.coordinates) for x in ccdcCrystal.molecule.atoms],
                       info = {'spacegroup': Spacegroup(*ccdcCrystal.spacegroup_number_and_setting)}
                       )
    aseCrystal.set_cell(list(ccdcCrystal.cell_lengths) + list(ccdcCrystal.cell_angles))
    return aseCrystal

def ccdcCrystalToRes(ccdcCrystal, fileName):
    ''' Untested '''
    from ccdc.io import _ResWriter
    rw = _ResWriter(fileName)
    rw.write(ccdcCrystal)
    
def ccdcCrystalsToCIF(listCrystals, fileName):
    ''' pass crystals and a name for the output file to this
        change name of each crystal with crystal.identifier 
        all cifs go into one fileName[example.cif] '''

    from ccdc.io import CrystalWriter
    
    #catch the error of passing one crystal not within list
    if type(listCrystals) != list:
        listCrystals = list(listCrystals)

    with CrystalWriter(fileName) as outWriter:
        for c in listCrystals:
            outWriter.write(c)

def resToASEAtoms(inFilename, spacegroup):
    ''' Must manually give a spacegroup (ie number or possible setting too??) '''

    from ase import Atoms
    from ase.spacegroup import Spacegroup as ASESpacegroup
    from multipoleFile import is_number

    def acceptableAtomicLine(line):
        ''' Each atom begins with a line like this: 
        the element and a coordinate '''

        if any([x in line for x in ['CELL', 'TITL', 'ZERR', 'LATT', 'SFAC']]):
            return False
        parts = line.split()

        # assume [number] [neighcrys string] [position] [stuff]
        if len(parts) > 4 and parts[0][0].isalpha() and all(map(is_number, parts[2:5])):
            return True
        else: 
            return False

    def earlyAlpha(l):
        ''' Take all alpha before something else (could use re) '''
        outString = ''
        for x in l:
            if x.isalpha():
                outString += x
            else:
                return outString

    atomSymbols, atomScaledPositions = [], []

    with open(inFilename, 'r') as inFile:
        for l in inFile.readlines():
            if 'CELL' in l:
                cell = map(float, l.split()[-6:])
            if acceptableAtomicLine(l):
                atomSymbols.append(earlyAlpha(l.split()[0]))
                atomScaledPositions.append(map(float, l.split()[2:5]))

    return Atoms(symbols = atomSymbols,
                 scaled_positions = atomScaledPositions,
                 cell = cell,
                 info = {'spacegroup': ASESpacegroup(spacegroup)})


def aseWriteCifNoWrap(fileobj, images, format='default'):
    """Write *images* to CIF file.
       DOES NOT WRAP ATOMS - THIS MAINTAINS MOLECULES """
    if isinstance(fileobj, basestring):
#        fileobj = paropen(fileobj, 'w')
        fileobj = open(fileobj, 'w')

    if hasattr(images, 'get_positions'):
        images = [images]

    for i, atoms in enumerate(images):
        fileobj.write('data_image%d\n' % i)

        a, b, c, alpha, beta, gamma = atoms.get_cell_lengths_and_angles()

        if format == 'mp':

            comp_name = atoms.get_chemical_formula(mode='reduce')
            sf = split_chem_form(comp_name)
            formula_sum = ''
            ii = 0
            while ii < len(sf):
                formula_sum = formula_sum + ' ' + sf[ii] + sf[ii + 1]
                ii = ii + 2

            formula_sum = str(formula_sum)
            fileobj.write('_chemical_formula_structural       %s\n' %
                          atoms.get_chemical_formula(mode='reduce'))
            fileobj.write('_chemical_formula_sum      "%s"\n' % formula_sum)

        fileobj.write('_cell_length_a       %g\n' % a)
        fileobj.write('_cell_length_b       %g\n' % b)
        fileobj.write('_cell_length_c       %g\n' % c)
        fileobj.write('_cell_angle_alpha    %g\n' % alpha)
        fileobj.write('_cell_angle_beta     %g\n' % beta)
        fileobj.write('_cell_angle_gamma    %g\n' % gamma)
        fileobj.write('\n')
        if atoms.pbc.all():
            fileobj.write('_symmetry_space_group_name_H-M    %s\n' % '"P 1"')
            fileobj.write('_symmetry_int_tables_number       %d\n' % 1)
            fileobj.write('\n')

            fileobj.write('loop_\n')
            fileobj.write('  _symmetry_equiv_pos_as_xyz\n')
            fileobj.write("  'x, y, z'\n")
            fileobj.write('\n')

        fileobj.write('loop_\n')

        if format == 'mp':
            fileobj.write('  _atom_site_type_symbol\n')
            fileobj.write('  _atom_site_label\n')
            fileobj.write('   _atom_site_symmetry_multiplicity\n')
            fileobj.write('  _atom_site_fract_x\n')
            fileobj.write('  _atom_site_fract_y\n')
            fileobj.write('  _atom_site_fract_z\n')
            fileobj.write('  _atom_site_occupancy\n')
        else:
            fileobj.write('  _atom_site_label\n')
            fileobj.write('  _atom_site_occupancy\n')
            fileobj.write('  _atom_site_fract_x\n')
            fileobj.write('  _atom_site_fract_y\n')
            fileobj.write('  _atom_site_fract_z\n')
            fileobj.write('  _atom_site_thermal_displace_type\n')
            fileobj.write('  _atom_site_B_iso_or_equiv\n')
            fileobj.write('  _atom_site_type_symbol\n')

        #THIS IS THE IMPORTANT CHANGE
        scaled = atoms.get_scaled_positions(wrap=False)
        no = {}
        for i, atom in enumerate(atoms):
            symbol = atom.symbol
            if symbol in no:
                no[symbol] += 1
            else:
                no[symbol] = 1
            if format == 'mp':
                fileobj.write(
                    '  %-2s  %4s  %4s  %7.5f  %7.5f  %7.5f  %6.1f\n' %
                    (symbol, symbol + str(no[symbol]), 1,
                     scaled[i][0], scaled[i][1], scaled[i][2], 1.0))
            else:
                fileobj.write(
                    '  %-8s %6.4f %7.5f  %7.5f  %7.5f  %4s  %6.3f  %s\n' %
                    ('%s%d' % (symbol, no[symbol]),
                     1.0,
                     scaled[i][0],
                     scaled[i][1],
                     scaled[i][2],
                     'Biso',
                     1.0,
                     symbol))
