LJ_TYPE_FOR_ELEMENT = {
    'O': 'O',
    'OM': 'O',
    'OA': 'O',
    'OE': 'O',
    'OW': 'O',
    'N': 'N',
    'NT': 'N',
    'NL': 'N',
    'NR': 'N',
    'NZ': 'N',
    'NE': 'N',
    'C': 'C',
    'CH0': 'C',
    'CH1': 'C',
    'CH2': 'C',
    'CH3': 'C',
    'CH4': 'C',
    'CH2r': 'C',
    'CR1': 'C',
    'HC': 'H',
    'H': 'H',
    'DUM': '',
    'S': '',
    'CU1+': 'CU',
    'CU2+': 'CU',
    'FE': 'FE',
    'ZN2+': 'ZN',
    'MG2+': 'MG',
    'CA2+': 'CA',
    'P,SI': 'P',
    'AR': 'AR',
    'F': 'F',
    'CL': 'CL',
    'BR': 'BR',
    'CMet': 'C',
    'OMet': 'O',
    'NA+': 'NA',
    'CL-': 'CL',
    'CChl': 'C',
    'CLChl': 'CL',
    'HChl': 'G',
    'SDmso': 'S',
    'CDmso': 'C',
    'ODmso': 'O',
    'CCl4': 'C',
    'CLCl4': 'CL',
    'FTFE': '',
    'CTFE': '',
    'CHTFE': '',
    'OTFE': '',
    'CUrea': 'C',
    'OUrea': 'O',
    'NUrea': 'N',
    'CH3p': 'C',
    'I': 'I',
    'CLOpt': 'CL',
    'B': 'B',
    'SE': 'SE',
    'HS14': 'H',
    'CLAro': 'CL',
    'BROpt': 'BR',
    'OEOpt': 'O',
    'NOpt': 'N',
    'CAro': 'C',
    'CPos': 'C',
    'NPri': 'N',
    'NTer': 'N',
    'OAlc': 'O',
    'SI': 'SI',
    'P': 'P',
}

S = ('>', '')
CHARGE = ('<', '')
F = ('', '.3f')

# Columns definition are inclusive on both ends
PDB_ATOM_INDEX_FIELD = (7, 11)
PDB_ATOM_NAME_FIELD = (13, 16)
PDB_COORD_FIELDS = ((31, 38), (39, 46),
                    (47, 54),
                    )# Taken from ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf, page 194
PDB_ELEMENT_FIELD = (77, 78)
PDB_CHARGE = (79, 80)

HETATM_SPECS = (
    (1, 6) + S,  # Record name
    PDB_ATOM_INDEX_FIELD + S,  # Atom serial number
    PDB_ATOM_NAME_FIELD + S,  # Atom name
    #    (17, 17) + S, # Alternate location indicator
    (18, 20) + S,  # Residue name
    (22, 22) + S,  # Chain identifier
    (23, 26) + S,  # Residue sequence number
    #    (27, 27) + S, # Code for insertion of residues
    PDB_COORD_FIELDS[0] + F,  # Orthogonal coordinates for X
    PDB_COORD_FIELDS[1] + F,  # Orthogonal coordinates for Y
    PDB_COORD_FIELDS[2] + F,  # Orthogonal coordinates for Z
    (55, 60) + S,  # Occupancy
    (61, 66) + S,  # Temperature factor
    PDB_ELEMENT_FIELD + S,  # Element symbol; right-justified
    PDB_CHARGE + CHARGE,  # Charge on the atom
)


def PDB_FORMAT_STR(pdb_specs):
    all_fields = []

    for i, fields in enumerate(pdb_specs):
        if i == 0:
            diff = 0
        else:
            diff = (fields[0] - 1) - (pdb_specs[i - 1][1])
        if diff != 0:
            all_fields.append(' ' * diff)

        all_fields.append(
            '{' + '{i}:{left_or_right}{len}{formatter}'.format(
                i=i,
                len=(fields[1] - fields[0] + 1),
                left_or_right=fields[2],
                formatter=fields[3],
            ) + '}'
        )
    return ''.join(all_fields)


PDB_TEMPLATE = PDB_FORMAT_STR(HETATM_SPECS)
