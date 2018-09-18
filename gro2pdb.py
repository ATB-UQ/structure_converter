import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from structure_converter.lib import PDB_TEMPLATE, LJ_TYPE_FOR_ELEMENT


def add_res_atom_numbers(gro_lines):
    res_atom_counter = 1
    N = len(gro_lines)
    gro_lines[0]["res_atom_number"] = res_atom_counter
    res_num = gro_lines[0]["res_num"]
    for i, l in enumerate(gro_lines[1:], start=1):
        if l["res_num"] == res_num:
            res_atom_counter += 1
            gro_lines[i]["res_atom_number"] = res_atom_counter
        elif l["res_num"] == res_num + 1:
            res_num += 1
            res_atom_counter = 1
            gro_lines[i]["res_atom_number"] = res_atom_counter
        else:
            raise Exception("Residue numbers not sequential:\n  {}\n  {}".format(gro_lines[i - 1], gro_lines[i]))


def read_gro(gro_file):
    def nm_to_A(x):
        return float(x)*10.
    columns = (
        ("res_num", (None, 5, int)),
        ("res_name", (5, 9, str)),
        ("atom_name", (9, 15, str)),
        ("atom_num", (15, 20, int)),
        ("x", (20, 28, nm_to_A)),
        ("y", (28, 36, nm_to_A)),
        ("z", (36, 44, nm_to_A)),
        ("vx", (44, 52, float)),
        ("vy", (52, 60, float)),
        ("vz", (60, 68, float)),
    )

    def parse_gro_line(l):
        return dict([(name, formatter(l[start:end].strip())) for name, (start, end, formatter) in columns])

    with open(gro_file, "r") as fh:
        gro_lines = [parse_gro_line(l) for l in fh.read().splitlines() if not l.startswith(";") and len(l.split()) > 7]
    add_res_atom_numbers(gro_lines)
    return gro_lines


def read_itp(itp_file):
    columns = (
        ("atom_num", int),
        ("atom_type", str),
        ("res_num", int),
        ("res_id", str),
        ("atom_name", str),
        ("ch_group", int),
        ("charge", float),
        ("mass", float),
    )

    def parse_itp_line(l):
        line_dict = dict([(name, formatter(col)) for (name, formatter), col in zip(columns, l.split())])
        return line_dict["atom_num"], line_dict

    with open(itp_file, "r") as fh:
        itp_str = fh.read()
    atom_lines = dict(
        [parse_itp_line(l)
         for l in itp_str.split("[ atoms ]")[1].split("[ bonds ]")[0].splitlines()
         if not l.startswith(";") and len(l.split()) > 7
         ]
    )
    assert len(set([l["res_id"] for l in atom_lines.values()])) == 1, "ITP file contains multiple res_id values"
    return atom_lines[1]["res_id"], atom_lines


def read_itp_files(itp_files):
    return dict([read_itp(f) for f in itp_files])


def gen_pdb(gro_lines, itp_dict):
    def pdb_line(index, name, residue_name, residue_number, x, y, z, element):
        return PDB_TEMPLATE.format(
            'ATOM  ',
            index,
            name,
            residue_name,
            '',
            residue_number,
            x,
            y,
            z,
            '',
            '',
            element,
            '',
        )

    pdb_lines = []
    for l in gro_lines:
        pdb_lines.append(
            pdb_line(
                l["atom_num"],
                l["atom_name"],
                l["res_name"],
                l["res_num"],
                l["x"],
                l["y"],
                l["z"],
                LJ_TYPE_FOR_ELEMENT[ itp_dict[l["res_name"]][l["res_atom_number"]]["atom_type"] ],
            )
        )
    return pdb_lines


def main(gro_file, itp_files):
    gro_lines = read_gro(gro_file)
    itp_dict = read_itp_files(itp_files)
    pdb_lines = gen_pdb(gro_lines, itp_dict)
    print("\n".join(pdb_lines))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:])
