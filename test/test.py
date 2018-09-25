import sys
from os.path import dirname, abspath, join
TEST_DIR = dirname(abspath(__file__))
sys.path.insert(0, dirname(dirname(TEST_DIR)))

from structure_converter.gro2pdb import run


def peptide_test():
    pdb_lines = run(
        join(TEST_DIR, "peptide.gro"),
        [join(TEST_DIR, "peptide.itp")],
        match_on_index=True,
    )
    print("\n".join(pdb_lines))


def test():
    pdb_lines = run(
        join(TEST_DIR, "test.gro"),
        [
            join(TEST_DIR, "GRP.itp"),
            join(TEST_DIR, "IPS.itp"),
            join(TEST_DIR, "CBP.itp"),
        ],
        match_on_index=False,
    )
    print("\n".join(pdb_lines))


if __name__=="__main__":
    test()
    peptide_test()