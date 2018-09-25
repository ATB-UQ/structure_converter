### Convert GRO format to PDB
```
usage: gro2pdb.py [-h] -g GRO_FILE -f ITP_FILES [ITP_FILES ...] [-i]

optional arguments:

    -h, --help          show this help message and exit
    -g GRO_FILE, --gro_file GRO_FILE
                        GRO file.
    -f ITP_FILES [ITP_FILES ...], --itp_files ITP_FILES [ITP_FILES ...]
                        GROMACS ITP file.
    -i, --index_based_matching
                        Match GRO and ITP files only on atom index rather than
                        index in residue.
```