#!/usr/bin/env python3

def search(args):
    """
    Search the PDB for proteins bound to molecules possessing defined fragments as substructures.

    This command will take your defined fragments and their corresponding PubChem substructure search results to find
    proteins in the PDB observed to bind these molecules. Search results are saved and used by the next step to either
    download the structures from the PDB FTP server or use a local copy of the PDB.

    Usage: bsff search <user_defined_dir> [options]

    Arguments:
      <user_defined_dir>      Path to project root directory

    Options:
      -p, --pdb                 Use the PDB smilesQuery REST API to find substructure containing ligands
    """
    from ..fragments import Fragments
    frag = Fragments(args['<user_defined_dir>'])

    if args['--pdb']:
        frag.substructure_search_through_pdb()
    else:
        frag.substructure_search_through_pubchem()