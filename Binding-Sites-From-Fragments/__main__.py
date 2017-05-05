#!/usr/bin/env python3

"""
Generate hypothetical ligand binding pockets for a given ligand based on
experimentally observed protein-ligand interactions in the PDB.

The proposed protocol is as follows:
1.  Specify target ligand
2.  Partition target ligand into fragments
3.  Search Pubmed/PDB for small molecules containing fragments as substructures
4.  For each identified fragment-containing small molecule, download all PDBs 
    with molecule bound
5.  Filter for good quality PDBs
6.  Align all small molecules based on fragments
7.  Cluster interactions and identify representative contacts
8.  Mix and score representative interactions

I am going to focus on identifying hydrogen-bonding contacts as we hypothesize
that these interactions will be contributing the majority of binding energy as 
well as binding specificity.

Usage:
    bsff generate_fragments <ligand> [options]

Arguments:
    generate_fragments
        Generate fragments for a given compound

    <ligand>
        By default, this is the name of the target ligand. This can be changed
        using the [ligand_input_format] option

Options:
    -f --ligand_input_format <format>
        Use a different input format for the ligand. [CID|name|smiles]

"""

if __name__ == '__main__':
    import docopt
    import numpy as np
    import Bio
    from .fragments import Fragments

    args = docopt.docopt(__doc__)

    frag = Fragments()

    if args['generate_fragments']:
        frag.generate_fragements_from_ligand(args['<ligand>'])
