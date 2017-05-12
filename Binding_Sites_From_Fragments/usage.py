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
    bsff search <user_defined_dir>
    bsff align <user_defined_dir>

Arguments:
    generate_fragments
        Generate fragments for a given compound
    
    search
        Search for fragment-containing ligands
        
    align
        Align fragments in fragment-containing small molecules
        
    <ligand>
        By default, this is the name of the target ligand. This can be changed
        using the [ligand_input_format] option
        
    <user_defined_dir>
        Directory defined by user containing PubChem search results

Options:
    -f --ligand_input_format <format>
        Use a different input format for the ligand. [CID|name|smiles]

"""
import docopt
import os
import sys
import pprint
from .fragments import Fragments
from .alignments import Alignments

def main():

    args = docopt.docopt(__doc__)

    if args['<user_defined_dir>']:
        working_directory = os.path.join(os.path.curdir, args['<user_defined_dir>'])
    else:
        working_directory = os.path.curdir

    frag = Fragments(working_directory)
    align = Alignments(working_directory)

    if args['generate_fragments']:
        frag.generate_fragements_from_ligand(args['<ligand>'])

    if args['search']:
        frag.search_for_fragment_containing_ligands()

    if args['align']:
        # For each fragment, align all fragment-containing ligands to fragment
        # Generate PDBs with aligned coordinate systems

        # Fragment_1, Fragment_2, ...
        for fragment in directory_check(os.path.join(working_directory, 'Fragment_PDB_Matches')):
            fragment_pdb = os.path.join(working_directory, 'Inputs')

            # Three-letter codes for fragment-containing compounds
            for fcc in directory_check(fragment):
                # Each PDB containing a fragment-containing compound
                for pdb in pdb_check(fcc):
                    ligand = os.path.basename(os.path.normpath(fcc))
                    # ligand_records = align.extract_atoms_and_connectivities(ligand, pdb)
                    ligand_records = align.extract_atoms_and_connectivities(ligand, pdb)
                    pprint.pprint(ligand_records)
                    sys.exit()



def directory_check(dir):
    for subdir in os.listdir(dir):
        path = os.path.join(dir, subdir)
        if os.path.isdir(path):
            yield path


def pdb_check(dir):
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if path.endswith('.pdb'):
            yield path
