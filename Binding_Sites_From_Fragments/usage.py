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

    if args['generate_fragments']:
        frag.generate_fragements_from_ligand(args['<ligand>'])

    if args['search']:
        frag.search_for_fragment_containing_ligands()

    if args['align']:
        # For each fragment, align all fragment-containing ligands to fragment
        # Generate PDBs with aligned coordinate systems

        # Import list of ligands to exclude from processing
        # Three-letter ligand codes (UPPER CASE)
        # e.g. IMD, so much randomly bound IMD everywhere
        exclude_ligand_list = []
        exclude_txt = os.path.join(working_directory, 'Inputs', 'Exclude_ligands.txt')
        if os.path.exists(exclude_txt):
            with open(exclude_txt, 'r') as exlude_ligands:
                exclude_ligand_list = [lig.strip() for lig in exlude_ligands]

        # Fragment_1, Fragment_2, ...
        for fragment in directory_check(os.path.join(working_directory, 'Fragment_PDB_Matches')):
            fragment_pdb = os.path.join(working_directory, 'Inputs')
            current_fragment = os.path.basename(fragment)

            # Create directory for processed PDBs
            processed_PDBs_path = os.path.join(working_directory, 'Transformed_Aligned_PDBs', current_fragment)
            os.makedirs(processed_PDBs_path, exist_ok=True)

            # Three-letter codes for fragment-containing compounds
            for fcc in directory_check(fragment):
                ligand = os.path.basename(os.path.normpath(fcc))

                # Check if ligand is in exclusion list
                if ligand not in exclude_ligand_list:

                    align = Alignments(working_directory, current_fragment, ligand, processed_PDBs_path)

                    # Each PDB containing a fragment-containing compound
                    for pdb in pdb_check(fcc):
                        pdbid = os.path.basename(os.path.normpath(pdb)) # ASDF.pdb

                        # Check if PDB has already been processed
                        rejected_list_path = os.path.join(processed_PDBs_path, 'rejected_PDBs.txt')
                        rejected_list = []
                        if os.path.exists(rejected_list_path):
                            with open(rejected_list_path, 'r') as rejected_PDBs:
                                rejected_list = [pdb.strip() for pdb in rejected_PDBs]

                        processed_dir = os.path.join(working_directory, 'Transformed_Aligned_PDBs', current_fragment)

                        if not processed_check(processed_dir, pdbid, rejected_list):

                            # Extract HETATM and CONECT records for the target ligand
                            # ligand_records, ligand_chain = align.extract_atoms_and_connectivities(ligand, pdb)
                            reject, ligand_records, ligand_chain = align.fetch_specific_ligand_record(pdb)

                            # Reject if no ligands with all atoms represented can be found for the given pdb/ligand combo
                            if reject:
                                with open(rejected_list_path, 'a+') as reject_list:
                                    reject_list.write('{}\n'.format(pdbid))

                            # Continue if PDB has not been processed, rejected, or excluded by the user
                            else:

                                # Mapping of fragment atoms to target ligand atoms
                                pdb_path = os.path.join(fragment_pdb, '{}.pdb'.format(current_fragment))
                                fragment_target_mapping, target_mol_PDB_Block = align.fragment_target_mapping(pdb_path, ligand_records)

                                # Determine translation vector and rotation matrix
<<<<<<< HEAD
                                transformation_matrix = align.determine_rotation_and_translation(fragment_target_mapping, pdb_path, target_mol_PDB_Block)
                                # print(transformation_matrix.getMatrix())
=======
                                transformation_matrix = align.determine_rotation_and_translation(fragment_target_mapping, pdb_path, ligand_records)
                                print(transformation_matrix.getMatrix())
>>>>>>> a2142eb5441dfd7b5d92f650c296a556012ada27

                                # Apply transformation to protein_ligand complex
                                align.apply_transformation(transformation_matrix, pdb, ligand, ligand_chain)

                        else:
                            print('{} exists!'.format(pdb))
            # sys.exit()

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

def processed_check(processed_dir, pdb, rejected_list):
    pdbid = pdb.split('.')[0].upper()

    in_reject_list = any([pdb == reject for reject in rejected_list])
    already_processed = os.path.exists(os.path.join(processed_dir, '{}_processed.pdb'.format(pdbid)))

    return any([in_reject_list, already_processed])
