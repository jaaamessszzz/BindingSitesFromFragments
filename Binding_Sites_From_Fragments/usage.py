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
    bsff new <compound_name>
    bsff generate_fragments <ligand> [options]
    bsff search <user_defined_dir>
    bsff align <user_defined_dir>
    bsff cluster <user_defined_dir> [options]
    bsff generate_motifs <user_defined_dir> [options]
    bsff bind <user_defined_dir> [options]

Arguments:
    new
        Sets up a new project with template input files and instructions
        
    generate_fragments
        Generate fragments for a given compound
    
    search
        Search for fragment-containing ligands
        
    align
        Align fragments in fragment-containing small molecules
    
    cluster
        Identify representitive residue contacts with each aligned fragment ensemble
        
    generate_motifs
        Generate motif residues from clusters specified in Inputs/motif_clusters.yml
    
    bind
        Generate hypothetical binding sites based on motif residues
        
    <ligand>
        By default, this is the name of the target ligand. This can be changed
        using the [ligand_input_format] option
        
    <user_defined_dir>
        Directory defined by user containing PubChem search results

Options:
    -f --ligand_input_format <format>
        Use a different input format for the ligand. [CID|name|smiles]
    
    -c --clusters
        Set number of clusters
        
    -d --distance_cutoff
        Distance cutoff in angstroms for residues to consider in clustering
        
    -w --weights
        Comma separated values for representative vector weights
        
    

"""
import docopt
import os
import sys
import pprint
import yaml
from .fragments import Fragments
from .alignments import Alignments, Align_PDB
from .clustering import Cluster
from .motifs import Generate_Motif_Residues, Generate_Binding_Sites
from .utils import *

def main():

    args = docopt.docopt(__doc__)

    if args['<user_defined_dir>']:
        working_directory = os.path.join(os.path.curdir, args['<user_defined_dir>'])
    else:
        working_directory = os.path.curdir

    if args['new']:
        print('Starting a new project for {}'.format(args['<compound_name>']))

    if args['generate_fragments']:
        frag = Fragments(working_directory)
        frag.generate_fragements_from_ligand(args['<ligand>'])

    if args['search']:
        frag = Fragments(working_directory)
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

                    prepare_ligand = Alignments(working_directory, current_fragment, ligand, processed_PDBs_path)

                    # Each PDB containing a fragment-containing compound
                    for pdb in pdb_check(fcc):
                        pdbid = os.path.basename(os.path.normpath(pdb))

                        # Check if PDB has already been processed
                        rejected_list_path = os.path.join(processed_PDBs_path, 'rejected_PDBs.txt')
                        rejected_list = []
                        if os.path.exists(rejected_list_path):
                            with open(rejected_list_path, 'r') as rejected_PDBs:
                                rejected_list = [pdb.strip() for pdb in rejected_PDBs]

                        processed_dir = os.path.join(working_directory, 'Transformed_Aligned_PDBs', current_fragment)

                        if not processed_check(processed_dir, pdbid, rejected_list):
                            # Set things up! Get ligands from Ligand Expo
                            # This is to avoid downloading ligands when all PDBs have already been processed
                            prepare_ligand.fetch_records()

                            align = Align_PDB(prepare_ligand)
                            # Extract HETATM and CONECT records for the target ligand
                            # ligand_records, ligand_chain = align.extract_atoms_and_connectivities(ligand, pdb)
                            reject = align.fetch_specific_ligand_record(pdb)

                            # Reject if no ligands with all atoms represented can be found for the given pdb/ligand combo
                            if reject:
                                with open(rejected_list_path, 'a+') as reject_list:
                                    reject_list.write('{}\n'.format(pdbid))

                            # Continue if PDB has not been processed, rejected, or excluded by the user
                            else:
                                # Mapping of fragment atoms to target ligand atoms
                                align.fragment_path = os.path.join(fragment_pdb, '{}.pdb'.format(current_fragment))
                                mapping_successful = align.fragment_target_mapping()

                                if not mapping_successful:
                                    with open(rejected_list_path, 'a+') as reject_list:
                                        reject_list.write('{}\n'.format(pdbid))
                                    continue

                                # Determine translation vector and rotation matrix
                                transformation_matrix = align.determine_rotation_and_translation()
                                # print(transformation_matrix.getMatrix())

                                # Apply transformation to protein_ligand complex
                                align.apply_transformation(transformation_matrix)

                        else:
                            print('{} exists!'.format(pdb))
    if args['cluster']:
        for fragment in directory_check(os.path.join(args['<user_defined_dir>'], 'Transformed_Aligned_PDBs')):
            # processed_PDBs_dir, distance_cutoff, number_of_clusters, weights

            # Set weights
            weights = args['--weights'] if args['--weights'] else [1,1,1,1]
            # Set distance cutoff
            distance_cutoff = args['--distance_cutoff'] if args['--distance_cutoff'] else 4
            # Set number of clusters
            number_of_clusters = args['--clusters'] if args['--clusters'] else 6

            cluster = Cluster(fragment, distance_cutoff, number_of_clusters, weights)
            # cluster.cluster_sklearn_agglomerative()
            cluster.cluster_scipy()
            if cluster.clusters is not None:
                cluster.generate_output_directories(args['<user_defined_dir>'], fragment)

    if args['generate_motifs']:
        motif_cluster_yaml = yaml.load(open(os.path.join(args['<user_defined_dir>'], 'Inputs', 'motif_clusters.yml'), 'r'))
        motifs = Generate_Motif_Residues(os.path.join(args['<user_defined_dir>'], 'Cluster_Results'), motif_cluster_yaml)
        motifs.generate_motif_residues()

    if args['bind']:
        motif_residue_bins = yaml.load(open(os.path.join(args['<user_defined_dir>'], 'Inputs', 'motif_residue_bins.yml'), 'r'))
        hypothetical_binding_sites = yaml.load(open(os.path.join(args['<user_defined_dir>'], 'Inputs', 'hypothetical_binding_sites.yml'), 'r'))
        bind = Generate_Binding_Sites(args['<user_defined_dir>'], motif_residue_bins, hypothetical_binding_sites)
        bind.generate_binding_sites()
