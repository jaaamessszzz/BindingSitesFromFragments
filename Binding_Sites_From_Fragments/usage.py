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
    bsff search <user_defined_dir>
    bsff align <user_defined_dir>
    bsff cluster <user_defined_dir> [options]
    bsff generate_motifs <user_defined_dir> [options]
    bsff prepare_motifs <user_defined_dir> [options]
    bsff bind_everything <user_defined_dir> [options]
    bsff generate_constraints <user_defined_dir> <score_cutoff>
    bsff generate_fragments <ligand> [options]
    bsff bind_by_hand <user_defined_dir> [options]
    bsff derp

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
    
    prepare_motifs
        Prepare motifs for conformer binding site generation
    
    bind_everything
        EVERYTHING.
    
    generate_constraints
        Generate constraint files for all binding site conformations under a given score cutoff
        
    bind_by_hand
        Generate hypothetical binding sites based on motif residues bins defined by the user
    
    <ligand>
        By default, this is the name of the target ligand. This can be changed
        using the [ligand_input_format] option
        
    <user_defined_dir>
        Directory defined by user containing PubChem search results

    <score_cutoff>
        Absolute value of dcore cutoff to use to generate constraint files. WILL BE CONVERTED TO NEGATIVE NUMBER,
        it's a pain to enter negative numbers in the command line... :(
        
Options:
    -c --clusters <clusters>
        Set number of clusters
        
    -d --distance_cutoff <distance_cutoff>
        Distance cutoff in angstroms for residues to consider in clustering
        
    -f --ligand_input_format <format>
        Use a different input format for the ligand. [CID|name|smiles]
        
    -r --rosetta_scores
        Calculate Rosetta score terms for unique residue-ligand interactions during the motif residue generation step
    
    -s --score_cutoff_option <score_cutoff_option>
        Score cutoff to use for generating constraint files
        
    -w --weights <weights>
        Comma separated values for representative vector weights

"""
import docopt
import os
import sys
import pprint
import yaml
from .fragments import Fragments
from .alignments import Fragment_Alignments
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
        project_root_dir = os.path.join(os.getcwd(), args['<compound_name>'])
        os.makedirs(project_root_dir)
        os.makedirs(os.path.join(project_root_dir, 'Inputs'))
        os.makedirs(os.path.join(project_root_dir, 'Inputs', 'Rosetta_Inputs'))
        os.makedirs(os.path.join(project_root_dir, 'Inputs', 'Fragment_Inputs'))
        os.makedirs(os.path.join(project_root_dir, 'Inputs', 'User_Inputs'))
        with open(os.path.join(project_root_dir, 'Inputs', 'Fragment_Inputs', 'Fragment_Inputs.csv'), 'w') as fragment_inputs:
            fragment_inputs.write('Fragment, SMILES_fragment')

        with open(os.path.join(project_root_dir, 'Inputs', 'User_Inputs', 'Hypothetical_Binding_Sites.yml'), 'w') as Hypothetical_Binding_Sites:
            Hypothetical_Binding_Sites.write('\n'.join([
                '# Hypothetical binding sites for {}'.format(args['<compound_name>']),
                '# Create dictionaries with lists of motif bins that will be combinatorially',
                '# combined to generate unique binding sites.',
                '# These motif cluster names should be defined in Motif_Residue_Bins.yml',
                '# Example:',
                '# Binding_Site_1:',
                '#  - Ring_Sammich_Top',
                '#  - Ring_Sammich_Bottom',
                '#  - Carboxyl_Bidentate',
                '#  - Phenol_H-Bond']
            ))

        with open(os.path.join(project_root_dir, 'Inputs', 'User_Inputs', 'Motif_Clusters.yml'), 'w') as Motif_Clusters:
            Motif_Clusters.write('\n'.join([
                '# Clusters to use for generating motif residues',
                '# Create dictionaries with lists of cluster indices for each fragment.',
                '# Representative motif residues will be pulled from these clusters.',
                '# {}'.format(args['<compound_name>']),
                '# Example:',
                '# Fragment_1:',
                '#   - 42',
                '#   - 24',
                '#   - 40',
                '#   - 41',]
            ))

        with open(os.path.join(project_root_dir, 'Inputs', 'User_Inputs', 'Motif_Residue_Bins.yml'), 'w') as Motif_Residue_Bins:
            Motif_Residue_Bins.write('\n'.join([
                '# Residue bins for {}'.format(args['<compound_name>']),
                '# Define groups of residues around the ligand that make a unique type of contact',
                '# Create dictionaries with lists of motif residue indices for each bin',
                '# Example:',
                '# Phenol_H-Bond:',
                '#   - 14',
                '#   - 15',
                '#   - 20',
                '#   - 23',
                '#   - 22']
            ))

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
        exclude_txt = os.path.join(working_directory, 'Inputs', 'User_Inputs', 'Exclude_Ligands.txt')
        if os.path.exists(exclude_txt):
            with open(exclude_txt, 'r') as exlude_ligands:
                exclude_ligand_list = [lig.strip() for lig in exlude_ligands]

        # Fragment_1, Fragment_2, ...
        for fragment in directory_check(os.path.join(working_directory, 'Fragment_PDB_Matches')):
            fragment_pdb = os.path.join(working_directory, 'Inputs', 'Fragment_Inputs')
            current_fragment = os.path.basename(fragment)

            # Create directory for processed PDBs
            processed_PDBs_path = os.path.join(working_directory, 'Transformed_Aligned_PDBs', current_fragment)
            os.makedirs(processed_PDBs_path, exist_ok=True)

            # Three-letter codes for fragment-containing compounds
            for fcc in directory_check(fragment):
                ligand = os.path.basename(os.path.normpath(fcc))

                # Check if ligand is in exclusion list
                if ligand not in exclude_ligand_list:

                    align = Fragment_Alignments(working_directory, ligand, processed_PDBs_path)

                    # Each PDB containing a fragment-containing compound
                    for pdb in pdb_check(fcc):
                        pdbid = os.path.basename(os.path.normpath(pdb))

                        # Check if PDB has already been processed
                        rejected_list_path = os.path.join(processed_PDBs_path, 'Rejected_PDBs.txt')
                        rejected_list = []
                        if os.path.exists(rejected_list_path):
                            with open(rejected_list_path, 'r') as rejected_PDBs:
                                rejected_list = [pdb.strip() for pdb in rejected_PDBs]

                        processed_dir = os.path.join(working_directory, 'Transformed_Aligned_PDBs', current_fragment)

                        if not processed_check(processed_dir, pdbid, rejected_list):
                            # Set things up! Get ligands from Ligand Expo
                            # This is to avoid downloading ligands when all PDBs have already been processed

                            # lazy try/except block to catch errors in retrieving pdbs... namely empty pdbs...
                            # I can't believe how unorganized everything is, damn.
                            try:
                                align.fetch_records()

                            except Exception as e:
                                print('{}: {}'.format(pdbid,e))
                                with open(rejected_list_path, 'a+') as reject_list:
                                    reject_list.write('{}\n'.format(pdbid))
                                continue

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
                                align.fragment_string = open(os.path.join(fragment_pdb, '{}.pdb'.format(current_fragment))).read()
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
                            print('{} has been processed!'.format(pdb))
    if args['cluster']:
        for fragment in directory_check(os.path.join(args['<user_defined_dir>'], 'Transformed_Aligned_PDBs')):
            # processed_PDBs_dir, distance_cutoff, number_of_clusters, weights

            # Set weights
            weights = [int(a) for a in args['--weights'].split()] if args['--weights'] else [1,1,1,1]
            # Set distance cutoff
            distance_cutoff = args['--distance_cutoff'] if args['--distance_cutoff'] else 4
            # Set number of clusters
            number_of_clusters = args['--clusters'] if args['--clusters'] else 6

            cluster = Cluster(fragment, distance_cutoff, number_of_clusters, weights)

            if len(cluster.pdb_object_list) > 2:
                cluster.cluster_scipy()

            if cluster.clusters is not None:
                cluster.generate_output_directories(args['<user_defined_dir>'], fragment)

    if args['generate_motifs']:
        motif_cluster_yaml = yaml.load(open(os.path.join(args['<user_defined_dir>'], 'Inputs', 'User_Inputs', 'Motif_Clusters.yml'), 'r'))
        # Generate motif residues for each ligand conformer
        motifs = Generate_Motif_Residues(args['<user_defined_dir>'], motif_cluster_yaml)
        motifs.generate_motif_residues()

    if args['prepare_motifs']:
        motif_cluster_yaml = yaml.load(open(os.path.join(args['<user_defined_dir>'], 'Inputs', 'User_Inputs', 'Motif_Clusters.yml'), 'r'))
        motifs = Generate_Motif_Residues(args['<user_defined_dir>'], motif_cluster_yaml)
        motifs.prepare_motifs_for_conformers()

    if args['bind_by_hand']:
        motif_residue_bins = yaml.load(open(os.path.join(args['<user_defined_dir>'], 'Inputs', 'User_Inputs', 'Motif_Residue_Bins.yml'), 'r'))
        hypothetical_binding_sites = yaml.load(open(os.path.join(args['<user_defined_dir>'], 'Inputs', 'User_Inputs', 'Hypothetical_Binding_Sites.yml'), 'r'))
        bind = Generate_Binding_Sites(args['<user_defined_dir>'], motif_residue_bins, hypothetical_binding_sites)
        bind.generate_binding_sites_by_hand()

    if args['bind_everything']:
        bind = Generate_Binding_Sites(args['<user_defined_dir>'])
        bind.calculate_energies_and_rank()
        bind.generate_binding_site_constraints(score_cutoff=float(args['--score_cutoff_option']) if args['--score_cutoff_option'] else -15)

    if args['generate_constraints']:
        bind = Generate_Binding_Sites(args['<user_defined_dir>'])
        bind.generate_binding_site_constraints(score_cutoff=-float(args['<score_cutoff>']))

    if args['derp']:
        print('hi')