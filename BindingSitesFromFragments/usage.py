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
    bsff bind_everything <user_defined_dir> <motif_size> [-s <score_cutoff_option>]
    bsff generate_constraints <user_defined_dir> <score_cutoff> [-m]
    bsff pull_constraints <user_defined_dir> <score_cutoff> [-n <number_to_pull>]
    bsff generate_fragments <ligand> [options]
    bsff dump_single_poses <user_defined_dir>
    bsff gurobi <user_defined_dir>
    bsff classic <user_defined_dir> <motif_size>
    bsff magic <user_defined_dir>
    bsff derp <user_defined_dir> <conformer> <gurobi_solutions_csv>

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

    generate_constraints
        Generate constraint files for all binding site conformations under a given score cutoff
    
    pull_constraints
        Pull constraint files with scores above a defined cutoff into a new directory
    
    dump_single_poses
        Dump single poses
        
    gurobi
        Use Gurobi to solve for best binding motifs given a feature reporter SQLITE3 DB generated with Rosetta
        
    classic
        Do everything from search through prepare_motifs using default parameters. Follow with generate_constraints and
        pull_constraints.
    
    magic
        Literally do everything. Will eventually poop complete constraint files for best binding motifs. 
        
    <ligand>
        By default, this is the name of the target ligand. This can be changed
        using the [ligand_input_format] option
        
    <user_defined_dir>
        Directory defined by user containing PubChem search results
    
    <motif_size>
        Number of motif residues to use in a binding motif
        
    <score_cutoff>
        Absolute value of score cutoff to use to generate constraint files. WILL BE CONVERTED TO NEGATIVE NUMBER,
        it's a pain to enter negative numbers in the command line... :(
    
    <number_to_pull>
        Limit to number of constraint files to pull with --number_to_pull option
        
Options:
    -c --clusters <clusters>
        Set number of clusters
        
    -d --distance_cutoff <distance_cutoff>
        Distance cutoff in angstroms for residues to consider in clustering
        
    -f --ligand_input_format <format>
        Use a different input format for the ligand. [CID|name|smiles]
        
    -m --secondary_matching
        Use secondary matching for constraint files
    
    -n --number_to_pull
        Limit the number of binding site constraint files to pull
        
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
import shutil
import pandas as pd
import appdirs

from .fragments import Fragments
from .alignments import Fragment_Alignments
from .clustering import Cluster
from .motifs import Generate_Motif_Residues, Generate_Binding_Sites, Generate_Constraints
from .utils import *

def main():
    # Import config
    bsff_config_path = os.path.join(os.path.dirname(__file__), '..', 'Additional_Files', 'bsff_config.yml')
    with open(bsff_config_path, 'r') as config_stream:
        bsff_config_dict = yaml.load(config_stream, Loader=yaml.Loader)

    # Interpret command line args
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

    if args['generate_fragments']:
        frag = Fragments(working_directory)
        frag.generate_fragements_from_ligand(args['<ligand>'])

    if args['search']:
        frag = Fragments(working_directory)
        frag.search_for_fragment_containing_ligands()

    if args['align']:
        alignment_monstrosity(working_directory, args)

    if args['cluster']:
        cluster(working_directory, args)

    if args['generate_motifs']:
        # Generate motif residues for each ligand conformer
        motifs = Generate_Motif_Residues(args['<user_defined_dir>'], config_dict=bsff_config_dict)
        motifs.generate_motif_residues()

    if args['prepare_motifs']:
        motifs = Generate_Motif_Residues(args['<user_defined_dir>'], config_dict=bsff_config_dict)
        motifs.prepare_motifs_for_conformers()

    if args['bind_everything']:

        bind = Generate_Binding_Sites(args['<user_defined_dir>'], config_dict=bsff_config_dict)
        bind.calculate_energies_and_rank(motif_size=int(args['<motif_size>']))
        bind.generate_binding_site_constraints(score_cutoff=-float(args['--score_cutoff_option']) if args['--score_cutoff_option'] else -10)

    if args['generate_constraints']:
        bind = Generate_Binding_Sites(args['<user_defined_dir>'], config_dict=bsff_config_dict)
        bind.generate_binding_site_constraints(score_cutoff=-float(args['<score_cutoff>']), secondary_matching=args['--secondary_matching'])

    if args['pull_constraints']:

        df_unsorted = pd.read_csv(os.path.join(args['<user_defined_dir>'],
                                               'Complete_Matcher_Constraints',
                                               '{}_scores_df.csv'.format(os.path.basename(os.path.normpath(args['<user_defined_dir>'])))
                                               )
                                  )

        df = df_unsorted.sort_values('total_score', ascending=True)
        df.reset_index(inplace=True)

        score_cutoff = -float(args['<score_cutoff>'])

        source_dir = os.path.join(args['<user_defined_dir>'], 'Complete_Matcher_Constraints', 'Constraint_Files')
        destination_dir = os.path.join(args['<user_defined_dir>'], 'Pulled_Constraints')
        os.makedirs(destination_dir, exist_ok=False)

        number_cutoff = int(args['<number_to_pull>']) if args['--number_to_pull'] else df.shape[0]

        for index, row in df.iterrows():
            if row['total_score'] < score_cutoff and index < number_cutoff:
                cst_file_name = '{}.cst'.format(row['description'][:-5])
                shutil.copy(os.path.join(source_dir, cst_file_name), os.path.join(destination_dir, cst_file_name))

    if args['dump_single_poses']:
        # Generate single poses
        motifs = Generate_Motif_Residues(args['<user_defined_dir>'], config_dict=bsff_config_dict)
        motifs.single_pose_cluster_residue_dump()
        
    if args['gurobi']:
        from .gurobi_scoring import score_with_gurobi
        gurobi = score_with_gurobi(args['<user_defined_dir>'], config_dict=bsff_config_dict)
        # gurobi.generate_feature_reporter_db()
        gurobi.consolidate_scores_better()
        gurobi.do_gurobi_things()

    if args['classic']:
        # Search
        frag = Fragments(working_directory)
        frag.search_for_fragment_containing_ligands()

        # Align
        alignment_monstrosity(working_directory, args)

        # Cluster
        cluster(working_directory, args)

        # Generate motifs
        motifs = Generate_Motif_Residues(args['<user_defined_dir>'], config_dict=bsff_config_dict)
        motifs.generate_motif_residues()

        # Prepare motifs
        motifs.prepare_motifs_for_conformers()

        # Bind everything
        bind = Generate_Binding_Sites(args['<user_defined_dir>'], config_dict=bsff_config_dict)
        bind.calculate_energies_and_rank(int(args['<motif_size>']))
        bind.generate_binding_site_constraints(score_cutoff=float(args['--score_cutoff_option']) if args['--score_cutoff_option'] else -10)

    if args['magic']:
        from .gurobi_scoring import score_with_gurobi

        # Search
        frag = Fragments(working_directory)
        frag.search_for_fragment_containing_ligands()

        # Align
        alignment_monstrosity(working_directory, args)

        # Cluster
        cluster(working_directory, args)

        # Generate single poses
        motifs = Generate_Motif_Residues(args['<user_defined_dir>'], config_dict=bsff_config_dict)
        motifs.single_pose_cluster_residue_dump()

        # Gurobi
        gurobi = score_with_gurobi(args['<user_defined_dir>'], config_dict=bsff_config_dict)
        gurobi.generate_feature_reporter_db()
        gurobi.consolidate_scores_better()
        gurobi.do_gurobi_things()

    if args['derp']:
        generate_constraints = Generate_Constraints(args['<user_defined_dir>'])
        generate_constraints.import_res_idx_map()
        generate_constraints.BW_constraints_from_gurobi_solutions(args['<gurobi_solutions_csv>'], args['<conformer>'])

# todo: all of this (poorly implemented) logic should be in the alignment class...
# todo: write a small function for adding rejected pdbs to the text file...
def alignment_monstrosity(working_directory, args, rmsd_cutoff=0.5):
    """
    Consequences of not thinking ahead...
    :param args: 
    :param rmsd_cutoff: fragment alignment RMSD cutoff, anything higher gets rejected 
    :return: 
    """
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
                        print("\n\nProcessing {}".format(pdbid))
                        # Set things up! Get ligands from Ligand Expo
                        # This is to avoid downloading ligands when all PDBs have already been processed

                        align = Fragment_Alignments(working_directory, ligand, processed_PDBs_path)

                        # lazy try/except block to catch errors in retrieving pdbs... namely empty pdbs...
                        # I can't believe how unorganized everything is, damn.

                        # Fetches ideal pdbs and list of PDB records where ligand is bound (for now)
                        try:
                            align.fetch_records()

                        except Exception as e:
                            print('{}: {}'.format(pdbid, e))
                            with open(rejected_list_path, 'a+') as reject_list:
                                reject_list.write('{}\n'.format(pdbid))
                            print('REJECTED - something went wrong with fetching the idealized target from LigandExpo')
                            continue

                        # Extract HETATM and CONECT records for the target ligand
                        # ligand_records, ligand_chain = align.extract_atoms_and_connectivities(ligand, pdb)

                        reject = align.extract_ligand_records(pdb)
                        # reject = align.fetch_specific_ligand_record(pdb)

                        # Reject if no ligands with all atoms represented can be found for the given pdb/ligand combo
                        if reject:
                            with open(rejected_list_path, 'a+') as reject_list:
                                reject_list.write('{}\n'.format(pdbid))
                            print('REJECTED - no target ligands were fully represented in the PDB')


                        # Continue if PDB has not been processed, rejected, or excluded by the user
                        else:
                            # Mapping of fragment atoms to target ligand atoms
                            align.fragment_string = open(
                                os.path.join(fragment_pdb, '{}.pdb'.format(current_fragment))).read()
                            mapping_successful = align.fragment_target_mapping()

                            if not mapping_successful:
                                with open(rejected_list_path, 'a+') as reject_list:
                                    reject_list.write('{}\n'.format(pdbid))
                                print('REJECTED - failed atom mapping between target and reference fragment')
                                continue

                            # Determine translation vector and rotation matrix
                            trgt_atom_coords, frag_atom_coords, transformation_matrix = align.determine_rotation_and_translation(current_fragment=current_fragment)

                            # Apply transformation to protein_ligand complex if rmsd if below cutoff
                            rmsd = prody.calcRMSD(frag_atom_coords, prody.applyTransformation(transformation_matrix, trgt_atom_coords))
                            print('RMSD of target onto reference fragment:\t{}'.format(rmsd))
                            if rmsd < rmsd_cutoff:
                                align.apply_transformation(transformation_matrix)
                            else:
                                with open(rejected_list_path, 'a+') as reject_list:
                                    reject_list.write('{}\n'.format(pdbid))
                                print('REJECTED - high RMSD upon alignment to reference fragment')

                    else:
                        print('{} has been processed!'.format(pdb))

def cluster(working_directory, args):
    """
    Consequences of not thinking ahead...
    :return: 
    """
    for fragment in directory_check(os.path.join(args['<user_defined_dir>'], 'Transformed_Aligned_PDBs')):
        # processed_PDBs_dir, distance_cutoff, weights

        # Set weights
        weights = [int(a) for a in args['--weights'].split()] if args['--weights'] else [1, 1, 1, 1]
        # Set distance cutoff
        distance_cutoff = args['--distance_cutoff'] if args['--distance_cutoff'] else 5

        cluster = Cluster(fragment, distance_cutoff, weights)

        if len(cluster.pdb_object_list) > 2:
            cluster.cluster_scipy()

        if cluster.clusters is not None:
            cluster.generate_output_directories(args['<user_defined_dir>'], fragment)
            cluster.automated_cluster_selection()
