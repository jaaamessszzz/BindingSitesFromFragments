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
    bsff dump_single_poses <user_defined_dir>
    bsff gurobi <user_defined_dir>
    bsff classic_gurobi_constraints <user_defined_dir> [iteration] <gurobi_solutions_csv_dir> [-t <tolerance>] [-s <samples>] [-g] [-j]
    bsff score <user_defined_dir>
    bsff assemble <user_defined_dir>

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

    -f --ligand_input_format <format>
        Use a different input format for the ligand. [CID|name|smiles]

    -g --greasy_sampling
        Increase the tolerance and sample number for greasy residues (ACFILMVWY)

    -j --json
        Generate a JSON file containing relevant constraint blocks in lieu of complete constraint files
        
    -m --secondary_matching
        Use secondary matching for constraint files
    
    -n --number_to_pull
        Limit the number of binding site constraint files to pull
        
    -r --rosetta_scores
        Calculate Rosetta score terms for unique residue-ligand interactions during the motif residue generation step

    -s --samples <samples>
        Set sample number for matcher constraint files

    -t --tolerances <tolerance> <samples>
        Set tolerances for matcher constraint files

    -w --weights <weights>
        Comma separated values for representative vector weights

    -x --score_cutoff_option <score_cutoff_option>
        Score cutoff to use for generating constraint files
"""

# --- Out-dated options --- #
# bsff generate_motifs <user_defined_dir> [options]
# bsff prepare_motifs <user_defined_dir> [options]
# bsff bind_everything <user_defined_dir> <motif_size> [-x <score_cutoff_option>]
# bsff generate_constraints <user_defined_dir> <score_cutoff> [-m]
# bsff pull_constraints <user_defined_dir> <score_cutoff> [-n <number_to_pull>]
# bsff generate_fragments <ligand> [options]
# bsff classic < user_defined_dir > < motif_size >
# bsff magic < user_defined_dir >
# bsff BW_gurobi_constraints < user_defined_dir > < conformer > < gurobi_solutions_csv >
# bsff derp < user_defined_dir > < conformer > < gurobi_solutions_csv >

__author__ = 'James Lucas'

import docopt
import yaml
import pandas as pd
import xmltodict
import urllib
import json

from .fragments import Fragments
from .alignments import Align_PDB_Factory
from .clustering import Cluster
from .motifs import Generate_Motif_Residues, Generate_Constraints, Generate_Fuzzball_using_PyRosetta
from .utils import *


def main():
    # Import config
    bsff_config_path = os.path.join(os.path.dirname(__file__), '..', 'Additional_Files', 'bsff_config.yml')
    with open(bsff_config_path, 'r') as config_stream:
        bsff_config_dict = yaml.load(config_stream, Loader=yaml.Loader)

    # Interpret command line args
    args = docopt.docopt(__doc__)

    # Get ligand three-letter code
    # ligand_code = args['<user_defined_dir>'][:3]

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
        with open(os.path.join(project_root_dir, 'Inputs', 'Fragment_Inputs', 'Fragment_inputs.csv'), 'w') as fragment_inputs:
            fragment_inputs.write('Fragment, SMILES_fragment')

    # if args['generate_fragments']:
    #     frag = Fragments(working_directory)
    #     frag.generate_fragements_from_ligand(args['<ligand>'])

    # todo: merge seach and align functions
    """
    20181222 - The old search implementation generated an intersection of the LigandExpo PDB ligand list and PubChem
    search results for each fragment, followed by an impossibly slow serial download of all relevant PDBs using PyPDB. 
    This was impossibly slow since it uses the PDB downloadFile.do tool. 
    
    Search will now dump all search results into a JSON that can be referenced by the alignment step. 
    Align will now either use the wwPDB FTP server to download structures or use locally available structures.
    """

    if args['search']:
        frag = Fragments(working_directory)
        frag.search_for_fragment_containing_ligands()

    if args['align']:
        use_cluster_pdb_database = False
        perform_alignments = Align_PDB_Factory(working_directory)
        perform_alignments.alignment_monstrosity(use_cluster_pdb_database=use_cluster_pdb_database)

    if args['cluster']:
        for fragment in directory_check(os.path.join(args['<user_defined_dir>'], 'Transformed_Aligned_PDBs')):
            # processed_PDBs_dir, distance_cutoff, weights

            # Set weights
            weights = [int(a) for a in args['--weights'].split()] if args['--weights'] else [1, 1, 1, 1]

            cluster = Cluster(fragment, weights)

            if len(cluster.pdb_object_list) > 2:
                cluster.cluster_scipy()

            if cluster.clusters is not None:
                cluster.generate_output_directories(args['<user_defined_dir>'], fragment)
                cluster.automated_cluster_selection()

    if args['dump_single_poses']:
        # Generate single poses
        motifs = Generate_Motif_Residues(args['<user_defined_dir>'], config_dict=bsff_config_dict)
        motifs.generate_all_the_fuzzballs()

    if args['gurobi']:
        from .gurobi_scoring import score_with_gurobi
        gurobi = score_with_gurobi(args['<user_defined_dir>'], config_dict=bsff_config_dict)
        gurobi.generate_feature_reporter_db()
        gurobi.consolidate_scores_better()
        # gurobi.do_gurobi_things()

    if args['classic_gurobi_constraints']:
        # todo: add option for number of constraints to generate
        # todo: add option for iterations
        generate_constraints = Generate_Constraints(args['<user_defined_dir>'])
        generate_constraints.import_res_idx_map()

        tolerance = int(args['--tolerances']) if args['--tolerances'] else 5
        samples = int(args['--samples']) if args['--samples'] else 1
        generate_constraints.conventional_constraints_from_gurobi_solutions(args['<gurobi_solutions_csv_dir>'],
                                                                            constraints_to_generate=5000,
                                                                            offset=0,
                                                                            iteration=args['iteration'],
                                                                            angle_dihedral_tolerance=tolerance,
                                                                            angle_dihedral_sample_number=samples,
                                                                            greasy_sampling=args['--greasy_sampling'],
                                                                            json_output=args['--json'])

    if args['score']:
        derp = Generate_Fuzzball_using_PyRosetta(args['<user_defined_dir>'])
        derp.score_motif_conformer_interactions()

    if args['assemble']:
        derp = Generate_Fuzzball_using_PyRosetta(args['<user_defined_dir>'])
        for conformer in pdb_check(os.path.join(args['<user_defined_dir>'], 'Inputs', 'Rosetta_Inputs'), base_only=True, conformer_check=True):
            conformer_name = conformer.split('.')[0]
            derp.assemble_defined_fuzzball(conformer_name)