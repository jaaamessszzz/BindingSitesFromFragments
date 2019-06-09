#!/usr/bin/env python3

import os
import re
import prody
from pprint import pprint

from ..solve import mc_crystallize, write_solutions_to_pdb
from ..motifs import Generate_Fuzzball_using_PyRosetta
from ..utils import pdb_check, find_conformer_and_constraint_resnums


def assemble(args):
    """
    Build a fuzzball using the best scoring side chain interactions with the defined ligand

    Usage:
      bsff assemble <user_defined_dir> [options]
      bsff assemble <user_defined_dir> iteration <selected_match_dir> [options]
      bsff assemble <user_defined_dir> existing <existing_complex_path> [options]

    Arguments:
      <user_defined_dir>              Path to project root directory
      existing                        Generate a fuzzball for an existing protein-ligand complex (e.g. from the PDB)
      <existing_complex_path>         Path to an existing pritein-ligand complex
      iteration                       Add motif residues from a previous iteration
      <selected_match_dir>            Path to directory containing selected matches using `prepare_match_iterations.py`

    Options:
      --previous_iteration=<dir>      Directory containing previous iteration fuzzballs
      --current_iteration=<dir>       Name for directory for current iteration fuzzballs
      --add_user_defined_motifs       Add any motif residues defined under `Inputs/Defined_Interactions` to the fuzzball
      --fuzzball_limit=<fuzz_limit>   Limit to the number of motif residues to be added to the fuzzball
      --hbond_limit=<hb_limit>        Limit to the number of hydrogen bonding residues to be added for each hydrogen
                                      bond donor/acceptor on the ligand
      -i=<index>, --index=<index>     Only generate fuzzballs with specified index
    """
    # Dump current iteration of fuzzballs into own directory
    fuzzball_dir = os.path.join(args['<user_defined_dir>'], 'Fuzzballs')
    iteration = len([x for x in os.listdir(fuzzball_dir) if os.path.isdir(os.path.join(fuzzball_dir, x))])

    if args['--current_iteration']:
        current_iteration_dir = os.path.join(fuzzball_dir, args['--current_iteration'])
    else:
        current_iteration_dir = os.path.join(fuzzball_dir, f'Iteration-{iteration}')

    # Make directory for Fuzzballs
    os.makedirs(current_iteration_dir, exist_ok=True)
    derp = Generate_Fuzzball_using_PyRosetta(args['<user_defined_dir>'], current_iteration_dir)

    # Map fuzzball iteration composition to struct_id
    if args['iteration']:
        """
        20190320 - Fuzzballs will be generated for each individual match to maximize the number of motif residues that 
        can be accommodated by a scaffold. 
        """

        # Find previous iteration Fuzzball dir
        previous_iteration_dirname = args['--previous_iteration'] if args['--previous_iteration'] else f'Iteration-{iteration - 1}'

        previous_iteration = os.path.join(fuzzball_dir, previous_iteration_dirname)
        if not os.path.exists(previous_iteration):
            raise Exception(f'Previous iteration directory {previous_iteration} does not exist!')

        # Get list of matches from <selected_match_dir>
        selected_matches = [match for match in pdb_check(args['<selected_match_dir>'])]
        if len(selected_matches) == 0:
            raise Exception('No matches were found in your selected_match_dir!')

        # Iteration stuff is too messy to turn into a function, have this jankiness instead
        match_indicies = [int(args['--index']) - 1] if args['--index'] else range(len(selected_matches))

        for match_index in match_indicies:

            def indicies_from_filename(file):
                """[2, 3, 4] from 38E_0001-1_2_3_4"""
                split_filename = file.split('-')[1].split('_')
                return [] if len(split_filename) == 1 else split_filename[1:]

            match_pdb = selected_matches[match_index]
            match_pdb_base = os.path.basename(os.path.normpath(match_pdb))

            # Get conformer and indicies for current constraint
            fuzzball_conformer, fuzzball_constraint_indicies = find_conformer_and_constraint_resnums(match_pdb_base)
            fuzzball_index = match_pdb_base.split('_')[-2]
            fuzzball_constraint = f'{fuzzball_conformer}-1_{"_".join([str(a) for a in fuzzball_constraint_indicies])}'

            # Get fuzzball from previous iteration
            previous_iteration_fuzzballs = [fuzz for fuzz in pdb_check(previous_iteration, base_only=True) if f'fuzz_{fuzzball_index}' in re.split('[.-]', fuzz)]

            pprint([fuzz.split('-') for fuzz in pdb_check(previous_iteration, base_only=True)])
            print(fuzzball_index)
            print(fuzzball_conformer)
            print(fuzzball_constraint_indicies)
            print(fuzzball_constraint)
            print(previous_iteration_fuzzballs)
            print(" ".join([str(a) for a in fuzzball_constraint_indicies]))

            assert len(previous_iteration_fuzzballs) == 1
            previous_iteration_fuzzball_file = previous_iteration_fuzzballs[0].split('.')[0]

            # Pull residues from previous iteration fuzzball
            previous_iteration_fuzzball = prody.loadAtoms(os.path.join(previous_iteration, f'{previous_iteration_fuzzball_file}.ag.npz'))
            current_defined_residues = previous_iteration_fuzzball.select(f'resnum {" ".join([str(a) for a in fuzzball_constraint_indicies])}')

            # Add residues to current fuzzball as defined residues
            iteration_dict = {'minimum_definition': current_defined_residues, 'match_path': match_pdb}
            fuzzball_path = derp.assemble_fuzzball(fuzzball_conformer, iteration=iteration_dict, iteration_name=fuzzball_constraint,
                                                   iteration_index=iteration, fuzzball_index=match_index)

            # Add fuzzball to FeatureReporter list
            with open(os.path.join(current_iteration_dir, 'fuzzball_list.txt'), 'a') as fuzzy_fuzz:
                fuzzy_fuzz.write(f'{fuzzball_path}\n')

            with open(os.path.join(current_iteration_dir, 'fuzzball_match_mapping.csv'), 'a') as fuzzy_fuzz:
                fuzzy_fuzz.write(f'{match_index},{os.path.basename(os.path.normpath(match_pdb))}\n')

    else:
        with open(os.path.join(current_iteration_dir, 'fuzzball_list.txt'), 'w') as fuzzy_fuzz:
            for conformer in pdb_check(derp.rosetta_inputs, base_only=True, conformer_check=True):
                conformer_name = conformer.split('.')[0]
                fuzzball_path = derp.assemble_fuzzball(conformer_name)
                fuzzy_fuzz.write(f'{fuzzball_path}\n')


def crystallize(args):
    """
    Build a fuzzball for an nucleated match

    Usage:
      bsff_clean crystallize <user_defined_dir> <match_dir> <match_fuzzball_dir> [options]

    Arguments:
      <user_defined_dir>        Path to project root directory
      <match_dir>               Path to a nucleated match
      <match_fuzzball_dir>      Path to fuzzball for current match

    Options:
      -b=<block_size>, --block=<block_size>                 If running on the cluster, number of trajectories per fuzzball
      -i=<positions_to_init>, --init=<positions_to_init>    Number of positions to initialize before MC
      -m=<match_index>, --match=<match_index>               Solve for a specific match
    """
    pprint(args)

    user_defined_dir = args['<user_defined_dir>']
    match_dir = args['<match_dir>']
    match_fuzzball_dir = args['<match_fuzzball_dir>']
    positions_to_init = args['--init'] if args['--init'] else 2

    # Ripped from mc_solve
    # todo: configure for both cluster and user-defined inputs
    task_id = 0
    block_count = None
    if os.environ.get("SGE_TASK_ID"):
        if args['--block']:
            block_size = int(args['--block'])
            task_id = int((int(os.environ["SGE_TASK_ID"]) - 1) / block_size)
            block_count = (int(os.environ["SGE_TASK_ID"]) - 1) - (task_id * block_size)
        else:
            # sge_task_id - 1 is to allow for easy indexing of lists...
            task_id = int(os.environ["SGE_TASK_ID"]) - 1

        print(os.environ)
        print(task_id)
    else:
        task_id = int(args['--match'])

    # Get match based on task id
    current_match = os.listdir(match_dir)[task_id]
    match_path = os.path.join(match_dir, current_match)

    # Get match_residue_map json for current match
    fuzzball_index = int(re.split('[_-]', current_match)[-2])

    mc_crystallize(user_defined_dir, match_path, match_fuzzball_dir, fuzzball_index, positions_to_init=positions_to_init, block_count=block_count)