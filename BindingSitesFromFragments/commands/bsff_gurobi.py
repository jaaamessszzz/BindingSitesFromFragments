#!/usr/bin/env python3

import os
from ..utils import pdb_check

def gurobi_setup(args, bsff_config_dict):
    """
    Build a fuzzball using the best scoring side chain interactions with the defined ligand
    Score fuzzballs with Rosetta's energy function and dump

    Usage:
      bsff gurobi_setup <user_defined_dir> <current_iteration_fuzzball_dir>

    Arguments:
      <user_defined_dir>                Path to project root directory
      <current_iteration_fuzzball_dir>  Directory containing fuzzballs for current iteration
    """
    from ..gurobi_scoring import score_with_gurobi

    # todo: look into making bb atoms virtual when generating FeatureReporter
    gurobi = score_with_gurobi(args['<user_defined_dir>'], args['<current_iteration_fuzzball_dir>'], config_dict=bsff_config_dict)
    gurobi.generate_feature_reporter_db()
    gurobi.consolidate_scores_better()


def mc_solve(args):
    """
    Produce binding motif solutions using a simulated annealing monte carlo protocol

    Usage:
      bsff mcsolve sge <user_defined_dir> <current_iteration_fuzzball_dir> <motif_size> [options]
      bsff mcsolve local <user_defined_dir> <fuzzball_pdb_path> <motif_size> [options]

    Arguments:
      <user_defined_dir>                        Path to project root directory
      <current_iteration_fuzzball_dir>          Directory containing fuzzballs for current iteration
      <fuzzball_pdb_path>                       Path to a single fuzzball PDB
      <motif_size>                              Target binding motif size

    Options:
      --block=<block_size>, -b <block_size>     If running on the cluster, number of trajectories per fuzzball (default: 1)
      --include_defined, -i                     Included defined residues
    """
    from ..solve import montecarlo_motif
    user_defined_dir = args['<user_defined_dir>']
    motif_size = int(args['<motif_size>'])
    include_defined = args['--include_defined']
    block_count = int(args['--block'])

    if args['sge']:
        task_id = 0
        fuzzball_dir = args['<current_iteration_fuzzball_dir>']
        fuzzballs = [fuzz for fuzz in pdb_check(fuzzball_dir)]

        if os.environ.get("SGE_TASK_ID"):
            if args['--block']:
                block_size = int(args['--block'])
                task_id = int((int(os.environ["SGE_TASK_ID"]) - 1) / block_size)
                block_count = (int(os.environ["SGE_TASK_ID"]) - 1) - (task_id * block_size)
            else:
                # sge_task_id - 1 is to allow for easy indexing of lists...
                task_id = int(os.environ["SGE_TASK_ID"]) - 1
        fuzzball_path = fuzzballs[task_id]

    if args['local']:
        fuzzball_path = args['<fuzzball_pdb_path>']

    montecarlo_motif(user_defined_dir, fuzzball_path, motif_size, block_count=block_count, include_defined=include_defined)
