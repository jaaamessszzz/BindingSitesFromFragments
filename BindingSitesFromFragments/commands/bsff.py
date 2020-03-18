#!/usr/bin/env python3
"""
Generate hypothetical ligand binding pockets for a given ligand based on experimentally observed protein-ligand
interactions in the PDB.

Usage:  bsff <command> [<args>...]

Commands listed in protocol order:

    new             Generate a new project
    search          Search the PDB for proteins bound compounds in PubChem search results
    align           Identify fragment substructures and align onto defined reference fragments
    cluster         Cluster side-chain interactions with defined fragments
    score           Use Rosetta to score interactions between side-chains and defined ligand
    assemble        Build a fuzzball using the best scoring side chain interactions with the defined ligand
    gurobi_setup    Solve for composite binding sites
    constraints     Generate RosettaMatch constraints for binding site solutions
    crystallize     Find full binding site definitions from a nucleated match

Options:
    -h --help
"""

__author__ = 'James Lucas'

import os
import sys
import yaml
from docopt import docopt


def main():
    # Import config
    bsff_config_path = os.path.join(os.path.dirname(__file__), '../..', 'Additional_Files', 'bsff_config.yml')
    with open(bsff_config_path, 'r') as config_stream:
        bsff_config_dict = yaml.load(config_stream, Loader=yaml.Loader)

    # Interpret command line args
    argv = sys.argv[1:]
    registered_commands = ['new',
                           'search',
                           'align',
                           'cluster',
                           'assemble',
                           'gurobi_setup',
                           'constraints',
                           'mcsolve',
                           'crystallize',
                           'benchmark',
                           'design'
                           ]

    if len(argv) == 0 or argv[0] not in registered_commands:
        args = docopt(__doc__)

    if argv[0] == 'new':
        from .bsff_new import new_project
        new_project(docopt(new_project.__doc__, argv=argv))

    if argv[0] == 'search':
        from .bsff_search import search
        search(docopt(search.__doc__, argv=argv))

    if argv[0] == 'align':
        from .bsff_align import align
        align(docopt(align.__doc__, argv=argv))

    if argv[0] == 'cluster':
        from .bsff_cluster import clustering
        clustering(docopt(clustering.__doc__, argv=argv))

    if argv[0] == 'assemble':
        from .bsff_assemble import assemble
        assemble(docopt(assemble.__doc__, argv=argv))

    if argv[0] == 'gurobi_setup':
        from .bsff_gurobi import gurobi_setup
        gurobi_setup(docopt(gurobi_setup.__doc__, argv=argv), bsff_config_dict=bsff_config_dict)

    if argv[0] == 'mcsolve':
        from .bsff_gurobi import mc_solve
        mc_solve(docopt(mc_solve.__doc__, argv=argv))

    if argv[0] == 'constraints':
        from .bsff_constraints import constraints
        constraints(docopt(constraints.__doc__, argv=argv))

    if argv[0] == 'crystallize':
        from .bsff_assemble import crystallize
        crystallize(docopt(crystallize.__doc__, argv=argv))

    if argv[0] == 'design':
        from .bsff_design import design
        design(docopt(design.__doc__, argv=argv))

    # If benchmark, pass things off to benchmark manager
    if argv[0] == 'benchmark':
        from .bsff_benchmark import bm_handler
        bm_handler(argv)


