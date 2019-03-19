#!/usr/bin/env python3

import os
from ..utils import pdb_check

def constraints(args):
    """
    Generate matcher constraints for composite binding sites found using Gurobi or simulated annealing Monte Carlo.

    Usage:
      bsff_clean constraints <user_defined_dir> <gurobi_solutions_csv_dir> [--iteration <iteration_count>]
        [--constraints-to-generate <target>] [--offset <number>] [--consolidate <fuzzball_dir>]
        [--tolerance <tolerance>] [--sample <samples>] [--greasy] [--json]

    Arguments:
      <user_defined_dir>            Path to project root directory
      <gurobi_solutions_csv_dir>    Path to directory containing binding site solutions (.csv files)

    Options:
      -c <fuzzball_dir>, --consolidate <fuzzball_dir>   Consolidate solutions if Monte Carlo approach was used to
                                                        generate binding site solutions
      -i <count>, --iteration <count>                   Match iteration, constraints generated for all solutions
                                                        Provide iteration number where fuzzballs can be found
      -t <tolerance>, --tolerance <tolerance>           Matcher samples per constraint
      -s <samples>, --sample <samples>                  Matcher tolerance (in degrees) per constraint
      -n <target>, --constraints-to-generate <target>   Number of constraint files to generate if not iterating
      -o <number>, --offset <number>                    Skip <number> best constraints by score
      -g, --greasy                                      +1 sampling for hydrophobic residues (ACFILMVWY)
      -j, --json                                        Dump constraints into a JSON file
    """
    from ..motifs import Generate_Constraints
    generate_constraints = Generate_Constraints(args['<user_defined_dir>'])

    tolerance = int(args['--tolerance']) if args['--tolerance'] else 5
    samples = int(args['--sample']) if args['--sample'] else 1
    constraints_to_generate = int(args['--constraints-to-generate']) if args['--constraints-to-generate'] else 10000
    offset = int(args['--offset']) if args['--offset'] else 0
    iteration = args['--iteration']

    generate_constraints.conventional_constraints_from_gurobi_solutions(args['<gurobi_solutions_csv_dir>'],
                                                                        constraints_to_generate=constraints_to_generate,
                                                                        offset=offset,
                                                                        iteration=iteration,
                                                                        angle_dihedral_tolerance=tolerance,
                                                                        angle_dihedral_sample_number=samples,
                                                                        greasy_sampling=args['--greasy'],
                                                                        json_output=args['--json'],
                                                                        consolidate_solutions=args['--consolidate'])