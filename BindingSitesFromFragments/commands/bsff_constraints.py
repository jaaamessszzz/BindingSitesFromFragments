#!/usr/bin/env python3

import os
from ..utils import pdb_check

def constraints(args):
    """
    Generate matcher constraints for composite binding sites found using Gurobi or simulated annealing Monte Carlo.

    Usage:
      bsff_clean constraints <user_defined_dir> <gurobi_solutions_csv_dir> <fuzzball_dir> [--iteration]
        [--constraints-to-generate <target>] [--offset <number>] [--consolidate]
        [--tolerance <tolerance>] [--sample <samples>] [--greasy] [--json]

    Arguments:
      <user_defined_dir>            Path to project root directory
      <gurobi_solutions_csv_dir>    Path to directory containing binding site solutions (.csv files)
      <fuzzball_dir>                Path to directory containing fuzzballs for current constraints

    Options:
      -c, --consolidate                                 Consolidate solutions if Monte Carlo approach was used to
                                                        generate binding site solutions
      -i, --iteration                                   Match iteration, constraints generated for all solutions
                                                        Provide iteration number where fuzzballs can be found
      -t <tolerance>, --tolerance <tolerance>           Matcher samples per constraint
      -s <samples>, --sample <samples>                  Matcher tolerance (in degrees) per constraint
      -n <target>, --constraints-to-generate <target>   Number of constraint files to generate if not iterating
      -o <number>, --offset <number>                    Skip <number> best constraints by score
      -g, --greasy                                      +1 sampling for hydrophobic residues (ACFILMVWY)
      -j, --json                                        Dump constraints into a JSON file
    """
    from ..motifs import Generate_Constraints
    chemical_component_identifier = os.path.basename(os.path.normpath(args['<user_defined_dir>']))[:3]
    generate_constraints = Generate_Constraints(chemical_component_identifier)

    tolerance = int(args['--tolerance']) if args['--tolerance'] else 5
    samples = int(args['--sample']) if args['--sample'] else 1
    constraints_to_generate = int(args['--constraints-to-generate']) if args['--constraints-to-generate'] else 100000
    offset = int(args['--offset']) if args['--offset'] else 0

    generate_constraints.conventional_constraints_from_gurobi_solutions(args['<user_defined_dir>'],
                                                                        args['<gurobi_solutions_csv_dir>'],
                                                                        args['<fuzzball_dir>'],
                                                                        constraints_to_generate=constraints_to_generate,
                                                                        offset=offset,
                                                                        iteration=args['--iteration'],
                                                                        angle_dihedral_tolerance=tolerance,
                                                                        angle_dihedral_sample_number=samples,
                                                                        greasy_sampling=args['--greasy'],
                                                                        json_output=args['--json'],
                                                                        consolidate_solutions=args['--consolidate'])