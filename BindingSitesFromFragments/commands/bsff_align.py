#!/usr/bin/env python3

import os

from ..alignments import Align_PDB_Factory

def align(args):
    """
    Identify fragment substructures and align onto defined reference fragments

    This step of the protocol takes the search results from the previous step and either downloads structures using the PDB
    FTP server or a local copy of the PDB.

    Usage: bsff align <user_defined_dir> [options]

    Arguments:
      <user_defined_dir>      Path to project root directory

    Options:
      --use_local_pdb_database=<path_to_database>, -d=<path_to_database>    Path to root pdb/ directory of local copy of PDB

    """
    working_directory = args['<user_defined_dir>']
    use_local_pdb_database = args['--use_local_pdb_database']

    # Verify local database exists
    if use_local_pdb_database and not os.path.exists(use_local_pdb_database):
        raise Exception(f'Local PDB root directory {use_local_pdb_database} does not exist!')

    perform_alignments = Align_PDB_Factory(working_directory)
    perform_alignments.alignment_monstrosity(use_local_pdb_database=use_local_pdb_database)