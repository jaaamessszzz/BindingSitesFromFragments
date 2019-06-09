#!/usr/bin/env python3
"""
Methods for benchmarking

Usage:
    bm fetch_ligands <ligand_csv>
    bm fragment <ligand_dir>

Arguments:

    fetch_ligands
        Generates .mol2 files as specified by an input .csv file

    fragment
        Uses eMolFrag to generate a non-redundant fragment pool from .mol2 files in <ligand_dir>

    <ligand_csv>
        .csv file specifying ligands (3-letter PDB codes)

    <ligand_dir>
        Directory containing .mol2 files for ligands

Options:

"""

import docopt
import os
import sys
import pprint
import re
import yaml
import shutil
import pandas as pd
import appdirs
import xmltodict
import urllib

from .fragments import Fragments
from .alignments import Fragment_Alignments
from .clustering import Cluster
from .motifs import Generate_Motif_Residues, Generate_Binding_Sites, Generate_Constraints
from .utils import *

def fetch_ligands(df):
    """

    :return:
    """

if __name__ == "__main__":
    args = docopt.docopt(__doc__)

    if args['fetch_ligands']:
        df = pd.read_csv(args[''])
        fetch_ligands()