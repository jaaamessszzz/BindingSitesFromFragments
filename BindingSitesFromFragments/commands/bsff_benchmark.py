#!/usr/bin/env python3

"""
Functions for benchmarking various components of the protocol.

Usage:  bsff benchmark <command> [<args>...]

Implemented benchmarks:

    fragment_representation         AUC measurement for fragment quality
    recover_fuzzball_contacts       Count the number and type of fuzzball contacts recovered in a set of matches

    contact_recovery_figs           Figures for Tanja...

Options:
    -h --help

"""

import os
import docopt

from ..utils import *
from ..benchmarks.rotamer_contact_recovery import *
from ..benchmarks.fragment_representation import *

def bm_handler(argv):
    """
    Perform the specified benchmark.

    :return:
    """

    registered_benchmarks = ['fuzzball_contact_recovery',
                             'fragment_representation',
                             'contact_recovery_figs',
                             ]
    current_argv = argv[1:]

    if len(current_argv) == 0 or current_argv[0] not in registered_benchmarks:
        docopt.docopt(__doc__)

    if current_argv[0] == 'fragment_representation':
        fragment_representation_bm(docopt.docopt(fragment_representation_bm.__doc__, argv=argv))

    if current_argv[0] == 'fuzzball_contact_recovery':
        fuzzball_contact_recovery_bm(docopt.docopt(fuzzball_contact_recovery_bm.__doc__, argv=argv))

    # Figures...
    if current_argv[0] == 'contact_recovery_figs':
        rotamer_recovery_figures_for_tanja(docopt.docopt(rotamer_recovery_figures_for_tanja.__doc__, argv=argv))