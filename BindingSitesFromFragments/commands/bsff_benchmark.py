#!/usr/bin/env python3

"""
Functions for benchmarking various components of the protocol.

Usage:  bsff benchmark <command> [<args>...]

Implemented benchmarks:

    fragment_representation         AUC measurement for fragment quality
    recover_fuzzball_contacts       Count the number and type of fuzzball contacts recovered in a set of matches
    binding_interface_recovery      Evaluate designs with varying special_rot weights
    binding_interface_figures       HMMER figures
    jsd_figures                     Profile similarity and native sequence recovery figures
    contact_recovery_figs           Figures for Tanja...
    design_metrics
    seqlogos

Options:
    -h --help

"""

import os
import docopt

from ..utils import *
from ..benchmarks.rotamer_contact_recovery import *
from ..benchmarks.fragment_representation import *
from ..benchmarks.binding_interface_recovery import *
from ..benchmarks.binding_interface_figures import *

def bm_handler(argv):
    """
    Perform the specified benchmark.

    :return:
    """

    registered_benchmarks = ['fuzzball_contact_recovery',
                             'fragment_representation',
                             'binding_interface_recovery',
                             ]

    registered_figures = ['contact_recovery_figs',
                          'binding_interface_figures',
                          'jsd_figures',
                          'design_metrics',
                          'seqlogos']

    registered_commands = registered_figures + registered_benchmarks
    current_argv = argv[1:]

    if len(current_argv) == 0 or current_argv[0] not in registered_commands:
        docopt.docopt(__doc__)

    if current_argv[0] == 'fragment_representation':
        fragment_representation_bm(docopt.docopt(fragment_representation_bm.__doc__, argv=argv))

    if current_argv[0] == 'fuzzball_contact_recovery':
        fuzzball_contact_recovery_bm(docopt.docopt(fuzzball_contact_recovery_bm.__doc__, argv=argv))

    if current_argv[0] == 'binding_interface_recovery':
        binding_interface_recovery_bm(docopt.docopt(binding_interface_recovery_bm.__doc__, argv=argv))

    # --- Figures... --- #

    if current_argv[0] == 'contact_recovery_figs':
        rotamer_recovery_figures_for_tanja(docopt.docopt(rotamer_recovery_figures_for_tanja.__doc__, argv=argv))

    if current_argv[0] == 'binding_interface_figures':
        binding_interface_figures(docopt.docopt(binding_interface_figures.__doc__, argv=argv))

    if current_argv[0] == 'jsd_figures':
        jsd_figures(docopt.docopt(jsd_figures.__doc__, argv=argv))

    if current_argv[0] == 'design_metrics':
        design_metric_histograms(docopt.docopt(design_metric_histograms.__doc__, argv=argv))

    if current_argv[0] == 'seqlogos':
        alignment_sequence_logo(docopt.docopt(alignment_sequence_logo.__doc__, argv=argv))

