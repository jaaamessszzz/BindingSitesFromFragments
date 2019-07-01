#!/usr/bin/env python3

"""
Benchmark recovery of fuzzball contacts in the context of a match

"""

import os
import re
import math
import pickle
from pprint import pprint

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt


import pyrosetta
from pyrosetta import rosetta

from ..utils import *
from ..design import generate_fuzzball_contact_rotamersets


def fuzzball_contact_recovery_bm(args):
    """
    Count recovered fuzzball contacts within the context of a match at different rotamer sampling levels

    Usage:
      bsff benchmark fuzzball_contact_recovery <user_defined_dir> <match_dir> <fuzzball_dir> [options]

    Arguments:
      <user_defined_dir>      Path to project root directory
      <match_dir>             Directory containing matches to benchmark
      <fuzzball_dir>          Directory containing fuzzballs for current match set

    Options:
      --dump_rotsets, -d                Dump a PDB for the RotamerSets generated for each extrachi/sample level
      --taskid=<taskid>,-t=<taskid>     Run benchmark for the <taskid> match in os.listdir(<match_dir>)

    """

    def perform_rotamer_recovery(match_pdb_path, dump_rotamerset_pdb=False):

        # Set up infrastructure
        rotamer_info_list = list()
        rotamer_stats = dict()  # Dict containing stats for all rotamers evaluated for viable rotamers

        # Get matched motif residue and ligand positions
        match_name = os.path.normpath(os.path.basename(match_pdb_path))
        conformer_name, fuzzball_resnums = find_conformer_and_constraint_resnums(match_name)

        # Get fuzzball index from match_pdb
        match_pdb_split = re.split('[-_]', match_name)
        fuzzball_index = match_pdb_split[match_pdb_split.index('fuzz') + 1]

        # Instantiate PyRosetta things
        my_options = [
            f"-extra_res_fa {os.path.join(user_defined_dir, 'Inputs', 'Rosetta_Inputs', f'{conformer_name}.params')}",
            "-preserve_header",
            "-extrachi_cutoff=0",
            "-mute core.chemical.ResidueType core.pack.task"
        ]
        pyrosetta.init(options=' '.join(my_options))
        sfxn = rosetta.core.scoring.get_score_function()
        sfxn_weights = sfxn.weights()

        # Load match into pose
        match_pose = rosetta.core.pose.Pose()
        rosetta.core.import_pose.pose_from_file(match_pose, match_pdb_path)
        ligand_resnum = match_pose.size()

        # --- Load viable scaffold positions and corresponding residue types --- #

        match_residue_map_temp = pickle.load(open(os.path.join(fuzzball_dir, f'fuzz_{fuzzball_index}-match_residue_map.pickle'), 'rb'))

        # Only consider non-CGP positions
        match_residue_map = dict()
        for position, residues in match_residue_map_temp.items():
            if match_pose.residue(int(position)).name1() not in ['C', 'G', 'P']:
                match_residue_map[int(position)] = residues

        # --- Generate rotamers at different sampling levels --- #

        for extrachi in range(1, 5):

            # Set extrachi options
            for chi_on in range(1, extrachi + 1):
                rosetta.basic.options.set_boolean_option(f'packing:ex{chi_on}', True)

            # Track
            rotamer_stats[extrachi] = dict()

            for sample_level in range(1, 8):

                # Set sample level options for active chi
                for chi_sample in range(1, extrachi + 1):
                    rosetta.basic.options.set_integer_option(f'packing:ex{chi_sample}:level', sample_level)

                # Sanity check
                print('\n\nCURRENT OPTIONS:')
                for chi_check in range(1, 5):
                    print(f'packing:ex{chi_check}', rosetta.basic.options.get_boolean_option(f'packing:ex{chi_check}'))
                    print(f'packing:ex{chi_check}:level', rosetta.basic.options.get_integer_option(f'packing:ex{chi_check}:level'))
                print('\n\n')

                # Generate Rotamers
                viable_rotamers, rotamer_trial_info = generate_fuzzball_contact_rotamersets(user_defined_dir, match_pdb_path, match_pose, sfxn, match_residue_map,
                                                                                            flag_special_rot=False, report_stats=True, RMSD_limit=1, dump_rotamerset_pdb=dump_rotamerset_pdb)
                rotamer_stats[extrachi][sample_level] = rotamer_trial_info

                # Match pose clone for rotamer scores
                match_pose_clone = match_pose.clone()
                sfxn(match_pose_clone)

                # Add information to a list of dicts (for making a dataframe later)
                for position in viable_rotamers:
                    for residue_type in viable_rotamers[position]:
                        for rotamer in viable_rotamers[position][residue_type]:

                            # Add rotamer to pose, score, get one body energies from energy graph
                            match_pose_clone.replace_residue(position, rotamer, False)
                            sfxn(match_pose_clone)
                            e_edge = match_pose_clone.energies().energy_graph()
                            current_edge = e_edge.find_energy_edge(position, ligand_resnum)

                            if current_edge is not None:
                                emap = current_edge.fill_energy_map()
                                bsff_score = sum([emap[rosetta.core.scoring.fa_rep] * 0.55,
                                                  emap[rosetta.core.scoring.fa_atr],
                                                  emap[rosetta.core.scoring.hbond_sc],
                                                  emap[rosetta.core.scoring.fa_sol],
                                                  emap[rosetta.core.scoring.fa_elec]
                                                  ])
                                total_score = current_edge.dot(sfxn_weights)
                            else:
                                print(f'No edge between {position} and {ligand_resnum}')
                                bsff_score = None
                                total_score = None

                            rotamer_fa_dun = float(match_pose_clone.energies().onebody_energies(position)[rosetta.core.scoring.fa_dun])

                            rotamer_data = {'extrachi': extrachi,
                                            'sample_level': sample_level,
                                            'position': position,
                                            'resname': residue_type,
                                            'fa_dun': rotamer_fa_dun,
                                            'bsff_score': bsff_score,
                                            'total_score': total_score,
                                            }

                            rotamer_info_list.append(rotamer_data)

        # if dump_rotamerset_pdb:
        #
        #     # --- Dump packer rotset (code adapted from design) --- #
        #
        #     # Ligand RS
        #     matched_ligand_rs = rosetta.core.select.residue_selector.ResidueIndexSelector(match_pose.size())
        #
        #     # Positions directly adjacent to ligand RS
        #     ligand_adjacent_rs = rosetta.core.select.residue_selector.NeighborhoodResidueSelector(matched_ligand_rs, 8, include_focus_in_subset=True)
        #
        #     # NATRO positions in pose RS
        #     natro_rs = rosetta.core.select.residue_selector.NotResidueSelector(ligand_adjacent_rs)
        #
        #     # Create and Populate Task Factory
        #
        #     task_factory = rosetta.core.pack.task.TaskFactory()
        #
        #     racaa = rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
        #     racaa.aas_to_keep('ADEFHIKLMNQRSTVWY')  # No CGP
        #     design_op = rosetta.core.pack.task.operation.OperateOnResidueSubset(racaa, ligand_adjacent_rs)
        #     task_factory.push_back(design_op)
        #
        #     natro_op = rosetta.core.pack.task.operation.OperateOnResidueSubset(
        #         rosetta.core.pack.task.operation.PreventRepackingRLT(), natro_rs)
        #     task_factory.push_back(natro_op)
        #
        #     fixed_ligand_op = rosetta.core.pack.task.operation.OperateOnResidueSubset(
        #         rosetta.core.pack.task.operation.PreventRepackingRLT(), matched_ligand_rs)
        #     task_factory.push_back(fixed_ligand_op)
        #
        #     # Extra rotamers
        #     extra_rotamers_op = rosetta.core.pack.task.operation.ExtraRotamersGeneric()
        #     extra_rotamers_op.ex1(True)
        #     extra_rotamers_op.ex2(True)
        #     extra_rotamers_op.ex1_sample_level(1)
        #     extra_rotamers_op.ex2_sample_level(1)
        #     task_factory.push_back(extra_rotamers_op)
        #
        #     # Create RotamerSets including fuzzball rotamers
        #
        #     design_packer_task = task_factory.create_task_and_apply_taskoperations(task_factory, match_pose)
        #
        #     sfxn.setup_for_packing(match_pose, design_packer_task.repacking_residues(), design_packer_task.designing_residues())
        #     packer_neighbor_graph = rosetta.core.pack.create_packer_graph(match_pose, sfxn, design_packer_task)
        #
        #     rotamer_sets = rosetta.core.pack.rotamer_set.RotamerSetsFactory.create_rotamer_sets(match_pose)
        #     rotamer_sets.set_task(design_packer_task)
        #     rotamer_sets.initialize_pose_for_rotsets_creation(match_pose)
        #     rotamer_sets.build_rotamers(match_pose, sfxn, packer_neighbor_graph)
        #
        #     rotamer_sets.dump_pdb(match_pose, f"{match_name.split('.')[0]}-packer_rotset.pdb")

        # Generate dataframe for vaiable rotamers at different sample levels
        df = pd.DataFrame(rotamer_info_list)
        df.to_csv(f"{match_name.split('.')[0]}-contact_recovery_benchmark.csv")

        # Keep track of sampling types, rotamers generated, and rotamers accepted
        with open(f"{match_name.split('.')[0]}-rotamer_statistics.pickle", 'wb') as rotamer_stats_pickle:
            pickle.dump(rotamer_stats, rotamer_stats_pickle)

    # Parse args
    match_directory = args['<match_dir>']
    user_defined_dir = args['<user_defined_dir>']
    fuzzball_dir = args['<fuzzball_dir>']
    dump_rotsets = args['--dump_rotsets']

    # Either run in series or in parallel
    if args['--taskid']:
        current_match = os.listdir(match_directory)[int(args['--taskid']) - 1]
        perform_rotamer_recovery(os.path.join(match_directory, current_match), dump_rotamerset_pdb=dump_rotsets)

    else:
        for match_pdb in pdb_check(match_directory):
            perform_rotamer_recovery(os.path.join(match_directory, match_pdb), dump_rotamerset_pdb=dump_rotsets)


def rotamer_recovery_figures_for_tanja(args):
    """
    Yeah...

    Usage:
      bsff benchmark contact_recovery_figs <results_dir> [options]

    Arguments:
      <results_dir>           Directory containing csvs generated by fuzzball_contact_recovery_bm()

    Options:
      --skip_figures, -s      Skip making figures
      --help, -h              Show this message.

    """

    results_dir = args['<results_dir>']
    template_dir = os.path.join(os.path.dirname(__file__), '../..', 'Additional_Files', 'templates')

    from jinja2 import Environment, FileSystemLoader

    # http://eosrei.net/articles/2015/11/latex-templates-python-and-jinja2-generate-pdfs
    latex_jinja_env = Environment(
        block_start_string='\BLOCK{',
        block_end_string='}',
        variable_start_string='\VAR{',
        variable_end_string='}',
        comment_start_string='\#{',
        comment_end_string='}',
        line_statement_prefix='%%',
        line_comment_prefix='%#',
        trim_blocks=True,
        autoescape=False,
        loader=FileSystemLoader(template_dir))

    template = latex_jinja_env.get_template('rotamer_recovery_bm.tex')

    # Get unique matches in results dir
    match_set = set([match[:match.rfind('-')] for match in os.listdir(results_dir)])
    print(match_set)

    for unique_match in match_set:
        if not args['--skip_figures']:
            # todo: add option for output dir
            figure_output_dir = generate_figures(results_dir, unique_match)

    # Pass jinja dict with match_names as key and list of figures as value
    jijna_figure_dict = {match: dict() for match in match_set}
    for key, value in jijna_figure_dict.items():
        key_figure_list = [fig for fig in os.listdir(figure_output_dir) if key in fig]
        jijna_figure_dict[key] = {re.split('[-.]', fig)[-2]: os.path.join(results_dir, fig) for fig in key_figure_list}

    pprint(jijna_figure_dict)

    with open(f'Complementary_RotamerSet-analysis.tex', 'w') as asdf:
        asdf.write(template.render(jijna_figure_dict=jijna_figure_dict))


def generate_figures(results_dir, unique_match, figure_output_dir='benchmark_figures'):
    """
    Make pretty figures using data generated during a fuzzball assembly iteration.

    :param results_dir: Directory containing contact_recovery_benchmark.csv and rotamer_statistics.pickle for <unique_match>.
    :param unique_match: Name of the nucleated match scaffold.
    """

    os.makedirs(figure_output_dir, exist_ok=True)
    viable_residue_df = pd.read_csv(os.path.join(results_dir, f'{unique_match}-contact_recovery_benchmark.csv'))
    all_rotamers_stats = pickle.load(open(os.path.join(results_dir, f'{unique_match}-rotamer_statistics.pickle'), 'rb'))

    # --- BSFF vs Actual REU --- #
    """
    Seaborn lmplot for each extrachi and sample level
    """
    if viable_residue_df.empty is False:

        sns_reu_comparison = sns.FacetGrid(viable_residue_df, col="extrachi", row="sample_level")
        sns_reu_comparison = sns_reu_comparison.map(plt.scatter, "bsff_score", "total_score", edgecolor="w")
        sns_reu_comparison.savefig(os.path.join(figure_output_dir, f'{unique_match}-reu_comparison.png'), format='png', dpi=300, bbox_inches='tight')

    # --- rotamer recovery for sample level + extrachi --- #
    """
    Heatmap for the time being...
    """
    if viable_residue_df.empty is False:

        for extrachi in all_rotamers_stats:
            for sample_level in all_rotamers_stats[extrachi]:

                total_trial = 0
                total_accepted = 0

                for position in all_rotamers_stats[extrachi][sample_level]:
                    for contact_residue in all_rotamers_stats[extrachi][sample_level][position]:
                        total_trial += all_rotamers_stats[extrachi][sample_level][position][contact_residue]['num_rotamers']
                        total_accepted += all_rotamers_stats[extrachi][sample_level][position][contact_residue][
                            'rotamers_accepted']

        rot_recovery_fig, rot_recovery_ax = plt.subplots()

        # Aggregate counts
        heatmap_counts = viable_residue_df.groupby(['extrachi', 'sample_level']).agg('count').pivot_table(index='extrachi',
                                                                                                          columns='sample_level',
                                                                                                          values='total_score')
        print(heatmap_counts)

        sns.heatmap(heatmap_counts,
                    cmap="Blues", cbar=False,
                    # norm=matplotlib.colors.LogNorm(vmin=heatmap_counts.min().min(), vmax=heatmap_counts.max().max()),
                    vmin=heatmap_counts.min().min(),
                    vmax=heatmap_counts.max().max(),
                    square=True,
                    annot=True, fmt='g', annot_kws={'color': 'black', 'fontsize': 'small'},
                    ax=rot_recovery_ax)
        rot_recovery_fig.savefig(os.path.join(figure_output_dir, f'{unique_match}-rot_recovery.png'), format='png', dpi=300, bbox_inches='tight')

    # --- fa_dun vs sample level + extrachi --- #
    """
    Look at fa_dun as extrachi/sampling increases
    """
    if viable_residue_df.empty is False:

        fa_dun_fig, fa_dun_ax = plt.subplots(figsize=(4, 12))
        viable_residue_df['extrachi_sample'] = [f'Extrachi {extrachi}, Sampling {sample}' for extrachi, sample in
                                                zip(viable_residue_df['extrachi'], viable_residue_df['sample_level'])]
        sns.violinplot(x='fa_dun', y='extrachi_sample', data=viable_residue_df, scale='width', inner='quartile',
                       ax=fa_dun_ax)
        fa_dun_fig.savefig(os.path.join(figure_output_dir, f'{unique_match}-fa_dun.png'), format='png', dpi=300, bbox_inches='tight')

    # --- Energy distribution of accepted residues --- #

    if viable_residue_df.empty is False:

        reu_fig, reu_ax = plt.subplots(figsize=(4, 12))
        viable_residue_df['extrachi_sample'] = [f'Extrachi {extrachi}, Sampling {sample}' for extrachi, sample in
                                                zip(viable_residue_df['extrachi'], viable_residue_df['sample_level'])]
        sns.violinplot(x='total_score', y='extrachi_sample', data=viable_residue_df, scale='width', inner='quartile',
                       ax=reu_ax)
        reu_fig.savefig(os.path.join(figure_output_dir, f'{unique_match}-reu_accepted_dist.png'), format='png', dpi=300, bbox_inches='tight')

    # --- RMSD vs energy (highlight Accepted/Rejected) --- #

    # Construct dataframe fro pickle dict
    df_dict_list = list()

    for extrachi in all_rotamers_stats:
        for sample_level in all_rotamers_stats[extrachi]:
            for position in all_rotamers_stats[extrachi][sample_level]:
                for contact_residue in all_rotamers_stats[extrachi][sample_level][position]:
                    for rotamer_info in all_rotamers_stats[extrachi][sample_level][position][contact_residue]['rotamer_info']:
                        row_dict = {'extrachi': extrachi,
                                    'sample_level': sample_level,
                                    'position': position,
                                    'contact_residue': contact_residue,
                                    'min_RMSD': min(rotamer_info[0]),
                                    'contact_REU': rotamer_info[2],
                                    'accepted': True if (min(rotamer_info[0]) < 1 and rotamer_info[2] < 5) else False,
                                    }

                        df_dict_list.append(row_dict)

    residue_rmsd_reu_df = pd.DataFrame(df_dict_list)
    select_accepted_group = True in residue_rmsd_reu_df['accepted']

    if viable_residue_df.empty is False:
        rmsd_reu_grid = sns.FacetGrid(residue_rmsd_reu_df.groupby('accepted').get_group(select_accepted_group), col="extrachi", row="sample_level")
        rmsd_reu_grid = rmsd_reu_grid.map(plt.scatter, "min_RMSD", "contact_REU", edgecolor="w")
    else:
        rmsd_reu_grid = sns.FacetGrid(residue_rmsd_reu_df, col="extrachi", row="sample_level", sharey=False)
        rmsd_reu_grid = rmsd_reu_grid.map(plt.hist, "min_RMSD", edgecolor="w")

    rmsd_reu_grid.savefig(os.path.join(figure_output_dir, f'{unique_match}-rmsd_reu_grid.png'), format='png', dpi=300, bbox_inches='tight')

    # --- FacetGrid Accepted/Rejected at each position --- #
    accepted_grid = sns.FacetGrid(residue_rmsd_reu_df, col="extrachi", row="sample_level", hue='accepted')
    accepted_grid = accepted_grid.map(plt.hist, "position", bins=list(residue_rmsd_reu_df['position'].unique()))
    accepted_grid.set(yticks=[10 ** i for i in range(1, 6)], yticklabels=[10 ** i for i in range(1, 6)])

    for x_idx, x_ax in enumerate(accepted_grid.axes):
        for y_idx, y_ax in enumerate(accepted_grid.axes[x_idx]):
            accepted_grid.axes[x_idx][y_idx].set_yscale('log')

    accepted_grid.add_legend()
    accepted_grid.savefig(os.path.join(figure_output_dir, f'{unique_match}-accepted_grid.png'), format='png', dpi=300, bbox_inches='tight')

    # --- How many contact modes are being fulfilled? --- #
    """
    I'll do this later... maybe...
    """

    return figure_output_dir
