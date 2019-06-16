#!/usr/bin/env python3

import os
import sys
import copy
import json
import math
import pickle
from pprint import pprint

import prody
import pandas as pd

from .utils import *

import pyrosetta
from pyrosetta import rosetta


def generate_constrained_backrub_ensemble(raw_match_path, matcher_constraint):
    """
    Generate a backrub ensemble around the hypothetical binding site with matcher constraints applied
    :return:
    :param raw_match_path: path to raw match output
    :param matcher_constraint: constraint file used to find match
    """
    pass


def parse_matcher_remarks(match_path):
    """
    Parse matcher remarks to return match positions
    :param match_pose: pose with matcher remark header
    :return:
    """

    motif_resnums = list()
    with open(match_path, 'r') as match:
        for line in match:
            split_remark = line.split()

            # Only parse matcher remarks
            if split_remark[0] == 'REMARK':
                if all([split_remark[:4] == ['REMARK', '666', 'MATCH', 'TEMPLATE'], split_remark[7:9] == ['MATCH', 'MOTIF']]):
                    motif_resnums.append(int(split_remark[11]))

    return motif_resnums


def generate_fuzzball_contact_rotamersets(user_defined_dir, match_path, match_pose, sfxn, match_residue_map, flag_special_rot=True, custom_taskop=None, rotset_limit=100, RMSD_limit=1, dump_rotamerset_pdb=False, report_stats=False):
    """
    Generate rotamers that recapitulate observed fuzzball contacts for each position in a nucleated match

    :param flag_special_rot: If true, flag rotamers as SPECIAL_ROT variants
    :param custom_taskop: list of task operations to apply to the PackerTask used to generate rotamers

    :return:
    """

    sfxn_weights = sfxn.weights()
    conformer_resnum = match_pose.size()  # Assumes single ligand appended to end of sequence

    # --- Find and store viable rotamers --- #

    viable_rotamers = dict()
    rotamer_stats = dict()

    # Setting things up is going to mess up the match pose, so use a clone
    match_pose_clone = match_pose.clone()
    sfxn(match_pose_clone)

    # --- Transform match pose clone onto fuzzball conformer --- #
    """Required for contact coordsets to make sense"""

    # Get ligand from match, always last residue
    match_pose_size = match_pose_clone.size()
    match_ligand = match_pose_clone.residue(match_pose_size)
    conformer_name, motif_resnums = find_conformer_and_constraint_resnums(os.path.basename(os.path.normpath(match_path)))
    motif_and_ligand_resnums = motif_resnums + [conformer_resnum]

    # Keep track of match positions and compatible residue identites
    # match_residue_map = {position: dict() for position in range(1, match_pose.size())}  # Assumes one ligand appended to end of sequence

    # Import conformer from pose
    ligand_conformer_path = os.path.join(user_defined_dir, 'Inputs', 'Rosetta_Inputs', f'{conformer_name}.pdb')
    fuzzball_ligand_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(fuzzball_ligand_pose, ligand_conformer_path)
    fuzzball_ligand = fuzzball_ligand_pose.residue(1)

    # Calculate rotation/translation by hand using first three atoms of ligand
    mobile_match = rosetta.numeric.xyzTransform_double_t(match_ligand.xyz(1), match_ligand.xyz(2), match_ligand.xyz(3))
    mobile_match_inverse = mobile_match.inverse()
    target_fuzzball = rosetta.numeric.xyzTransform_double_t(fuzzball_ligand.xyz(1), fuzzball_ligand.xyz(2), fuzzball_ligand.xyz(3))

    ligand_rotation = target_fuzzball.R * mobile_match_inverse.R
    ligand_translation = target_fuzzball.R * mobile_match_inverse.t + target_fuzzball.t

    # Apply transformation
    match_pose_clone.apply_transform_Rx_plus_v(ligand_rotation, ligand_translation)

    # --- All other operations --- #

    # Mutate all non-motif residues within 10A from ligand to ALA, interferes with RotamerSet generation
    ligand_residue_selector = rosetta.core.select.residue_selector.ChainSelector('X')
    neighborhood_selector = rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand_residue_selector, 10, False)
    neighborhood_selector_bool = neighborhood_selector.apply(match_pose_clone)
    neighborhood_residues_resnums = rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
    positions_to_consider = list(set(neighborhood_residues_resnums) - set(motif_and_ligand_resnums))

    mutate = rosetta.protocols.simple_moves.MutateResidue()
    mutate.set_res_name('ALA')

    for position in positions_to_consider:
        if match_pose_clone.residue(position).name3() not in ['GLY', 'PRO'] and 'disulfide' not in match_pose_clone.residue(position).name():
            mutate.set_target(position)
            mutate.apply(match_pose_clone)

    # Build RotamerSets for each extrachi/sample level
    if dump_rotamerset_pdb:
        all_rotamersets = rosetta.core.pack.rotamer_set.RotamerSetsFactory.create_rotamer_sets(match_pose_clone)
        task_factory = rosetta.core.pack.task.TaskFactory()

        # NATRO positions TaskOp
        rotamer_candidates_rs = rosetta.core.select.residue_selector.ResidueIndexSelector(','.join([str(i) for i in match_residue_map.keys()]))
        natro_rs = rosetta.core.select.residue_selector.NotResidueSelector(rotamer_candidates_rs)
        natro_op = rosetta.core.pack.task.operation.OperateOnResidueSubset(
            rosetta.core.pack.task.operation.PreventRepackingRLT(), natro_rs)
        task_factory.push_back(natro_op)

        rotamersets_packer_task = task_factory.create_task_and_apply_taskoperations(match_pose_clone)

        all_rotamersets.set_task(rotamersets_packer_task)

    # Generate rotamers at each position
    for position in match_residue_map.keys():

        position_rotamer_list = list()
        rotamer_stats[position] = dict()

        if dump_rotamerset_pdb:
            current_rotamerset = rosetta.core.pack.rotamer_set.RotamerSetFactory.create_rotamer_set(match_pose_clone)

        # Keep rotamers that are compatible with minimal binding motif
        for contact_residue in match_residue_map[position]:

            # print(f'Considering position {position}: {contact_residue}')
            possible_contact_geometries = match_residue_map[position][contact_residue]

            # --- Prepare viable rotamers for each position --- #

            # Define packertask using neighborhood_selector
            packer_task = rosetta.core.pack.task.TaskFactory.create_packer_task(match_pose)
            packer_task.initialize_from_command_line()

            # Get boolean vector for packable positions and apply to packer task
            packable_positions = rosetta.utility.vector1_bool()
            packable_position_list = [True if i == position else False for i in range(1, match_pose.size())]
            for bool_value in packable_position_list:
                packable_positions.append(bool_value)
            packer_task.restrict_to_residues(packable_positions)

            # Only build rotamers for residues with Hbond donors/acceptors
            restrict_CAAs = rosetta.core.pack.task.operation.RestrictAbsentCanonicalAAS(position, rosetta.utility.vector1_bool(20))
            restrict_CAAs.keep_aas(contact_residue)
            restrict_CAAs.apply(match_pose, packer_task)

            packer_neighbor_graph = rosetta.core.pack.create_packer_graph(match_pose, sfxn, packer_task)

            match_rotamer_set = rosetta.core.pack.rotamer_set.RotamerSetFactory.create_rotamer_set(match_pose_clone)
            match_rotamer_set.set_resid(position)
            match_rotamer_set.build_rotamers(match_pose_clone, sfxn, packer_task, packer_neighbor_graph, use_neighbor_context=False)

            if match_rotamer_set.num_rotamers() <= 1 and match_rotamer_set.rotamer(1).name1() != contact_residue:
                continue

            # print(f'Comparing {match_rotamer_set.num_rotamers()} rotamers against {len(possible_contact_geometries)} contact modes')

            rotamer_stats[position][contact_residue] = dict()
            rotamer_stats[position][contact_residue]['num_rotamers'] = match_rotamer_set.num_rotamers()
            rotamer_info = list()
            rotamers_accepted = 0

            for rotamer in range(1, match_rotamer_set.num_rotamers() + 1):

                # Place residue before applying to pose!!!!
                # Rotamers need to be transformed back onto the backbone of the input pdb!!!
                trail_rotamer = match_rotamer_set.rotamer(rotamer)

                if dump_rotamerset_pdb:
                    trail_rotamer.place(match_pose.residue(position), match_pose.conformation_ptr())
                    current_rotamerset.add_rotamer(trail_rotamer)

                trail_rotamer.place(match_pose_clone.residue(position), match_pose_clone.conformation_ptr())
                match_pose_clone.replace_residue(position, trail_rotamer, False)

                # Evaluate RMSD to possible_contact_geometries
                contact_RMSDs = list()

                sad_atom_in_rotamer = False
                for contact_mode in possible_contact_geometries:

                    # Get contact atom coords using atom names
                    try:
                        rotamer_contact_coords = [list(match_pose_clone.residue(position).xyz(atom)) for atom in contact_mode[0]]
                    except Exception as e:
                        # print(f'Skipping {contact_mode[0]}: contains sad atom.')
                        sad_atom_in_rotamer = True
                        continue

                    contact_RMSDs.append(prody.calcRMSD(np.asarray(contact_mode[1]), np.asarray(rotamer_contact_coords)))

                if sad_atom_in_rotamer:
                    continue

                # Continue if current rotamer does not have <{RMSD_limit}A RMSD with any contact mode
                # todo: make RMSD satisfaction optional
                if min(contact_RMSDs) > RMSD_limit:
                    rotamer_info.append((contact_RMSDs, None, None))
                    continue

                # Evaluate possible clashes (fa_rep) with motif residues and ligand
                sfxn(match_pose_clone)
                edges = match_pose_clone.energies().energy_graph()

                motif_fa_rep = list()

                for motif in motif_and_ligand_resnums:
                    current_edge = edges.find_energy_edge(position, motif)
                    if current_edge is not None:
                        current_edge.fill_energy_map()
                        motif_fa_rep.append(current_edge[rosetta.core.scoring.fa_rep])

                # Get score for current rotamer against ligand
                current_edge = edges.find_energy_edge(position, conformer_resnum)
                rotamer_ligand_reu = current_edge.dot(sfxn_weights) if current_edge is not None else 0

                if all([len(motif_fa_rep) > 0, max(motif_fa_rep) < 10, rotamer_ligand_reu <= 5]):

                    new_rotamer = trail_rotamer.clone()
                    if flag_special_rot:
                        new_rotamer = rosetta.core.pose.add_variant_type_to_residue(new_rotamer, rosetta.core.chemical.SPECIAL_ROT, match_pose_clone)

                    # Place residue before applying to pose!!!!
                    # Rotamers need to be transformed back onto the backbone of the input pdb!!!
                    new_rotamer.place(match_pose.residue(position), match_pose.conformation_ptr())

                    position_rotamer_list.append(new_rotamer)
                    rotamers_accepted += 1

                rotamer_info.append((contact_RMSDs, max(motif_fa_rep), rotamer_ligand_reu))

            # print(f'{rotamers_accepted} of {match_rotamer_set.num_rotamers()} rotamers accepted')
            rotamer_stats[position][contact_residue]['rotamer_info'] = rotamer_info
            rotamer_stats[position][contact_residue]['rotamers_accepted'] = rotamers_accepted

        if len(position_rotamer_list) > 0:
            if position not in viable_rotamers.keys():
                viable_rotamers[position] = dict()
            viable_rotamers[position][contact_residue] = position_rotamer_list

        if dump_rotamerset_pdb:
            current_moltresid = all_rotamersets.resid_2_moltenres(position)
            all_rotamersets.set_explicit_rotamers(current_moltresid, current_rotamerset)

    if dump_rotamerset_pdb:
        current_extrachi = len([rosetta.basic.options.get_boolean_option(f'packing:ex{i}') for i in range(1,5) if rosetta.basic.options.get_boolean_option(f'packing:ex{i}') is True])
        current_sample_level = rosetta.basic.options.get_integer_option(f'packing:ex{current_extrachi}:level')

        if current_extrachi <= 2 and current_sample_level <= 3:
            match_name = os.path.normpath(os.path.basename(match_path))

            # todo: figure out why this doesn't work... problem with CONECT records...
            # all_rotamersets.dump_pdb(match_pose_clone, f"{match_name.split('.')[0]}-extrachi_{current_extrachi}-sampling_{current_sample_level}.pdb")

            all_rotamers_pose = pyrosetta.pose_from_sequence('A')

            for position in match_residue_map.keys():
                position_rotset = all_rotamersets.rotamer_set_for_residue(position)
                for rot in range(1, position_rotset.num_rotamers() + 1):
                    all_rotamers_pose.append_residue_by_jump(position_rotset.rotamer(rot), 1)
            all_rotamers_pose.dump_pdb(f"{match_name.split('.')[0]}-extrachi_{current_extrachi}-sampling_{current_sample_level}.pdb")

    if report_stats:
        return viable_rotamers, rotamer_stats
    else:
        return viable_rotamers


def fuzzball_composition_design(user_defined_dir, match_path, match_residue_map, design_config_json, nstruct=10):
    """
    Perform design using Vikram's AA_Composition score term, biasing toward rotamers that recapitulate contacts
    observed in the iteration fuzzball.
    :return:
    """

    # Get matched motif residue and ligand positions
    conformer_name, fuzzball_resnums = find_conformer_and_constraint_resnums(os.path.normpath(os.path.basename(match_path)))

    # --- Initiate PyRosetta and Score Function -- #

    my_options = [f"-extra_res_fa {os.path.join(user_defined_dir, 'Inputs', 'Rosetta_Inputs', f'{conformer_name}.params')}",
                  "-mute core.conformation core.chemical core.pack.task",
                  '-ex1 -ex2 -ex3 -ex4 -extrachi_cutoff 0',
                  '-packing:ex1:level 7 -packing:ex2:level 7 -packing:ex3:level 7 -packing:ex4:level 7']
    pyrosetta.init(options=' '.join(my_options))

    # Create match pose
    match_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(match_pose, match_path)

    # Normal scorefunction for generating rotamers
    sfxn = rosetta.core.scoring.get_score_function()

    # Add defined_rotamer scoreterm
    # sfxn.set_weight(rosetta.core.scoring.defined_rotamer, 1)

    # --- Set up Annealer for design --- #

    # Load design json
    design_json_info = json.load(open(design_config_json, 'r'))

    # Load viable scaffold positions and corresponding residue types
    # todo: make sure backrub ensemble structures also have matcher remarks added
    match_residue_map = pickle.load(open(match_residue_map, 'r'))

    # todo: write "layer selector" for two contact layer deep design shell

    # --- Residue Selectors --- #

    # Ligand
    matched_ligand_rs = rosetta.core.select.residue_selector.ChainSelector('X')

    # Matched motif residues
    matched_motif_residues = parse_matcher_remarks(match_pose)
    matched_motif_rs = rosetta.core.select.residue_selector.ResidueIndexSelector(','.join(matched_motif_residues))

    # User-defined design positions
    design_positions = [str(index) for index in design_json_info['design_residue_list']]
    design_position_rs = rosetta.core.select.residue_selector.ResidueIndexSelector(','.join(design_positions))

    # Packing shell around design/matched residues
    relevent_residue_rs = rosetta.core.select.residue_selector.OrResidueSelector()
    relevent_residue_rs.add_residue_selector(matched_motif_rs)
    relevent_residue_rs.add_residue_selector(design_position_rs)

    packing_shell_rs = rosetta.core.select.residue_selector.NeighborhoodResidueSelector(relevent_residue_rs, 8, include_focus_in_subset=True)

    # NATRO positions in pose
    natro_rs = rosetta.core.select.residue_selector.NotResidueSelector(packing_shell_rs)

    natro_or_designable_rs = rosetta.core.select.residue_selector.OrResidueSelector()
    natro_or_designable_rs.add_residue_selector(design_position_rs)
    natro_or_designable_rs.add_residue_selector(natro_rs)

    repack_only_rs = rosetta.core.select.residue_selector.NotResidueSelector(natro_or_designable_rs)

    # --- Create and Populate Task Factory --- #

    task_factory = rosetta.core.pack.task.TaskFactory()

    racaa = rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
    racaa.aas_to_keep('ADEFHIKLMNQRSTVWY') # No CGP
    design_op = rosetta.core.pack.task.operation.OperateOnResidueSubset(racaa, design_position_rs)
    task_factory.push_back(design_op)

    repack_op = rosetta.core.pack.task.operation.OperateOnResidueSubset(rosetta.core.pack.task.operation.RestrictToRepackingRLT(), repack_only_rs)
    task_factory.push_back(repack_op)

    natro_op = rosetta.core.pack.task.operation.OperateOnResidueSubset(rosetta.core.pack.task.operation.PreventRepackingRLT(), natro_rs)
    task_factory.push_back(natro_op)

    fixed_ligand_op = rosetta.core.pack.task.operation.OperateOnResidueSubset(rosetta.core.pack.task.operation.PreventRepackingRLT(), matched_ligand_rs)
    task_factory.push_back(fixed_ligand_op)

    # Extra rotamers
    extra_rotamers_op = rosetta.core.pack.task.operation.ExtraRotamersGeneric()
    extra_rotamers_op.ex1(True)
    extra_rotamers_op.ex2(True)
    extra_rotamers_op.ex1_sample_level(1)
    extra_rotamers_op.ex2_sample_level(1)
    task_factory.push_back(extra_rotamers_op)

    # --- Create RotamerSets including fuzzball rotamers --- #

    viable_rotamers = generate_fuzzball_contact_rotamersets(user_defined_dir, match_path, match_pose, sfxn, match_residue_map, flag_special_rot=True)

    # Turn off ex3 and ex4 after generating fuzzball contact rotamers
    rosetta.basic.options.set_boolean_option('packing:ex3', False)
    rosetta.basic.options.set_boolean_option('packing:ex4', False)

    # --- Create and apply filters --- #

    # todo: how the fuck do I add filters

    # --- Perform Design --- #

    design_packer_task = task_factory.create_task_and_apply_taskoperations(task_factory, match_pose)

    sfxn.setup_for_packing(match_pose, design_packer_task.repacking_residues(), design_packer_task.designing_residues())
    packer_neighbor_graph = rosetta.core.pack.create_packer_graph(match_pose, sfxn, design_packer_task)

    rotamer_sets = rosetta.core.pack.rotamer_set.RotamerSetsFactory.create_rotamer_sets(match_pose)
    rotamer_sets.set_task(design_packer_task)
    rotamer_sets.initialize_pose_for_rotsets_creation(match_pose)
    rotamer_sets.build_rotamers(match_pose, sfxn, packer_neighbor_graph)

    for position in viable_rotamers:
        position_rotamer_set = rotamer_sets.rotamer_set_for_residue(position)

        # Add fuzzball rotamers to the appropriate rotamer_set in rotamer_sets
        if int(position_rotamer_set.resid()) == position:

            for residue_type in viable_rotamers[position]:
                for fuzz_rotamer in viable_rotamers[position][residue_type]:
                    position_rotamer_set.add_rotamer_into_existing_group(fuzz_rotamer)

            position_rotamer_set.update_rotamer_offsets()

    sfxn.setup_for_packing_with_rotsets(match_pose, rotamer_sets)
    rotamer_sets.prepare_sets_for_packing(match_pose, sfxn)
    ig = rosetta.core.pack.interaction_graph.InteractionGraphFactory.create_and_initialize_annealing_graph(design_packer_task, rotamer_sets, match_pose, sfxn, packer_neighbor_graph)

    rosetta.core.pack.pack_rotamers.pack_rotamers_run(match_pose, design_packer_task, rotamer_sets, ig)
    ig.clean_up_after_packing(match_pose)
    sfxn(match_pose)

    # --- Write design to file --- #


# todo: add new energy method for applying bonus to flagged rotamer variants
# todo: figure out how to write a design guidance term using this implementation

from rosetta.core.scoring.methods import ContextIndependentOneBodyEnergy

@rosetta.EnergyMethod()
class FuzzballDesignMethod(ContextIndependentOneBodyEnergy):
    """
    Bias design toward contacts observed in the fuzzball
    """
    def __init__(self):
        ContextIndependentOneBodyEnergy.__init__(self, self.creator())

    def residue_energy(self, res, pose, emap):
        emap.get().set(self.scoreType, -1.0)