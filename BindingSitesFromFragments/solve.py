#!/usr/bin/env python3

import os
import sys
import copy
import math
import pickle
import random
import shutil
from pprint import pprint

import prody
import pandas as pd

from .utils import *

import pyrosetta
from pyrosetta import rosetta


def montecarlo_motif(user_defined_dir, fuzzball_path, motif_size, block_count=None, include_defined=False, solutions=100000):
    """
    Use PyRosetta and a simulated annealing monte carlo method to find binding site solutions.
    This method for finding binding motif solutions DOES NOT TAKE CONSTRAINTS INTO CONSIDERATION. This means motif
    residues may clash (residue-residue 2b REU > 0) and there is little control over binding motif composition i.e.
    hydrogen bonding residues will be seldomly included. Use Gurobi if this is necessary.

    :param user_defined_dir: Path to project root directory
    :param fuzzball_path: Path to directory containing fuzzballs for current iteration
    :param motif_size: Number of residues (including ligand) to include in binding motif solutions
    :param block_count: Number of trajectories per conformer
    """
    # --- Init variables --- #

    rosetta_inputs = os.path.join(user_defined_dir, 'Inputs', 'Rosetta_Inputs')
    conformer = os.path.basename(fuzzball_path).split('-')[0]  # 38E_0001

    # --- Init PyRosetta --- #

    my_options = [f"-extra_res_fa {os.path.join(rosetta_inputs, f'{conformer}.params')}",
                  "-mute core.conformation core.chemical"]
    pyrosetta.init(options=' '.join(my_options))

    # Custom score function using select 2b energies
    bsff_sfxn = rosetta.core.scoring.ScoreFunction()
    bsff_sfxn.set_weight(rosetta.core.scoring.fa_atr, 1)
    bsff_sfxn.set_weight(rosetta.core.scoring.fa_rep, 0.55)
    bsff_sfxn.set_weight(rosetta.core.scoring.fa_sol, 1)
    bsff_sfxn.set_weight(rosetta.core.scoring.fa_elec, 1)
    bsff_sfxn.set_weight(rosetta.core.scoring.hbond_sc, 1)

    # Load fuzzball into PyRosetta
    fuzzball_pose = rosetta.core.import_pose.pose_from_file(fuzzball_path)
    fuzzball_size = fuzzball_pose.size()  # Includes ligand
    print(f'\nFuzzball size: {fuzzball_size}\n')

    # Load fuzzball into prody and get defined residues
    fuzzball_base, fuzzball_pdb = os.path.split(fuzzball_path)
    fuzzball_agnpz = f'{fuzzball_pdb.split(".")[0]}.ag.npz'

    if include_defined:
        fuzzball_prody = prody.loadAtoms(os.path.join(fuzzball_base, fuzzball_agnpz))
        fuzzball_defined_residues = list(set(fuzzball_prody.select('defined 1').getResnums()))

    # --- Init motif pose --- #

    motif_pose = rosetta.core.pose.Pose()
    ligand = fuzzball_pose.residue(1).clone()
    motif_pose.append_residue_by_bond(ligand)

    if include_defined:
        for residue_index in fuzzball_defined_residues:
            defined_residue = fuzzball_pose.residue(residue_index)
            current_jump = motif_pose.size()
            motif_pose.append_residue_by_jump(defined_residue, current_jump)

        print(f'Initiated initial motif with defined residues {" ".join([str(a) for a in fuzzball_defined_residues])}')

    defined_motif_size = motif_pose.size()
    defined_motif_resnums = [i for i in range(1, defined_motif_size + 1)]

    minimum_base_score = bsff_sfxn(motif_pose)

    # Count hbonding residues in minimum motif
    defined_hbond_motifs = list()

    e_edge = motif_pose.energies().energy_graph()
    for residue in range(2, motif_pose.size() + 1):
        get_energy_edge = e_edge.find_energy_edge(1, residue)
        current_edge = get_energy_edge.fill_energy_map()
        hbond_sc = current_edge.get(rosetta.core.scoring.hbond_sc)
        if hbond_sc is not None and hbond_sc < 0:
            defined_hbond_motifs.append(residue)

    motif_residues_to_init = motif_size - defined_motif_size
    current_residues_in_fuzzball = dict()  # {motif_position: fuzzball_motif_index}


    def get_new_residue(current_residues_in_fuzzball):
        residue_from_fuzzball = int(random.randint(2, fuzzball_size))
        if residue_from_fuzzball not in current_residues_in_fuzzball.values():
            return residue_from_fuzzball
        else:
            return get_new_residue(current_residues_in_fuzzball)

    for i in range(motif_residues_to_init):
        new_residue_index = get_new_residue(current_residues_in_fuzzball)
        random_motif_from_fuzzball = fuzzball_pose.residue(new_residue_index).clone()

        current_jump = motif_pose.size()
        motif_pose.append_residue_by_jump(random_motif_from_fuzzball, current_jump)
        current_residues_in_fuzzball[current_jump + 1] = new_residue_index

    print(f'Initiated initial motif with random residues {" ".join([str(a) for a in current_residues_in_fuzzball])}')

    bsff_sfxn(motif_pose)

    # --- Get hydrogen-bonding residues in fuzzball --- #

    hbond_donor_list = list()
    bsff_sfxn(fuzzball_pose)
    e_edge = fuzzball_pose.energies().energy_graph()

    for residue in range(2 + defined_motif_size, fuzzball_pose.size() + 1):
        get_energy_edge = e_edge.find_energy_edge(1, residue)

        if get_energy_edge is None:
            continue

        current_edge = get_energy_edge.fill_energy_map()
        hbond_sc = current_edge.get(rosetta.core.scoring.hbond_sc)
        if hbond_sc is not None and hbond_sc < 0:
            hbond_donor_list.append(residue)

    hbond_residue_set = set(hbond_donor_list)

    # --- Perform Simulated Annealing MC --- #

    # Set up infrastructure
    solution_list = list()
    solution_set = set()
    temperatures = [50, 10, 5, 1.5, 0.9, 0.6, 0.3]

    for kT in temperatures:
        print(f'Temperature: {kT}')
        mc = rosetta.protocols.moves.MonteCarlo(motif_pose, bsff_sfxn, kT)

        for trial in range(2000 * fuzzball_size):

            trial_pose = motif_pose.clone()

            # Get random key:value pair and replace residue in motif
            current_motif_position, current_fuzzball_motif = random.choice(list(current_residues_in_fuzzball.items()))
            new_residue_index = get_new_residue(current_residues_in_fuzzball)
            random_motif_from_fuzzball = fuzzball_pose.residue(new_residue_index).clone()

            trial_pose.replace_residue(current_motif_position, random_motif_from_fuzzball, False)

            # Only store new and improved solutions
            accepted = mc.boltzmann(trial_pose)
            if accepted and mc.last_accepted_score() < minimum_base_score:

                # todo: Apply additional constraints here
                # Add HBond donors/acceptors
                # if not len(hbond_residue_set & set(defined_motif_resnums + current_residues_in_fuzzball.keys())) - len(defined_hbond_motifs) >= 1:
                #     continue
                
                # Record trial                
                motif_pose = trial_pose.clone()

                current_residues_in_fuzzball[current_motif_position] = new_residue_index
                new_solution_indicies = tuple(sorted(current_residues_in_fuzzball.values()))
                if new_solution_indicies not in solution_set:
                    solution_set.add(new_solution_indicies)
                    residue_indices_cell = [1] + list(new_solution_indicies)
                    if include_defined:
                        residue_indices_cell += fuzzball_defined_residues
                    solution_list.append({'Residue_indicies': residue_indices_cell,
                                          'Obj_score': mc.last_accepted_score(),
                                          })

            # Update minimum base score 
            if trial % solutions == 0:
                sorted_scores = sorted([a['Obj_score'] for a in solution_list])

                # todo: clean this up later...
                if len(sorted_scores) > solutions:
                    minimum_base_score = sorted_scores[solutions]
                elif len(sorted_scores) > 0:
                    minimum_base_score = sorted_scores[-1]

                print(f'New minimum base score: {minimum_base_score}')

    df = pd.DataFrame(solution_list)
    df_head = df.sort_values(by=['Obj_score'], ascending=True)[:solutions]
    df_head['Conformer'] = fuzzball_pdb.split('.')[0]

    if block_count is None:
        df_csv_name = f'{fuzzball_pdb.split(".")[0]}-{motif_size - 1}_residue_solutions.csv'
    else:
        df_csv_name = f'{fuzzball_pdb.split(".")[0]}-{motif_size - 1}_residue_solutions-{block_count}.csv'

    df_head.to_csv(df_csv_name, index=False)


def mc_crystallize(user_defined_dir, match_path, fuzzball_dir, fuzzball_index, positions_to_init=2, block_count=None):
    """
    Use PyRosetta and a simulated annealing monte carlo method to find binding site solutions in the context of a match.
    :return:
    """

    # Get matched motif residue and ligand positions
    conformer_name, fuzzball_resnums = find_conformer_and_constraint_resnums(os.path.normpath(os.path.basename(match_path)))
    motif_resname_resnum = determine_matched_residue_positions(match_path)

    # --- Initiate PyRosetta and Score Function -- #

    my_options = [f"-extra_res_fa {os.path.join(user_defined_dir, 'Inputs', 'Rosetta_Inputs', f'{conformer_name}.params')}",
                  "-mute core.conformation core.chemical",
                  '-ex1 -ex2 -ex3 -extrachi_cutoff 0 -ex1_sample_level 7 -ex2_sample_level 7 -ex3_sample_level 7']
    pyrosetta.init(options=' '.join(my_options))

    # Create match pose
    match_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(match_pose, match_path)

    # Custom score function using select 2b energies, consider hbond_bb_sc now that binding site is nucleated
    bsff_sfxn = rosetta.core.scoring.ScoreFunction()
    bsff_sfxn.set_weight(rosetta.core.scoring.fa_atr, 1)
    bsff_sfxn.set_weight(rosetta.core.scoring.fa_rep, 0.55)
    bsff_sfxn.set_weight(rosetta.core.scoring.fa_sol, 1)
    bsff_sfxn.set_weight(rosetta.core.scoring.fa_elec, 1)
    bsff_sfxn.set_weight(rosetta.core.scoring.hbond_sc, 1)
    bsff_sfxn.set_weight(rosetta.core.scoring.hbond_bb_sc, 1)
    bsff_sfxn(match_pose)

    sfxn_weights = bsff_sfxn.weights()

    # Normal scorefunction for generating rotamers
    sfxn = rosetta.core.scoring.get_score_function()

    # --- Load viable scaffold positions and corresponding residue types --- #

    match_residue_map_temp = pickle.load(open(os.path.join(fuzzball_dir, f'fuzz_{fuzzball_index}-match_residue_map.pickle'), 'r'))

    # Only consider non-CGP positions
    match_residue_map = dict()
    for position, residues in match_residue_map_temp.items():
        if match_pose.residue(int(position)).name1() not in ['C', 'G', 'P']:
            match_residue_map[int(position)] = residues

    conformer_resnum = match_pose.size()  # Assumes single ligand appended to end of sequence
    motif_resnums = [resnum for resname, resnum in motif_resname_resnum]
    motif_and_ligand_resnums = motif_resnums + [conformer_resnum]

    # DEBUGGING
    print(motif_and_ligand_resnums)
    pprint(match_residue_map)

    # --- Prepare viable rotamers for each position --- #

    # Define packertask using neighborhood_selector
    packer_task = rosetta.core.pack.task.TaskFactory.create_packer_task(match_pose)

    # Get boolean vector for packable positions and apply to packer task
    packable_positions = rosetta.utility.vector1_bool()
    packable_position_list = [True if i in match_residue_map.keys() else False for i in range(1, match_pose.size())]
    for bool_value in packable_position_list:
        packable_positions.append(bool_value)

    # DEBUGGING
    print(packable_positions)

    packer_task.restrict_to_residues(packable_positions)

    # Only build rotamers for residues with Hbond donors/acceptors
    for position, residue_types in match_residue_map.items():
        restrict_CAAs = rosetta.core.pack.task.operation.RestrictAbsentCanonicalAAS(position, rosetta.utility.vector1_bool(20))
        restrict_CAAs.keep_aas(''.join(residue_types))
        restrict_CAAs.apply(match_pose, packer_task)

    packer_neighbor_graph = rosetta.core.pack.create_packer_graph(match_pose, sfxn, packer_task)

    # --- Find and store viable rotamers --- #

    viable_rotamers = dict()

    # Setting things up is going to mess up the match pose, so use a clone
    match_pose_clone = match_pose.clone()

    # Mutate all non-motif residues within 8A from ligand to ALA, interferes with RotamerSet generation
    ligand_residue_selector = rosetta.core.select.residue_selector.ChainSelector('X')
    neighborhood_selector = rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand_residue_selector, 8, False)
    neighborhood_selector_bool = neighborhood_selector.apply(match_pose_clone)
    neighborhood_residues_resnums = rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
    positions_to_consider = list(set(neighborhood_residues_resnums) - set(motif_and_ligand_resnums))

    mutate = rosetta.protocols.simple_moves.MutateResidue()
    mutate.set_res_name('ALA')

    for position in positions_to_consider:
        if match_pose_clone.residue(position).name3() not in ['GLY', 'PRO'] and 'disulfide' not in match_pose_clone.residue(position).name():
            mutate.set_target(position)
            mutate.apply(match_pose_clone)

    # Generate rotamers at each position
    for position in match_residue_map.keys():

        match_rotamer_set = rosetta.core.pack.rotamer_set.RotamerSetFactory.create_rotamer_set(match_pose_clone)
        match_rotamer_set.set_resid(position)
        match_rotamer_set.build_rotamers(match_pose_clone, sfxn, packer_task, packer_neighbor_graph, use_neighbor_context=False)

        position_rotamer_list = list()

        # Keep rotamers that are compatible with minimal binding motif
        print(position, match_rotamer_set.num_rotamers())

        for rotamer in range(1, match_rotamer_set.num_rotamers() + 1):

            match_pose_clone.replace_residue(position, match_rotamer_set.rotamer(rotamer), False)

            # Evaluate possible clashes (fa_rep) with motif residues and ligand
            bsff_sfxn(match_pose_clone)
            edges = match_pose_clone.energies().energy_graph()

            motif_fa_rep = list()

            for motif in motif_and_ligand_resnums:
                current_edge = edges.find_energy_edge(position, motif)
                if current_edge is not None:
                    current_edge.fill_energy_map()
                    motif_fa_rep.append(current_edge[rosetta.core.scoring.fa_rep])

            if len(motif_fa_rep) > 0 and max(motif_fa_rep) < 10:

                # Get score for current rotamer against ligand
                current_edge = edges.find_energy_edge(position, conformer_resnum)
                rotamer_ligand_reu = current_edge.dot(sfxn_weights) if current_edge is not None else 0

                if rotamer_ligand_reu <= 0:
                    position_rotamer_list.append(match_rotamer_set.rotamer(rotamer))

        if len(position_rotamer_list) > 0:
            viable_rotamers[position] = position_rotamer_list

    # --- Perform packing with viable rotamers --- #
    """Writing our own Monte Carlo implementation since we score different parts of the pose during each trial"""

    def get_new_position(current_trial):
        match_position = np.random.choice(list(viable_rotamers.keys()))
        if match_position not in current_trial.values():
            return match_position
        else:
            return get_new_position(current_trial)

    def boltzmann(pose, positions, temp):
        """Accept/reject pose scoring current positions"""

        selector = rosetta.core.select.residue_selector.ResidueIndexSelector()
        for res in positions:
            selector.append_index(res)

        current_motif_score = bsff_sfxn.get_sub_score(pose, selector.apply(pose))
        delta_e = current_motif_score - previous_accepted
        boltzmann_factor = min(40, max(-40, -delta_e / temp))  # Copying Rosetta, does this to prevent overflow error
        probability = math.exp(boltzmann_factor)

        accepted = True

        if probability < 1:
            accepted = False if np.random.uniform() >= probability else True

        return accepted, current_motif_score

    # Set up infrastructure
    solution_list = list()
    solution_set = set()

    temperatures = [1.8, 1.5, 1.2, 0.9, 0.6, 0.3]
    previous_accepted = 9999

    residue_index_selector = rosetta.core.select.residue_selector.ResidueIndexSelector()
    for res in motif_and_ligand_resnums:
        residue_index_selector.append_index(res)

    base_match_score = bsff_sfxn.get_sub_score(match_pose, residue_index_selector.apply(match_pose))
    print(f'Base binding motif score:\t{base_match_score}')

    total_rotamers = sum([len(rotamers) for position, rotamers in viable_rotamers.items()])

    # Initialize initial motif with +2 random positions/rotamers
    # Keep track of {position: rotamer_index} for expanded motifs
    current_trial = dict()

    for i in range(positions_to_init):
        new_position = get_new_position(current_trial)
        new_rotamer = int(np.random.random_integers(low=1, high=len(viable_rotamers[new_position]))) - 1
        current_trial[new_position] = new_rotamer

    # Apply initial rotamers at selected positions to pose
    for position, rotamer in current_trial.items():
        match_pose.replace_residue(position, viable_rotamers[position][rotamer], False)

    last_accepted = copy.deepcopy(current_trial)

    # Perform trials
    for kT in temperatures:
        print(f'Temperature: {kT}')

        for trial in range(2500 * total_rotamers):

            # --- Add/remove/replace --- #

            action = np.random.choice(['add', 'remove', 'replace'], p=[0.1, 0.05, 0.85])

            if action == 'replace':
                # Remove an existing position/rotamer
                current_position, current_rotamer = random.choice(list(current_trial.items()))
                del current_trial[current_position]

                # Add new rotamer at random position
                new_position = random.choice(list(viable_rotamers.keys()))
                new_rotamer = np.random.random_integers(low=1, high=len(viable_rotamers[new_position])) - 1
                current_trial[new_position] = new_rotamer

                match_pose.replace_residue(new_position, viable_rotamers[new_position][new_rotamer], False)

            elif action == 'add':
                # if len(current_trial) > 10:
                #     current_position, current_rotamer = random.choice(list(current_trial.items()))
                #     match_pose.replace_residue(current_position, viable_rotamers[current_position][current_rotamer], False)
                # else:

                new_position = get_new_position(current_trial)
                new_rotamer = np.random.random_integers(low=1, high=len(viable_rotamers[new_position])) - 1
                current_trial[new_position] = new_rotamer
                match_pose.replace_residue(new_position, viable_rotamers[new_position][new_rotamer], False)

            else:
                current_position, current_rotamer = random.choice(list(current_trial.items()))
                if len(current_trial) > 1:
                    del current_trial[current_position]

                # Add a new residue if current_trial gets too small
                else:
                    new_position = get_new_position(current_trial)
                    new_rotamer = np.random.random_integers(low=1, high=len(viable_rotamers[new_position])) - 1
                    current_trial[new_position] = new_rotamer
                    match_pose.replace_residue(new_position, viable_rotamers[new_position][new_rotamer], False)

            # --- Evaluate --- #

            new_positions = sorted(list(current_trial.keys()))
            positions_to_score = motif_and_ligand_resnums + new_positions
            current_accepted, current_score = boltzmann(match_pose, positions_to_score, kT)

            if current_accepted:

                # print(f'Accepted:\t{current_score}\t{positions_to_score}')
                # todo: enforce hbond constraints

                new_solution = tuple([(position, rotamer) for position, rotamer in sorted(current_trial.items(), key=lambda x:x[1])])

                if new_solution not in solution_set:
                    solution_set.add(new_solution)
                    solution_list.append({'rotamers': new_solution,
                                          'Obj_score': current_score,
                                          })
                last_accepted = copy.deepcopy(current_trial)
                previous_accepted = current_score

            else:
                current_trial = copy.deepcopy(last_accepted)

    # DEBUGGING
    pprint(solution_list)

    write_solutions_to_pdb(match_path, solution_list, viable_rotamers, use_cluster_scratch=True, block_count=block_count, max_score=base_match_score)

    return solution_list


def write_solutions_to_pdb(match_path, solution_list, viable_rotamers, output_dir=os.getcwd(), use_cluster_scratch=True, block_count=0, max_score=0, max_solutions=10000):
    """
    Write binding site solutions to PDB
    :return:
    """

    def get_scaffold_from_match(match):
        """ Get scaffold name from match"""
        match_split = match.split('_')
        right_slice_index = match_split.index(ligand_resname)
        return '_'.join(match_split[4:right_slice_index])

    # Get conformer name and nucleated match positions
    conformer_name, motif_resnums = find_conformer_and_constraint_resnums(os.path.normpath(os.path.basename(match_path)))
    ligand_resname = conformer_name.split('_')[0]

    match_dir, match_name = os.path.split(match_path)
    match_scaffold = get_scaffold_from_match(match_name)
    match_nucleated_fuzz_indicies = match_name.split('-')[1]
    fuzzball_index = int(re.split('[_-]', match_name)[-2])

    pdb_output_dir = f'{match_name.split(".")[0]}-{block_count}'
    write_pdb_dir = os.path.join(os.environ.get("TMPDIR"), pdb_output_dir) if use_cluster_scratch else os.path.join(output_dir, pdb_output_dir)
    os.makedirs(write_pdb_dir, exist_ok=True)

    # Create match pose
    match_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(match_pose, match_path)

    # Iterate through solutions and output PDB
    curated_solution_list = sorted([solution for solution in solution_list if solution['Obj_score'] < max_score], key=lambda x: x['Obj_score'])[:max_solutions]

    for solution in curated_solution_list:

        if solution['Obj_score'] < max_score:
            match_pose_output = match_pose.clone()

            # Add solutions to match
            for position, rotamer in solution['rotamers']:
                match_pose_output.replace_residue(position, viable_rotamers[position][rotamer], False)

            # Write remarks to pose in style of matches
            all_remarks = rosetta.core.io.Remarks()
            match_positions = motif_resnums + [position for position, rotamer in solution['rotamers']]

            for index, match_postion in enumerate(match_positions, start=1):
                current_solution_remarks = rosetta.core.io.RemarkInfo()
                current_solution_remarks.num = 666
                current_solution_remarks.value = rosetta.protocols.toolbox.match_enzdes_util.assemble_remark_line('X', ligand_resname, 0, match_pose.pdb_info().chain(match_postion), match_pose.residue(match_postion).name3(), match_postion, index, 1)
                all_remarks.append(current_solution_remarks)

            match_pose_output.pdb_info().remarks(all_remarks)

            # Dump pose
            match_resname_resnum = ''.join([f'{match_pose_output.residue(position).name1()}{position}' for position in match_positions])
            match_name = f'UM_1_{match_resname_resnum}_{match_scaffold}_{conformer_name}-{match_nucleated_fuzz_indicies}-fuzz_{fuzzball_index}-xtal.pdb'
            rosetta.core.io.dump_pdb(match_pose_output, os.path.join(write_pdb_dir, match_name))

    if use_cluster_scratch:
        shutil.move(write_pdb_dir, output_dir)
