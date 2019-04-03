#!/usr/bin/env python3

import os
import sys
import random
import pprint

import prody
import pandas as pd
from .utils import *

import pyrosetta
from pyrosetta import rosetta

def montecarlo_motif(user_defined_dir, fuzzball_path, motif_size, block_count=None):
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

    # Load fuzzball into prody and get defined residues
    fuzzball_base, fuzzball_pdb = os.path.split(fuzzball_path)
    fuzzball_agnpz = f'{fuzzball_pdb.split(".")[0]}.ag.npz'
    fuzzball_prody = prody.loadAtoms(os.path.join(fuzzball_base, fuzzball_agnpz))
    fuzzball_defined_residues = list(set(fuzzball_prody.select('defined 1').getResnums()))

    # --- Init motif pose --- #

    motif_pose = rosetta.core.pose.Pose()
    ligand = fuzzball_pose.residue(1).clone()
    motif_pose.append_residue_by_bond(ligand)
    for residue_index in fuzzball_defined_residues:
        defined_residue = fuzzball_pose.residue(residue_index)
        current_jump = motif_pose.size()
        motif_pose.append_residue_by_jump(defined_residue, current_jump)

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

    print(f'Initiated initial motif with defined residues {" ".join([str(a) for a in fuzzball_defined_residues])}')

    def get_new_residue(current_residues_in_fuzzball):
        residue_from_fuzzball = int(random.randint(motif_pose.size() + 1, fuzzball_size))
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

    for residue in range(defined_motif_size, fuzzball_pose.size() + 1):
        get_energy_edge = e_edge.find_energy_edge(1, residue)
        current_edge = get_energy_edge.fill_energy_map()
        hbond_sc = current_edge.get(rosetta.core.scoring.hbond_sc)
        if hbond_sc is not None and hbond_sc < 0:
            hbond_donor_list.append(residue)

    hbond_residue_set = set(hbond_donor_list)

    # --- Perform Simulated Annealing MC --- #

    # Set up infrastructure
    solution_list = list()
    solution_set = set()
    temperatures = [1.5, 1.2, 0.9, 0.6, 0.3]

    for kT in temperatures:
        print(f'Temperature: {kT}')
        mc = rosetta.protocols.moves.MonteCarlo(motif_pose, bsff_sfxn, kT)

        for trial in range(2500 * fuzzball_size):

            # Get random key:value pair and replace residue in motif
            current_motif_position, current_fuzzball_motif = random.choice(list(current_residues_in_fuzzball.items()))
            new_residue_index = get_new_residue(current_residues_in_fuzzball)
            random_motif_from_fuzzball = fuzzball_pose.residue(new_residue_index).clone()

            motif_pose.replace_residue(current_motif_position, random_motif_from_fuzzball, False)

            if mc.boltzmann(motif_pose) and mc.last_accepted_score() < minimum_base_score:

                # todo: Apply additional constraints here
                # Add HBond donors/acceptors
                # if not len(hbond_residue_set & set(defined_motif_resnums + current_residues_in_fuzzball.keys())) - len(defined_hbond_motifs) >= 1:
                #     continue

                current_residues_in_fuzzball[current_motif_position] = new_residue_index
                new_solution_indicies = tuple(sorted(current_residues_in_fuzzball.values()))
                if new_solution_indicies not in solution_set:
                    solution_set.add(new_solution_indicies)
                    solution_list.append({'Residue_indicies': [1] + list(new_solution_indicies) + fuzzball_defined_residues,
                                          'Obj_score': mc.last_accepted_score(),
                                          })

    df = pd.DataFrame(solution_list)
    df['Conformer'] = fuzzball_pdb.split('.')[0]

    if block_count is None:
        df_csv_name = f'{fuzzball_pdb.split(".")[0]}-{motif_size - 1}_residue_solutions.csv'
    else:
        df_csv_name = f'{fuzzball_pdb.split(".")[0]}-{motif_size - 1}_residue_solutions-{block_count}.csv'

    df.to_csv(df_csv_name, index=False)

