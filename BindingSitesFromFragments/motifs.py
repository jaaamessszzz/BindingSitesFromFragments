#!/usr/bin/env python3

import io
import re
import os
import copy
import json
import collections
from pprint import pprint
from ast import literal_eval

import yaml
import prody
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import rdFMCS


from .utils import *

# --- Silence ProDy --- #
prody.confProDy(verbosity='none')

class Generate_Constraints(object):

    three_to_one = {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K',
                    'LEU':'L', 'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V',
                    'TRP':'W', 'TYR':'Y'}

    def __init__(self, user_defined_dir):
        self.user_defined_dir = user_defined_dir
        self.single_pose_dir = os.path.join(self.user_defined_dir, 'Motifs', 'Residue_Ligand_Interactions', 'Single_Poses')

        # DEPRECIATED
        self.fragment_cluster_list = yaml.load(open(os.path.join(self.user_defined_dir, 'Inputs', 'User_Inputs', 'Motif_Clusters.yml'), 'r'))
        self.res_idx_map_df = False
        self.reference_pose_dir = os.path.join(self.user_defined_dir, 'Motifs', 'Reference_Pose')

    # DEPRECIATED
    def import_res_idx_map(self):
        self.res_idx_map_df = pd.read_csv(os.path.join(self.single_pose_dir, 'residue_index_mapping.csv'),
                                          usecols=['residue_index', 'source_conformer', 'source_pdb', 'struct_id']
                                          )

    def determine_constraint_atoms_from_single_pose(self, single_pose_prody, current_conformer, current_residue_index, verbose=False):
        """
        Use single pose to determine the three atoms on the ligand and residue to use for matcher constraints

        :param single_pose_prody: prody object of single pose for current conformer
        :param current_conformer: name of current conformer e.g. MEH_0001
        :param current_residue_index: residue number for current residue in single pose
        :return:
        """
        residue_prody = single_pose_prody.select('resnum {} and not hydrogen'.format(current_residue_index)).copy()
        residue_index_row = self.res_idx_map_df.loc[(self.res_idx_map_df['residue_index'] == current_residue_index) &
                                                    (self.res_idx_map_df['source_conformer'] == current_conformer)]
        residue_split = re.split('-|\.', residue_index_row['source_pdb'].values[0])

        # For representative residues or cluster-sourced
        # EDIT: if any of the MySQL crap doesn't work anymore, this is why.
        source_pdb = residue_index_row['source_pdb'].values[0]
        residue_index = int(residue_split[2].split('_')[1])
        fragment = residue_split[0]
        cluster = residue_split[1]
        resname = residue_prody.getResnames()[0]

        constraints_dict = self.determine_constraint_atoms(residue_prody, current_conformer, fragment=fragment)

        constraints_dict['source_pdb'] = source_pdb
        constraints_dict['residue_index'] = residue_index
        constraints_dict['fragment'] = fragment
        constraints_dict['cluster'] = cluster
        constraints_dict['resname'] = resname

        return constraints_dict

    def determine_constraint_atoms(self, residue_prody, conformer_prody, fragment=None, verbose=False):
        """
        Determine the three atoms in the ligand and the residue that will be used to calculate ideal values for the DOFs
        required for matcher.

        :param residue_prody: prody object of residue
        :param conformer_prody: prody object of current conformer
        :param fragment: string of current fragment (e.g. "Fragment_1")
        :return: 
        """

        ligand_code = conformer_prody.getResnames()[0]

        residue_index_atom_map = {atom.getIndex(): atom.getName() for atom in residue_prody.select('not hydrogen')}
        residue_atom_index_map = {v: k for k, v in residue_index_atom_map.items()}
        ligand_index_atom_map = {atom.getIndex(): atom.getName() for atom in conformer_prody.select('not hydrogen')}

        # Load residue and ligand as RDKit Mol objects
        RD_residue = RDKit_Mol_from_ProDy(residue_prody, removeHs=True)
        RD_ligand = RDKit_Mol_from_ProDy(conformer_prody, removeHs=False)

        # Calculate closest atom-atom contacts and two additional atoms for determining bond torsions and angles
        # NOTE: Contact distance and indicies are for residue and ligand with hydrogens stripped!

        # Select source fragment atoms from current conformer is fragment is provided, else use entire ligand
        if fragment:
            fragment_atoms_prody = prody.parsePDB(os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs', f'{fragment}.pdb'))
            fragment_atom_names = fragment_atoms_prody.select('not hydrogen').getNames()
            ligand_atoms_for_determining_contacts = conformer_prody.select(f'resname {ligand_code} and name {" ".join(fragment_atom_names)}')

        else:
            ligand_atoms_for_determining_contacts = conformer_prody

        contact_distance, residue_index_low, ligand_index_low = minimum_contact_distance(residue_prody, ligand_atoms_for_determining_contacts, return_indices=True)

        residue_contact_atom = residue_prody.select('not hydrogen')[residue_index_low]
        ligand_contact_atom = ligand_atoms_for_determining_contacts.select('not hydrogen')[ligand_index_low]

        # Select two additional atoms downstream from contact atom for meaningful constraints... programmatically...
        # Using RDKit to determine neighboring atoms, removes Hydrogens by default

        # RESIDUE CONSTRAINT ATOM INDICES
        # Special cases: O (O>C>CA), N (N>CA>C), C (C>CA>N), CA (CA>C>O)
        if residue_contact_atom.getName() in ['C', 'CA', 'CB', 'N', 'O', 'OXT']:
            if residue_contact_atom.getName() == 'CB':
                residue_second_atom = residue_atom_index_map['CA']
                residue_third_atom = residue_atom_index_map['C']
            elif residue_contact_atom.getName() == 'O':
                residue_second_atom = residue_atom_index_map['C']
                residue_third_atom = residue_atom_index_map['CA']
            elif residue_contact_atom.getName() == 'N':
                residue_second_atom = residue_atom_index_map['CA']
                residue_third_atom = residue_atom_index_map['C']
            elif residue_contact_atom.getName() == 'C':
                residue_second_atom = residue_atom_index_map['CA']
                residue_third_atom = residue_atom_index_map['N']
            elif residue_contact_atom.getName() == 'CA':
                residue_second_atom = residue_atom_index_map['C']
                residue_third_atom = residue_atom_index_map['O']
            elif residue_contact_atom.getName() == 'OXT':
                residue_second_atom = residue_atom_index_map['C']
                residue_third_atom = residue_atom_index_map['CA']

            else:
                raise Exception('wut.')

        else:
            residue_second_atom = self._determine_next_residue_constraint_atom(residue_contact_atom.getIndex(), RD_residue, residue_prody)
            residue_third_atom = self._determine_next_residue_constraint_atom(residue_second_atom, RD_residue, residue_prody)

        ligand_second_atom, ligand_third_atom = self._determine_ligand_constraint_atoms(ligand_contact_atom.getIndex(), RD_ligand, ligand_atoms_for_determining_contacts)

        if verbose:

            print('RESIDUE - FIRST ATOM: {} {}'.format(residue_contact_atom.getIndex(), residue_index_atom_map[residue_contact_atom.getIndex()]))
            print('RESIDUE - SECOND ATOM: {} {}'.format(residue_second_atom, residue_index_atom_map[residue_second_atom]))
            print('RESIDUE - THIRD ATOM: {} {}'.format(residue_third_atom, residue_index_atom_map[residue_third_atom]))

            # LIGAND CONSTRAINT ATOM INDICES
            # So as far as choosing ligand constraints go, I think it's safe to chose whatever atoms as long as the
            # terminal atom has >1 neighbor. The rationale being this will propagate atom selection toward less
            # flexible parts of the ligand... hopefully...

            print('LIGAND - FIRST ATOM: {} {}'.format(ligand_contact_atom.getIndex(), ligand_index_atom_map[ligand_contact_atom.getIndex()]))
            print('LIGAND - SECOND ATOM: {} {}'.format(ligand_second_atom, ligand_index_atom_map[ligand_second_atom]))
            print('LIGAND - THIRD ATOM: {} {}'.format(ligand_third_atom, ligand_index_atom_map[ligand_third_atom]))

        return {'contact_distance': float(contact_distance),
                'ligand': {'atom_names': [ligand_index_atom_map[ligand_contact_atom.getIndex()],
                                          ligand_index_atom_map[ligand_second_atom],
                                          ligand_index_atom_map[ligand_third_atom]
                                          ],
                           'atom_indices': [ligand_contact_atom.getIndex(),
                                            ligand_second_atom,
                                            ligand_third_atom]
                           },
                'residue': {'atom_names': [residue_index_atom_map[residue_contact_atom.getIndex()],
                                           residue_index_atom_map[residue_second_atom],
                                           residue_index_atom_map[residue_third_atom]
                                           ],
                            'atom_indices': [residue_contact_atom.getIndex(),
                                             residue_second_atom,
                                             residue_third_atom]
                            }
                }

    def _determine_next_residue_constraint_atom(self, current_atom_index, RD_residue, residue_prody):
        """
        Grabs the next residue atom down the line for residue constraint records
        :return: 
        """
        # Okay so I can't figure out for now how to retrieve atom names from RDKit atom objects,,,
        # I'm going to create a dict mapping atom indices to atom names from the prody object
        # Prody and RDKit both seem to derive atom indices in order of the PDB file
        residue_index_atom_map = {atom.getIndex(): atom.getName() for atom in residue_prody.select('not hydrogen')}
        residue_atom_index_map = {v: k for k, v in residue_index_atom_map.items()}

        # Setting up hierarchy for residue atom selections - want selections to propagate toward main chain
        # So.... H>Z>E>D>G>B>A
        # Okay so starting at any given atom I should be able to determine connectivity and how many atoms out I am from CA
        atom_hierarchy = ['H', 'Z', 'E', 'D', 'G', 'B', 'A']

        residue_contact_atom_neighbors = [atom.GetIdx() for atom in
                                          RD_residue.GetAtomWithIdx(int(current_atom_index)).GetNeighbors() if
                                          atom.GetSymbol != 'H']

        if len(residue_contact_atom_neighbors) == 1:
            next_atom_index = residue_contact_atom_neighbors[0]
            return next_atom_index


        else:
            first_atom_ID = residue_index_atom_map[current_atom_index]
            # Special case for Prolines. CG > N causes issues with atom[1]... is there another fix? Yes. But it's 12am
            second_atom_IDs = [residue_index_atom_map[atom] for atom in residue_contact_atom_neighbors if len(residue_index_atom_map[atom]) > 1]

            for atom in second_atom_IDs:
                if atom[1] == atom_hierarchy[atom_hierarchy.index(first_atom_ID[1]) + 1]:
                    next_atom_index = residue_atom_index_map[atom]
                    return next_atom_index

            raise Exception('Somehow you\'ve run out of neighbors for your contact atom...')

    def _determine_ligand_constraint_atoms(self, current_atom_index, RD_ligand, ligand_fragment_prody):
        """
        Grabs the next ligand atoms down the line for ligand constraint records
        :param current_atom_index: index of ligand contact atom in RD_ligand
        :param RD_ligand: RDKit molecule of FULL ligand
        :param ligand_fragment_prody: prody of ligand FRAGMENT from current conformer
        :return: 
        """
        ligand_index_atom_map = {atom.getIndex(): atom.getName() for atom in ligand_fragment_prody.select('not hydrogen')}
        ligand_contact_atom_neighbors = [atom for atom in
                                         RD_ligand.GetAtomWithIdx(int(current_atom_index)).GetNeighbors() if atom.GetIdx() in ligand_index_atom_map.keys()]

        # Check all atoms that are adjacent to ligand contact atom
        for second_atom in ligand_contact_atom_neighbors:

            # Only investigate potential second constraint atoms with more than one neighbor
            if len([atom.GetIdx() for atom in second_atom.GetNeighbors()]) > 1:
                ligand_second_atom_neighbors = [atom for atom in second_atom.GetNeighbors() if atom.GetSymbol() != 'H' and atom.GetIdx() in ligand_index_atom_map.keys()]
                second_atom_neighbors_set = set(
                    [ligand_index_atom_map[atom.GetIdx()] for atom in ligand_second_atom_neighbors if atom.GetIdx() in ligand_index_atom_map.keys()])

                # Check all atoms that are adjacent to the potential second constraint atom
                for third_atom in ligand_second_atom_neighbors:
                    third_atom_neighbors_set = set(
                        [ligand_index_atom_map[atom.GetIdx()] for atom in third_atom.GetNeighbors() if
                         atom.GetSymbol() != 'H' and atom.GetIdx() in list(ligand_index_atom_map.keys())])
                    common_two_three = third_atom_neighbors_set & second_atom_neighbors_set

                    # Only investigate potential third constraint atoms with more than one neighbor, not including second atom
                    if len(third_atom_neighbors_set - common_two_three) > 1 and third_atom.GetIdx() != current_atom_index:
                        ligand_second_atom = second_atom.GetIdx()
                        ligand_third_atom = third_atom.GetIdx()
                        return ligand_second_atom, ligand_third_atom

        # NOTE!
        # The following code was added to only select ligand constraint atoms from source fragments!!!
        # To revert back to the old (worse?) method, remove the following and switch ligand_fragment_prody to full ligand

        # If the code reaches this point, then all neighbors of the second fragment atom are dead-ends
        # That, or the second fragment atom is a dead-end
        for second_atom in ligand_contact_atom_neighbors:
            # ligand_second_atom_neighbors = [atom for atom in second_atom.GetNeighbors() if atom.GetSymbol() != 'H' and atom.GetIdx() in ligand_index_atom_map.keys()]
            ligand_second_atom = second_atom.GetIdx()

            if len([atom.GetIdx() for atom in second_atom.GetNeighbors()]) > 1:
                ligand_third_atom = list(set([q for q in ligand_index_atom_map.keys()]) - set([current_atom_index, ligand_second_atom]))[0]
                return ligand_second_atom, ligand_third_atom

        # If the code reaches this point, then the second fragment atom is a dead end
        # The fragment contains 3 atoms or is contacted at the center of a branched fragment
        # Just make the constraints a triangle with some neighbors
        ligand_third_atom = list(set([atom.GetIdx() for atom in ligand_contact_atom_neighbors]) - set([ligand_second_atom]))[0]
        return ligand_second_atom, ligand_third_atom

    def generate_residue_ligand_constraints(self):
        """
        Generate matcher constraint for a given residue-ligand interaction using information from clusters
        res1 = residue
        res2 = ligand
        :return: 
        """

        # For each conformer and all of its associated motif residues...
        for single_pose_pdb_path in pdb_check(self.single_pose_dir):
            single_pose_prody = prody.parsePDB(single_pose_pdb_path)
            current_conformer = os.path.basename(os.path.normpath(single_pose_pdb_path)).split('-')[0]

            # For each motif residue...
            # todo: this will need to be updated once we start thinking about DNA/RNA interactions
            for residue_index in range(2, single_pose_prody.getHierView().numResidues() + 1):

                single_constraint_path = os.path.join(self.user_defined_dir, 'Motifs', 'Single_Constraints')
                os.makedirs(single_constraint_path, exist_ok=True)

                constraint_atoms_dict, constraint_block = self.generate_single_constraint_block(single_pose_prody, current_conformer, residue_index)

                with open(os.path.join(single_constraint_path,
                                       '{}-{}.cst'.format(os.path.basename(os.path.normpath(current_conformer)),
                                                          constraint_atoms_dict['residue_index'])),
                          'w') as constraint_file:
                    constraint_file.write('\n'.join(constraint_block))

    def generate_single_constraint_block(self, single_pose_prody, current_conformer, residue_index,
                                         distance_tolerance_d=0.5, angle_A_tolerance_d=5, angle_B_tolerance_d=5,
                                         torsion_A_tolerance_d=5, torsion_AB_tolerance_d=5, torsion_B_tolerance_d=5,
                                         torsion_constraint_sample_number=1, angle_constraint_sample_number=1,
                                         distance_constraint_sample_number=0, use_default_tolerances=True, greasy_sampling=False):
        """
        Generate a single constraint block for one residue-ligand interaction
        :return: 
        """
        constraint_atoms_dict = self.determine_constraint_atoms_from_single_pose(single_pose_prody, current_conformer, residue_index)
        ligand = os.path.basename(os.path.normpath(self.user_defined_dir))[:3]
        ligand_prody = single_pose_prody.select('resname {}'.format(ligand)).copy()
        residue_prody = single_pose_prody.select('resnum {}'.format(residue_index)).copy()

        #########################################################################################
        # Get ideal distance, angle, and torsion values from representative motif residue prody #
        #########################################################################################

        ideal_distance = constraint_atoms_dict['contact_distance']
        # 'angle_A' is the angle Res1:Atom2 - Res1:Atom1 - Res2:Atom1
        ideal_angle_A = prody.calcAngle(
            ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][1])),
            ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][0])),
            residue_prody.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][0]))
            )
        # 'angle_B' is the angle Res1:Atom1 - Res2:Atom1 - Res2:Atom2
        ideal_angle_B = prody.calcAngle(
            ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][0])),
            residue_prody.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][0])),
            residue_prody.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][1]))
            )
        # 'torsion_A' is the dihedral Res1:Atom3 - Res1:Atom2 - Res1:Atom1 - Res2:Atom1
        ideal_torsion_A = prody.calcDihedral(
            ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][2])),
            ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][1])),
            ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][0])),
            residue_prody.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][0]))
        )
        # 'torsion_AB' is the dihedral Res1:Atom2 - Res1:Atom1 - Res2:Atom1 - Res2:Atom2
        ideal_torsion_AB = prody.calcDihedral(
            ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][1])),
            ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][0])),
            residue_prody.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][0])),
            residue_prody.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][1]))
        )
        ideal_torsion_B = prody.calcDihedral(
            ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][0])),
            residue_prody.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][0])),
            residue_prody.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][1])),
            residue_prody.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][2]))
        )

        if use_default_tolerances == False:
            # todo: get rid of this, imminently obsolete
            ##########################################################
            # Import cluster residues for cluster-derived tolerances #
            ##########################################################

            # Import residues from clusters
            cluster_residue_prody_list = []
            fragment_cluster_path = os.path.join(self.user_defined_dir, 'Cluster_Results',
                                                 constraint_atoms_dict['fragment'])

            for cluster_residue in pdb_check(fragment_cluster_path, base_only=True):
                if cluster_residue.split('-')[0] == constraint_atoms_dict['cluster']:
                    cluster_residue_prody = prody.parsePDB(os.path.join(fragment_cluster_path, cluster_residue))
                    residue_type = cluster_residue_prody.getResnames()[0]

                    cluster_residue_prody_list.append(cluster_residue_prody)

            #################################################################################################
            # Get tolerance from clusters; I'm going to try making the tolerance +/- 1 SD of cluster values #
            #################################################################################################

            if len(cluster_residue_prody_list) > 3:
                distance_list = [
                    prody.calcDistance(
                        ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][0])),
                        motif.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][0])))
                    for motif in cluster_residue_prody_list]
                distance_tolerance_SD = np.std(distance_list)

                # 'angle_A' is the angle Res1:Atom2 - Res1:Atom1 - Res2:Atom1
                angle_A_list = [
                    prody.calcAngle(
                        ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][1])),
                        ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][0])),
                        motif.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][0])))
                    for motif in cluster_residue_prody_list]
                angle_A_tolerance_SD = np.std(angle_A_list)

                # 'angle_B' is the angle Res1:Atom1 - Res2:Atom1 - Res2:Atom2
                angle_B_list = [
                    prody.calcAngle(
                        ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][0])),
                        motif.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][0])),
                        motif.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][1])))
                    for motif in cluster_residue_prody_list]
                angle_B_tolerance_SD = np.std(angle_B_list)

                # 'torsion_A' is the dihedral Res1:Atom3 - Res1:Atom2 - Res1:Atom1 - Res2:Atom1
                torsion_A_list = [
                    prody.calcDihedral(
                        ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][2])),
                        ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][1])),
                        ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][0])),
                        motif.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][0])),
                        ) for motif in cluster_residue_prody_list]
                torsion_A_tolerance_SD = np.std(torsion_A_list)

                # 'torsion_AB' is the dihedral Res1:Atom2 - Res1:Atom1 - Res2:Atom1 - Res2:Atom2
                torsion_AB_list = [
                    prody.calcDihedral(
                        ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][1])),
                        ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][0])),
                        motif.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][0])),
                        motif.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][1])),
                        ) for motif in cluster_residue_prody_list]
                torsion_AB_tolerance_SD = np.std(torsion_AB_list)

                # 'torsion_B' is the dihedral Res1:Atom1 - Res2:Atom1 - Res2:Atom2 - Res2:Atom3
                torsion_B_list = [
                    prody.calcDihedral(
                        ligand_prody.select('index {}'.format(constraint_atoms_dict['ligand']['atom_indices'][0])),
                        motif.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][0])),
                        motif.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][1])),
                        motif.select('index {}'.format(constraint_atoms_dict['residue']['atom_indices'][2])),
                        ) for motif in cluster_residue_prody_list]
                torsion_B_tolerance_SD = np.std(torsion_B_list)

                # Set lower and upper bounds for tolerances
                distance_tolerance = 0.5 if distance_tolerance_SD > 0.5 else distance_tolerance_SD
                angle_A_tolerance = (360 / (
                    2 * angle_constraint_sample_number + 1) * angle_constraint_sample_number) if angle_A_tolerance_SD > 120 else angle_A_tolerance_SD
                angle_B_tolerance = (360 / (
                    2 * angle_constraint_sample_number + 1) * angle_constraint_sample_number) if angle_B_tolerance_SD > 120 else angle_B_tolerance_SD
                torsion_A_tolerance = (360 / (
                    2 * torsion_constraint_sample_number + 1) * torsion_constraint_sample_number) if torsion_A_tolerance_SD > 120 else torsion_A_tolerance_SD
                torsion_B_tolerance = (360 / (
                    2 * torsion_constraint_sample_number + 1) * torsion_constraint_sample_number) if torsion_B_tolerance_SD > 120 else torsion_B_tolerance_SD
                torsion_AB_tolerance = (360 / (
                    2 * torsion_constraint_sample_number + 1) * torsion_constraint_sample_number) if torsion_AB_tolerance_SD > 120 else torsion_AB_tolerance_SD

            else:
                distance_tolerance = distance_tolerance_d
                angle_A_tolerance = angle_A_tolerance_d
                angle_B_tolerance = angle_B_tolerance_d
                torsion_A_tolerance = torsion_A_tolerance_d
                torsion_AB_tolerance = torsion_AB_tolerance_d
                torsion_B_tolerance = torsion_B_tolerance_d

        else:
            distance_tolerance = distance_tolerance_d
            angle_A_tolerance = angle_A_tolerance_d
            angle_B_tolerance = angle_B_tolerance_d
            torsion_A_tolerance = torsion_A_tolerance_d
            torsion_AB_tolerance = torsion_AB_tolerance_d
            torsion_B_tolerance = torsion_B_tolerance_d

        # Updated this for backbone contacts to allow for any residue identity
        if constraint_atoms_dict['residue']['atom_names'][0] in ['C', 'CA', 'N', 'O', 'OXT']:
            residue_resname = 'ACDEFHIKLMNQRSTVWY'
            residue_tag = 'residue1'
        else:
            residue_resname = residue_prody.getResnames()[0]
            residue_tag = 'residue3'
        ligand_resname = ligand_prody.getResnames()[0]

        # Increase tolerance/sampling by 5/1 for greasy residues (ACFILMVWY)
        if greasy_sampling and residue_resname in ['ALA', 'CYS', 'PHE', 'ILE', 'LEU', 'MET', 'VAL', 'TRP', 'TYR']:
            angle_A_tolerance += 5
            angle_B_tolerance += 5
            torsion_A_tolerance += 5
            torsion_AB_tolerance += 5
            torsion_B_tolerance += 5
            torsion_constraint_sample_number += 1
            angle_constraint_sample_number += 1

        constraint_block = [
            '  TEMPLATE::   ATOM_MAP: 1 atom_name: {}'.format(' '.join(constraint_atoms_dict['ligand']['atom_names'])),
            '  TEMPLATE::   ATOM_MAP: 1 residue3: {}\n'.format(ligand_resname),
            '  TEMPLATE::   ATOM_MAP: 2 atom_name: {}'.format(' '.join(constraint_atoms_dict['residue']['atom_names'])),
            '  TEMPLATE::   ATOM_MAP: 2 {}: {}\n'.format(residue_tag, residue_resname),
            '  CONSTRAINT:: distanceAB: {0:7.2f} {1:6.2f} {2:6.2f}       0  {3:3}'.format(
                ideal_distance, distance_tolerance, 100, distance_constraint_sample_number),
            '  CONSTRAINT::    angle_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                float(ideal_angle_A), angle_A_tolerance, 100, angle_constraint_sample_number),
            '  CONSTRAINT::    angle_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                float(ideal_angle_B), angle_B_tolerance, 100, angle_constraint_sample_number),
            '  CONSTRAINT::  torsion_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                float(ideal_torsion_A), torsion_A_tolerance, 100,
                torsion_constraint_sample_number),
            '  CONSTRAINT::  torsion_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                float(ideal_torsion_B), torsion_B_tolerance, 100,
                torsion_constraint_sample_number),
            '  CONSTRAINT:: torsion_AB: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                float(ideal_torsion_AB), torsion_AB_tolerance, 100,
                torsion_constraint_sample_number)
            ]

        return constraint_atoms_dict, constraint_block

    def BW_constraints_from_gurobi_solutions(self, gurobi_solutions_csv, conformer, use_solution_residues_only=True):
        """
        Generate xml constraint files for Gurobi solutions using Brian Weitzner's combinatorial matcher.
        
        :param gurobi_solutions_csv: path to csv containing solutions generated by gurobi
        :param conformer: string specifying current conformer e.g. MEH_0001
        :param use_solution_residues_only: if true, only use cluster residues that are members of solutions. Otherwise,
        use all residues from clusters in matcher file
        :return: 
        """
        #################################
        # Import all solutions from csv #
        #################################

        gurobi_solutions = pd.read_csv(gurobi_solutions_csv, usecols=['Residue_indicies'])
        pprint.pprint(gurobi_solutions)

        # todo: ALSO NEED TO CONSIDER RESIDUE IDENTITY!!!
        # todo: alter res_idx_map_df so that all of this information is already available in the csv
        for index, row in self.res_idx_map_df.iterrows():
            pdb_name_split = re.split('-|_', row['source_pdb'])
            cluster_pdb = '-'.join(row['source_pdb'].split('-')[1:])
            self.res_idx_map_df.loc[index, 'fragment'] = pdb_name_split[1]
            self.res_idx_map_df.loc[index, 'cluster'] = pdb_name_split[3]
            self.res_idx_map_df.loc[index, 'resname'] = prody.parsePDB(os.path.join(self.user_defined_dir,
                                                                                    'Cluster_Results',
                                                                                    'Fragment_{}'.format(pdb_name_split[1]),
                                                                                    cluster_pdb)
                                                                       ).getResnames()[0]

        #########################################################
        # Consolidate solutions based on residue source cluster #
        #########################################################

        # Make set([source clusters]) E all solutions
        cluster_list_set = set()
        # Keep track of residues used in solutions
        cluster_residue_set = set()

        # For each solution...
        for solution in gurobi_solutions['Residue_indicies']:
            # Use residue_index_mapping to find source cluster for all residues excluding index 1 (ligand)
            relevant_rows = self.res_idx_map_df.loc[self.res_idx_map_df['residue_index'].isin(literal_eval(solution)[1:]) & (self.res_idx_map_df['source_conformer'] == conformer)]

            # Add tuple of (fragment, cluster, resname) to cluster_list_set
            cluster_list_set.add(tuple([(row['fragment'], row['cluster'], row['resname']) for index, row in relevant_rows.iterrows()]))

            # if use_solution_residues_only, keep track of residues (df index in self.res_idx_map_df)
            if use_solution_residues_only:
                for index, row in relevant_rows.iterrows():
                    cluster_residue_set.add(index)

        print("Solutions contain {} unique cluster solutions consisted of {} residues".format(len(cluster_list_set), len(cluster_residue_set)))

        ############################################################
        # Generate matcher xml for each unique cluster combination #
        ############################################################

        single_pose_prody = prody.parsePDB(os.path.join(self.single_pose_dir, '{}-single_pose.pdb'.format(conformer)))
        fragment_cluster_groups = self.res_idx_map_df.groupby(['fragment', 'cluster', 'resname'])
        bw_gurobi_constraints_dirpath = os.path.join(self.user_defined_dir, 'BW_Matcher_Constraints')
        os.makedirs(bw_gurobi_constraints_dirpath, exist_ok=True)

        # todo: would make it faster if I stuck all the residue_row_indices in a dict an just did lookups...
        # For each unqiue combination of cluster solutions
        for solution_tuples in sorted(list(cluster_list_set)):
            # Create constraint file
            contraint_file_name = '-'.join(['_'.join([toop[0], toop[1], self.three_to_one[toop[2]]]) for toop in solution_tuples]) + '.xml'
            bw_constraint_file = open(os.path.join(bw_gurobi_constraints_dirpath, contraint_file_name), 'w')

            print(contraint_file_name)

            # For each cluster within a solution
            for solution_tuple in solution_tuples: # (fragment, cluster)

                bw_constraint_file.write('<MatcherConstraint>\n')

                # Get upstream and downstream three-letter IDs
                downstream_residue = os.path.basename(os.path.normpath(self.user_defined_dir))[:3]
                upstream_residue = solution_tuple[2]

                # Get residues that will constitute a constraint block
                if use_solution_residues_only:
                    residue_row_indices = list(set(fragment_cluster_groups.get_group(solution_tuple).index) & cluster_residue_set)
                else:
                    residue_row_indices = list(fragment_cluster_groups.get_group(solution_tuple).index)

                # Get concensus contacts for each residue in cluster
                residue_contact_list = []
                ligand_contact_list = []

                for cluster_residue in residue_row_indices:
                    residue_index = fragment_cluster_groups.get_group(solution_tuple).loc[cluster_residue, 'residue_index']
                    constraint_atoms_dict = self.determine_constraint_atoms_from_single_pose(single_pose_prody, conformer, residue_index, verbose=False)
                    residue_contact_list.append(tuple(constraint_atoms_dict['residue']['atom_names']))
                    ligand_contact_list.append(tuple(constraint_atoms_dict['ligand']['atom_names']))

                # Get most frequently selected matcher atoms
                residue_contact_count = collections.Counter(residue_contact_list)
                ligand_contact_count = collections.Counter(ligand_contact_list)

                residue_contact_atoms = residue_contact_count.most_common()[0][0]
                ligand_contact_atoms = ligand_contact_count.most_common()[0][0]

                # Add to constraint block
                bw_constraint_file.write('\t<UpstreamResidue atom1="{0}" atom2="{1}" atom3="{2}" name="{3}"/>\n'.format(residue_contact_atoms[0], residue_contact_atoms[1], residue_contact_atoms[2], upstream_residue))
                bw_constraint_file.write('\t<DownstreamResidue atom1="{0}" atom2="{1}" atom3="{2}" name="{3}"/>\n'.format(ligand_contact_atoms[0], ligand_contact_atoms[1], ligand_contact_atoms[2], downstream_residue))

                # For each residue in the cluster...
                for cluster_residue in residue_row_indices:
                    residue_index = fragment_cluster_groups.get_group(solution_tuple).loc[cluster_residue, 'residue_index']

                    # ideal_distance = constraint_atoms_dict['contact_distance']
                    ideal_distance = prody.calcDistance(
                        single_pose_prody.select('resnum {} and name {}'.format(1, ligand_contact_atoms[0])), # constraint_atoms_dict['ligand']['atom_indices'][0])),
                        single_pose_prody.select('resnum {} and name {}'.format(residue_index, residue_contact_atoms[0])), # constraint_atoms_dict['residue']['atom_indices'][0]))
                    )

                    # 'angle_A' is the angle Res1:Atom2 - Res1:Atom1 - Res2:Atom1
                    ideal_angle_A = prody.calcAngle(
                        single_pose_prody.select('resnum {} and name {}'.format(1, ligand_contact_atoms[1])), #constraint_atoms_dict['ligand']['atom_indices'][1])),
                        single_pose_prody.select('resnum {} and name {}'.format(1, ligand_contact_atoms[0])), #constraint_atoms_dict['ligand']['atom_indices'][0])),
                        single_pose_prody.select('resnum {} and name {}'.format(residue_index, residue_contact_atoms[0])), #constraint_atoms_dict['residue']['atom_indices'][0]))
                    )
                    # 'angle_B' is the angle Res1:Atom1 - Res2:Atom1 - Res2:Atom2
                    ideal_angle_B = prody.calcAngle(
                        single_pose_prody.select('resnum {} and name {}'.format(1, ligand_contact_atoms[0])), #constraint_atoms_dict['ligand']['atom_indices'][0])),
                        single_pose_prody.select('resnum {} and name {}'.format(residue_index, residue_contact_atoms[0])), #constraint_atoms_dict['residue']['atom_indices'][0])),
                        single_pose_prody.select('resnum {} and name {}'.format(residue_index, residue_contact_atoms[1])), #constraint_atoms_dict['residue']['atom_indices'][1]))
                    )
                    # 'torsion_A' is the dihedral Res1:Atom3 - Res1:Atom2 - Res1:Atom1 - Res2:Atom1
                    ideal_torsion_A = prody.calcDihedral(
                        single_pose_prody.select('resnum {} and name {}'.format(1, ligand_contact_atoms[2])), #constraint_atoms_dict['ligand']['atom_indices'][2])),
                        single_pose_prody.select('resnum {} and name {}'.format(1, ligand_contact_atoms[1])), #constraint_atoms_dict['ligand']['atom_indices'][1])),
                        single_pose_prody.select('resnum {} and name {}'.format(1, ligand_contact_atoms[0])), #constraint_atoms_dict['ligand']['atom_indices'][0])),
                        single_pose_prody.select('resnum {} and name {}'.format(residue_index, residue_contact_atoms[0])), #constraint_atoms_dict['residue']['atom_indices'][0]))
                    )
                    # 'torsion_AB' is the dihedral Res1:Atom2 - Res1:Atom1 - Res2:Atom1 - Res2:Atom2
                    ideal_torsion_AB = prody.calcDihedral(
                        single_pose_prody.select('resnum {} and name {}'.format(1, ligand_contact_atoms[1])), #constraint_atoms_dict['ligand']['atom_indices'][1])),
                        single_pose_prody.select('resnum {} and name {}'.format(1, ligand_contact_atoms[0])), #constraint_atoms_dict['ligand']['atom_indices'][0])),
                        single_pose_prody.select('resnum {} and name {}'.format(residue_index, residue_contact_atoms[0])), #constraint_atoms_dict['residue']['atom_indices'][0])),
                        single_pose_prody.select('resnum {} and name {}'.format(residue_index, residue_contact_atoms[1])), #constraint_atoms_dict['residue']['atom_indices'][1]))
                    )
                    # 'torsion_B' is the dihedral Res1:Atom1 - Res2:Atom1 - Res2:Atom2 - Res2:Atom3
                    ideal_torsion_B = prody.calcDihedral(
                        single_pose_prody.select('resnum {} and name {}'.format(1, ligand_contact_atoms[0])), #constraint_atoms_dict['ligand']['atom_indices'][0])),
                        single_pose_prody.select('resnum {} and name {}'.format(residue_index, residue_contact_atoms[0])), #constraint_atoms_dict['residue']['atom_indices'][0])),
                        single_pose_prody.select('resnum {} and name {}'.format(residue_index, residue_contact_atoms[1])), #constraint_atoms_dict['residue']['atom_indices'][1])),
                        single_pose_prody.select('resnum {} and name {}'.format(residue_index, residue_contact_atoms[2])), #constraint_atoms_dict['residue']['atom_indices'][2]))
                    )

                    # todo: add options for tolerances and other goodies
                    # Write to constraint file
                    bw_constraint_file.write('\t<Combination>\n')
                    bw_constraint_file.write('\t\t<!--{}-->\n'.format(fragment_cluster_groups.get_group(solution_tuple).loc[cluster_residue, 'source_pdb']))
                    bw_constraint_file.write('\t\t<DistanceAB x0="{0}" xtol="0.00" k="0" periodicity="0.00" noSamples="0"/>\n'.format(ideal_distance[0]))
                    bw_constraint_file.write('\t\t<AngleA x0="{0}" xtol="0.00" k="0" periodicity="360.00" noSamples="0"/>\n'.format(ideal_angle_A[0]))
                    bw_constraint_file.write('\t\t<AngleB x0="{0}" xtol="0.00" k="0" periodicity="360.00" noSamples="0"/>\n'.format(ideal_angle_B[0]))
                    bw_constraint_file.write('\t\t<TorsionA x0="{0}" xtol="0.00" k="0" periodicity="360.00" noSamples="0"/>\n'.format(ideal_torsion_A[0]))
                    bw_constraint_file.write('\t\t<TorsionAB x0="{0}" xtol="0.00" k="0" periodicity="360.00" noSamples="0"/>\n'.format(ideal_torsion_AB[0]))
                    bw_constraint_file.write('\t\t<TorsionB x0="{0}" xtol="0.00" k="0" periodicity="360.00" noSamples="0"/>\n'.format(ideal_torsion_B[0]))
                    bw_constraint_file.write('\t</Combination>\n')

                bw_constraint_file.write('</MatcherConstraint>\n')

            bw_constraint_file.close()

    def conventional_constraints_from_gurobi_solutions(self, gurobi_solutions_csv_dir, constraints_to_generate=1000,
                                                       offset=0, angle_dihedral_tolerance=5, angle_dihedral_sample_number=1,
                                                       iteration=False, greasy_sampling=False, json_output=False):
        """
        Generates conventional matcher constraints for a given Gurobi solution
        
        :param gurobi_solutions_csv_dir: path to directory containing csv solutions generated by gurobi
        :param constraints_to_generate: number of (top ranking) constraint files to generate from gurobi solutions. Ignores `constraints_to_generate` parameter.
        :param offset: number of top-ranking constraints to skip before starting to generate constraint files
        :param angle_dihedral_tolerance: tolerance value for matcher angle and dihedral constraints
        :param angle_dihedral_sample_number: number of samples within tolerance for matcher angle and dihedral constraints
        :param iteration: if iteration, then generate constraints for all solutions without consolidating scores.
        :param greasy_sampling: +5 tolerance and +1 sampling for hydrophobic residues
        :return: 
        """

        def _populate_constraint_json():
            """
            Populate JSON with all relevant constraint blocks...

            Constraint_JSON[Conformer][Residue_Index] = `constraint_block_string`

            :return:
            """
            for index, row in conformer_df.iterrows():
                index_list = literal_eval(row['Residue_indicies'])

                # Add constraint blocks to new constraint file for each residue in motif minus ligand
                for residue_index in index_list[1:]:

                    if residue_index not in previously_calculated_constraint_blocks.keys():
                        print('Calculating constraint block for {} - {}'.format(current_conformer, residue_index))
                        constraint_atoms_dict, constraint_block = self.generate_single_constraint_block(
                            conformer_fuzzball, current_conformer, residue_index,
                            angle_A_tolerance_d=angle_dihedral_tolerance, angle_B_tolerance_d=angle_dihedral_tolerance,
                            torsion_A_tolerance_d=angle_dihedral_tolerance,
                            torsion_AB_tolerance_d=angle_dihedral_tolerance,
                            torsion_B_tolerance_d=angle_dihedral_tolerance,
                            torsion_constraint_sample_number=angle_dihedral_sample_number,
                            angle_constraint_sample_number=angle_dihedral_sample_number,
                            greasy_sampling=greasy_sampling)
                        previously_calculated_constraint_blocks[residue_index] = '\n'.join(constraint_block)

            json_dict[current_conformer] = previously_calculated_constraint_blocks

        def _generate_constraint_files():
            """
            Generate a constraint file...
            :return:
            """
            # For each unique solution for a given conformer
            for index, row in conformer_df.iterrows():
                index_list = literal_eval(row['Residue_indicies'])

                # Open new constraint file
                constraint_filename = '{}-{}.cst'.format(current_conformer, '_'.join([str(a) for a in index_list]))
                current_constraint_file = open(os.path.join(gurobi_constraints_path, constraint_filename), 'w')

                # Add constraint blocks to new constraint file for each residue in motif minus ligand
                for residue_index in index_list[1:]:

                    if residue_index in previously_calculated_constraint_blocks.keys():
                        constraint_block = previously_calculated_constraint_blocks[residue_index]
                    else:
                        print('Calculating constraint block for {} - {}'.format(current_conformer, residue_index))
                        constraint_atoms_dict, constraint_block = self.generate_single_constraint_block(
                            conformer_fuzzball, current_conformer, residue_index,
                            angle_A_tolerance_d=angle_dihedral_tolerance, angle_B_tolerance_d=angle_dihedral_tolerance,
                            torsion_A_tolerance_d=angle_dihedral_tolerance,
                            torsion_AB_tolerance_d=angle_dihedral_tolerance,
                            torsion_B_tolerance_d=angle_dihedral_tolerance,
                            torsion_constraint_sample_number=angle_dihedral_sample_number,
                            angle_constraint_sample_number=angle_dihedral_sample_number,
                            greasy_sampling=greasy_sampling)
                        previously_calculated_constraint_blocks[residue_index] = constraint_block

                    current_constraint_file.write('CST::BEGIN\n')
                    current_constraint_file.write('\n'.join(constraint_block))
                    current_constraint_file.write('\n')
                    current_constraint_file.write('CST::END\n')

                # Close file
                current_constraint_file.close()

                # Write motif to PDB
                # binding_motif = conformer_fuzzball.select('resnum {}'.format(' '.join([str(a) for a in index_list])))
                if not iteration:
                    binding_motif = conformer_fuzzball.select('resnum 1')

                    for resnum in index_list[1:]:
                        binding_motif = binding_motif + conformer_fuzzball.select('resnum {}'.format(resnum))

                    motif_pdb_filename = '{}-{}.pdb'.format(current_conformer, '_'.join([str(a) for a in index_list]))
                    prody.writePDB(os.path.join(gurobi_motif_path, motif_pdb_filename), binding_motif)

        # Set things up
        gurobi_things_path = os.path.join(self.user_defined_dir, 'Gurobi_Constraints')
        gurobi_constraints_path = os.path.join(gurobi_things_path, 'Constraints')
        gurobi_motif_path = os.path.join(gurobi_things_path, 'Motif_PDBs')

        os.makedirs(gurobi_things_path, exist_ok=True)
        os.makedirs(gurobi_constraints_path, exist_ok=True)
        os.makedirs(gurobi_motif_path, exist_ok=True)

        # Import all solutions from csv
        gurobi_solutions = pd.DataFrame(columns=['Obj_score', 'Residue_indicies', 'Conformer'])

        for solution_set in os.listdir(gurobi_solutions_csv_dir):
            if solution_set.endswith('.csv'):
                # todo: only try to open csv files!!! Ran into issue with hidden files...
                temp_solution_df = pd.read_csv(os.path.join(gurobi_solutions_csv_dir, solution_set), usecols=['Obj_score', 'Residue_indicies', 'Conformer'])
                gurobi_solutions = gurobi_solutions.append(temp_solution_df, ignore_index=True)

        if not iteration:
            gurobi_solutions = gurobi_solutions.sort_values(by=['Obj_score'], ascending=True).head(n=constraints_to_generate + offset).tail(constraints_to_generate)

        if json_output:
            json_dict = {}

        # Group dataframe by conformer
        for current_conformer, conformer_df in gurobi_solutions.groupby(['Conformer']):

            # Open pdb only once!
            conformer_fuzzball = prody.parsePDB(os.path.join(self.single_pose_dir, '{}-single_pose.pdb'.format(current_conformer)))

            # Keep track of previously calculated constraint blocks
            previously_calculated_constraint_blocks = {}

            # Generate constraint files or JSON
            _populate_constraint_json() if json_output else _generate_constraint_files()

        if json_output:
            json_output_path = os.path.join(gurobi_constraints_path, '{}-constraint_blocks.json'.format(self.user_defined_dir))
            with open(json_output_path, 'w') as muh_json:
                json.dump(json_dict, muh_json)

    @staticmethod
    def calculate_transformation_matrix(reference_conformer_atoms, target_conformer_atoms):
        """
        Caluclate transformation matrix by hand...
        When applied to target, transformation matrix transforms target onto reference.

        :param reference_conformer_atoms: target (in ProDy speak)
        :param target_conformer_atoms: mobile (in ProDy speak)
        :return:
        """
        # Calculate centroid for selected atoms
        reference_centroid = np.sum(reference_conformer_atoms, axis=0) / len(reference_conformer_atoms)
        target_centroid = np.sum(target_conformer_atoms, axis=0) / len(target_conformer_atoms)

        # Generate centered coordinate sets
        reference_centered = np.asarray(
            [atom - reference_centroid for atom in reference_conformer_atoms])
        target_centered = np.asarray([atom - target_centroid for atom in target_conformer_atoms])

        # Calculate Covariance Matrix
        covariance_matrix = np.dot(reference_centered.T, target_centered)

        # Singular Value Decomposition of covariance matrix
        V, s, Ut = np.linalg.svd(covariance_matrix)

        # Calculate determinant of Ut.T and V.T
        det = np.linalg.det(np.dot(Ut.T, V.T))
        det_matrix = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, det]])

        # Calculate Rotation Matrix
        rotation_matrix = np.dot(np.dot(Ut.T, det_matrix), V.T)

        # Translation vector
        translation_vector = np.dot(rotation_matrix, -(reference_centroid[np.newaxis].T)) + target_centroid[
            np.newaxis].T

        return prody.Transformation(rotation_matrix, np.ravel(translation_vector))

class Generate_Fuzzball_using_PyRosetta(object):
    """
    20190109 - Just did a bunch of refactoring to improve everything leading up to the fuzzball generation step. Most
    pertinent to the implementation of this new class is that clusters are now written to disk as .ag.npz files
    where all residues are labelled with their source using the setData() method.

    PyRosetta will allow for fine-tune control of fuzzball composition based on REU and other relevant metrics (e.g.
    total fraction of H-bonding residues, fraction of H-bonding residues for each donor/acceptor on the target ligand).
    In addition, PyRosetta will filter out all residues that do not make favorable interactions with the ligand before
    scoring with the FeatureReporter. The previous implementation could not filter these residues out based on geometry
    alone, and most of the data in the FeatureReporter database was garbage as a result.

    Clusters are composed of residues that have passed relatively crude filters. These clusters will be loaded one at a
    time by PyRosetta as a compromise between loading residues individually (loading over and over...) and loading all
    residues at once (unnecessary residue-residue edge scores). All we want is REU for each ligand-residue contact.
    From there, we can take residues that pass our REU cut-off and compose the fuzzball as we please.

    Logic:

    load reference conformer (conformer used to generate user-defined fragments) as prody object
    load all cluster prody objects into a dict

    for each conformer:

        load conformer as prody object

        for all clusters from fragments:

            for each residue in cluster:

                determine contact atoms on ligand and residue for reference conformer
                calculate transformation for ligand contact atoms between reference conformer and current conformer
                apply transformation to residue

            load transformed residues and current conformer as Pose in PyRosetta (prody > PDBStream > PyRosetta)
            score pose
            save relevant metrics into a dataframe (contact REU, hbond energy, hbond contact atom on ligand)

        decide which residues to use based on metrics from dataframe
        pull residues into a single prody

        write prody to disk as PDB


    """

    resnames = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
                'MET', 'ASP', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']

    def __init__(self, user_defined_dir, config_dict=None):

        # --- Paths --- #
        self.user_defined_dir = user_defined_dir
        self.user_config = config_dict
        self.rosetta_inputs = os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs')
        self.defined_interactions = os.path.join(self.user_defined_dir, 'Inputs', 'Defined_Interactions')

        # --- Other Variables --- #
        self.chemical_component_identifier = os.path.normpath(os.path.basename(self.user_defined_dir))[:3]
        self.reference_conformer_prody = None  # Initialized in _init_load_reference_conformer_and_clusters()
        self.cluster_dict = None  # Initialized in _init_load_reference_conformer_and_clusters()
        self.hbond_atom_set = None  # Initialized in _init_get_ligand_hbond_atoms()

        self.motif_residue_attributes_df = None

        # todo: fix this jankiness
        self.generate_constraints = Generate_Constraints(self.user_defined_dir)
        
        # --- Initialization Functions --- #
        self._init_load_reference_conformer_and_clusters()
        self._init_get_ligand_hbond_atoms()

    def _init_load_reference_conformer_and_clusters(self):
        """
        Loads everything I need to memory for generating fuzzballs
        :return:
        """
        # Load reference conformer as prody object
        # Assumes the reference conformer is {self.chemical_component_identifier}_0001.pdb
        self.reference_conformer_name = f'{self.chemical_component_identifier}_0001'
        path_to_reference_conformer = os.path.join(self.rosetta_inputs, f'{self.reference_conformer_name}.pdb')
        self.reference_conformer_prody = prody.parsePDB(path_to_reference_conformer)

        cluster_results_dir = os.path.join(self.user_defined_dir, 'Cluster_Results')
        temp_cluster_dict = dict()

        print('\nLoading clusters...\n')
        for fragment in os.listdir(cluster_results_dir):
            temp_cluster_dict[fragment] = dict()
            fragment_cluster_dir = os.path.join(cluster_results_dir, fragment)

            for cluster in pdb_check(fragment_cluster_dir, extension='.ag.npz'):
                cluster_number = int(re.split('[_.]', os.path.normpath(os.path.basename(cluster)))[1])
                temp_cluster_dict[fragment][cluster_number] = prody.loadAtoms(os.path.join(cluster))

        self.cluster_dict = temp_cluster_dict

    def _init_get_ligand_hbond_atoms(self):
        """
        Use RDKit to get the names of ligand atoms that can mediate hydrogen bond interactions.

        :return:
        """
        ligand_mol = RDKit_Mol_from_ProDy(self.reference_conformer_prody)

        # --- Hydrogen bond donor/acceptor definitions --- #

        # Pulled from rdkit/Chem/Lipinski.py
        # Lipinski module will report number of hydrogen bond donor/acceptors, but not their identities

        # 2 definitions adapted from those in the Gobbi Paper
        #  NOTE: if you want traditional Lipinski numbers, you
        #  should use NOCounts (below) instead of HAcceptor
        #
        HDonorSmarts = Chem.MolFromSmarts('[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]')
        # changes log for HAcceptorSmarts:
        #  v2, 1-Nov-2008, GL : fix amide-N exclusion; remove Fs from definition
        HAcceptorSmarts = Chem.MolFromSmarts('[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0,o,s;+0])]')

        donors = ligand_mol.GetSubstructMatches(HDonorSmarts)
        acceptors = ligand_mol.GetSubstructMatches(HAcceptorSmarts)

        # Pull atom names from ProDy using indicies (should be unchanged... knock on wood)
        # No reason to discriminate between donors/acceptors at the moment, all go into one set
        self.hbond_atom_set = set()

        for hbd in donors:
            self.hbond_atom_set.add(self.reference_conformer_prody[hbd[0]].getName())

        for hba in acceptors:
            self.hbond_atom_set.add(self.reference_conformer_prody[hba[0]].getName())

        print(f'Hydrogen bond donors/acceptors identified: {" ".join(self.hbond_atom_set)}\n')

    def score_motif_conformer_interactions(self):
        """
        Performs the bulk of work for generating fuzzballs
        :return:
        """

        import pyrosetta
        from pyrosetta import rosetta

        # todo: decide whether to save a bunch of dataframes or consolidate... dataframe is currently written only after all residues for all conformers are scored
        # List of dicts for residue information
        fuzzball_list_of_dicts = list()

        # For each conformer
        for conformer in pdb_check(self.rosetta_inputs):

            # Load conformer as Prody object
            conformer_prody = prody.parsePDB(conformer)

            # Iterate through self.cluster_dict (for loops all the way down...)
            for fragment in self.cluster_dict:

                # Save transformed clusters in a dict()
                transformed_cluster_dict = dict()

                for cluster in self.cluster_dict[fragment]:

                    current_cluster = self.cluster_dict[fragment][cluster]

                    # Make copy of conformer for current cluster
                    current_conformer_and_transformed_residues = conformer_prody.copy()

                    # --- Transform cluster residues onto target conformer --- #

                    for residue in current_cluster.iterResidues():

                        # Skip residues if essential backbone atoms are missing
                        if len({'C', 'CA', 'N', 'O'} & set([atom.getName() for atom in residue])) != 4:
                            continue

                        transformed_motif, constraint_atoms_dict = self.transform_residue_about_current_conformer(conformer_prody, residue, fragment=fragment)

                        # Add residue to current_conformer_and_transformed_residues collection for scoring
                        if transformed_motif is not None:
                            current_conformer_and_transformed_residues += transformed_motif

                    # Save all transformed cluster residues to transformed_cluster_dict
                    if current_conformer_and_transformed_residues.numResidues() > 1:
                        transformed_cluster_dict[cluster] = current_conformer_and_transformed_residues
                    else:
                        print(f'\nAll motif residues from {fragment} Cluster {cluster} were rejected!\n')

                # --- Score all clusters for current fragment using PyRosetta and record scores/attributes --- #

                # Initialize PyRosetta with stuff
                conformer_name = os.path.basename(os.path.normpath(conformer)).split('.')[0]
                current_params_path = os.path.join(self.rosetta_inputs, f'{conformer_name}.params')

                # core.conformation.Conformation muted, OXT is missing from all motif residues
                my_options = [f"-extra_res_fa {current_params_path}",
                              "-mute core.conformation core.chemical"]

                pyrosetta.init(options=' '.join(my_options))
                sfxn = rosetta.core.scoring.get_score_function()

                for transformed_cluster in transformed_cluster_dict:

                    transformed_cluster_prody = transformed_cluster_dict[transformed_cluster]

                    transformed_residues_only_hv = transformed_cluster_prody.select(f'not resname {self.chemical_component_identifier}').getHierView()
                    ligand_prody = transformed_cluster_prody.select(f'resname {self.chemical_component_identifier} and not hydrogen')

                    # Convert prody of transformed cluster into string
                    # import_pose.pose_from_pdb_stream() uses import_pose.pose_from_pdbstring() anyways and doesn't require BS options
                    transformed_cluster_stream = io.StringIO()
                    prody.writePDBStream(transformed_cluster_stream, transformed_cluster_prody)

                    # Load transformed cluster into PyRosetta
                    transformed_cluster_pose = rosetta.core.pose.Pose()
                    rosetta.core.import_pose.pose_from_pdbstring(transformed_cluster_pose, transformed_cluster_stream.getvalue())
                    sfxn(transformed_cluster_pose)

                    # Assert only the ligand was removed from selection...
                    if transformed_residues_only_hv.numResidues() + 1 != transformed_cluster_pose.size():
                        print("Hold up. What.")
                        print(transformed_residues_only_hv.numResidues() + 1, transformed_cluster_pose.size())

                        for res in transformed_residues_only_hv.iterResidues():
                            print(res)

                        for i in range(1, transformed_cluster_pose.size() + 1):
                            print(transformed_cluster_pose.residue(i))

                        raise SystemExit

                    e_edge = transformed_cluster_pose.energies().energy_graph()

                    # Record attributes for all residues in a transformed cluster
                    for index, residue_prody in enumerate(transformed_residues_only_hv.iterResidues(), start=2):

                        residue_attribute_dict = dict()

                        # Identification attributes
                        residue_attribute_dict['conformer'] = conformer_name
                        residue_attribute_dict['fragment'] = residue_prody.copy().getData('fragment_id')[0]
                        residue_attribute_dict['cluster'] = residue_prody.copy().getData('cluster')[0]
                        residue_attribute_dict['contact_source'] = residue_prody.copy().getData('contact_source')[0]
                        residue_attribute_dict['resname'] = residue_prody.getResname()

                        # Energy-related attributes
                        what = e_edge.find_energy_edge(index, 1).fill_energy_map()

                        total_score = sum([what[rosetta.core.scoring.fa_rep] * 0.55,
                                           what[rosetta.core.scoring.fa_atr],
                                           what[rosetta.core.scoring.hbond_sc],
                                           what[rosetta.core.scoring.fa_sol],
                                           what[rosetta.core.scoring.fa_elec]
                                           ]
                                          )

                        residue_attribute_dict['total_score'] = total_score

                        hbond_sc_score = what[rosetta.core.scoring.hbond_sc]
                        residue_attribute_dict['hbond_sc'] = hbond_sc_score

                        # Report contact atom in self.hbond_atom_set if hbondsc < 0, else None
                        ligand_prody_stripped = ligand_prody.select('not hydrogen')
                        residue_prody_stripped = residue_prody.select('not hydrogen')

                        ligand_residue_distance_matrix = prody.buildDistMatrix(ligand_prody_stripped.getCoords(),
                                                                               residue_prody_stripped.getCoords())

                        if hbond_sc_score < 0:

                            # Identify atom on residue mediating contact (residue_contact_atom)
                            residue_contact_index = np.argmin(np.min(ligand_residue_distance_matrix, axis=0))  # Row (Residue)

                            # Identify which atom in self.hbond_atom_set is closest to residue_contact_atom
                            ligand_contact_distances = ligand_residue_distance_matrix[:, residue_contact_index]

                            defined_contacts = [(atom.getName(), dist) for atom, dist in zip(ligand_prody_stripped, ligand_contact_distances) if atom.getName() in self.hbond_atom_set]
                            defined_contacts.sort(key=lambda x: x[1])
                            defined_contact_atom = defined_contacts[0][0]

                            residue_attribute_dict['ligand_contact'] = defined_contact_atom

                        else:
                            ligand_contact_index = np.argmin(np.min(ligand_residue_distance_matrix, axis=1))  # Column (Ligand)
                            residue_attribute_dict['ligand_contact'] = ligand_prody_stripped[ligand_contact_index].getName()

                        fuzzball_list_of_dicts.append(residue_attribute_dict)

        temp_df = pd.DataFrame(fuzzball_list_of_dicts)
        temp_df.drop_duplicates(subset=['conformer', 'contact_source', 'hbond_sc', 'ligand_contact', 'resname', 'total_score'], keep='first', inplace=True)
        self.motif_residue_attributes_df = temp_df

        # todo: should really be dumping into SQLite3 database...
        self.motif_residue_attributes_df.to_csv(os.path.join(self.user_defined_dir, 'motif_residue_attributes.csv'), index=False)

    def transform_residue_about_current_conformer(self, conformer_prody, residue, fragment=None, reference_conformer=None):
        """
        Transforms a given motif residue about a given conformer to maintain contact geometries with reference fragment
        as observed in PDB. (How to English?)

        :param conformer_prody: prody object of target conformer
        :param residue: prody object of motif residue to transform about conformer_prody
        :param fragment: Name of current Fragment as string ('Fragment_1'). Uses entire ligand to determine constraint
        atoms if None.
        :param reference_conformer:
        :return:
        """
        # --- RIPPED FORM GENERATE_MOTIF_RESIDUES.GENERATE_FUZZBALL() --- #

        deepcopy_residue = residue.copy()

        # Check whether current residue can be parsed using RDKit, skip if returns None
        if RDKit_Mol_from_ProDy(deepcopy_residue) is None:
            print(f'RDKit was unable to parse the current residue: {deepcopy_residue.getData("contact_source")[0]}')
            return None

        if not reference_conformer:
            reference_conformer = self.reference_conformer_prody

        # Get fragment atoms for transformation of mobile onto reference fragment
        # Use reference conformer to get constraint atoms! Residues were pulled from PDB using fragments in the same
        # coordinate frame as the reference conformer.
        constraint_atoms_dict = self.generate_constraints.determine_constraint_atoms(deepcopy_residue,
                                                                                     reference_conformer,
                                                                                     fragment=fragment, verbose=False)
        ligand_constraint_atoms = constraint_atoms_dict['ligand']['atom_names']

        # todo: manually append atom coords to a list, select goes by order in file and may not be the same for reference and target
        # Select atoms from reference ligand conformer (conformer #1 e.g. MEH_0001)
        reference_conformer_coord_list = [reference_conformer.select(f'name {atom}').getCoords()[0] for atom in ligand_constraint_atoms]
        reference_conformer_atoms = np.asarray(reference_conformer_coord_list)

        # Select atoms from current conformer
        target_conformer_coord_list = [conformer_prody.select(f'name {atom}').getCoords()[0] for atom in ligand_constraint_atoms]
        target_conformer_atoms = np.asarray(target_conformer_coord_list)

        prody.writePDB('reference_conformer_atoms', reference_conformer.select('name {}'.format(' '.join(ligand_constraint_atoms))))
        prody.writePDB('target_conformer_atoms', conformer_prody.select('name {}'.format(' '.join(ligand_constraint_atoms))))

        transformation_matrix = Generate_Constraints.calculate_transformation_matrix(reference_conformer_atoms, target_conformer_atoms)
        return prody.applyTransformation(transformation_matrix, deepcopy_residue), constraint_atoms_dict

    def define_fuzzball_composition(self, conformer):
        """
        Define fuzzball composition based on importance.

        Add user-defined residues
        Add top scoring packing residues
        Add residues specifically for defined hbonding atoms
        Evenly distribute residues around ligand based on contact atom

        :param conformer: name of current conformer for which residues are being selected
        :param top_greasy: number of top scoring hydrophobic residues to include in the nucleation fuzzball before
        adding hydrogen bonding residues only.
        :return:
        """
        # Import previously generated motif_residue_attributes_df, else generate it
        # todo: update this when I migrate to an SQLite3 database
        motif_residue_attributes_csv = os.path.join(self.user_defined_dir, 'motif_residue_attributes.csv')
        if os.path.exists(motif_residue_attributes_csv):
            print('motif_residue_attributes.csv found...\n')
            self.motif_residue_attributes_df = pd.read_csv(motif_residue_attributes_csv)

        else:
            print('motif_residue_attributes.csv not found, generating it now...\n')
            self.score_motif_conformer_interactions()

        # --- Decide fuzzball composition --- #
        """
        20190117.
        I'm going to decide this on the fly for now, but I need to make sure that hydrophobic residues do not dominate
        the set of solutions that Gurobi provides. I think it is important (at least for the nucleation step) to include
        only a small subset of the top-scoring hydrophobic residues and then populate the rest of the fuzzball with my
        hydrogen bonding residues. This will make the inital fuzzball smaller, but should provide more viable solutions
        for iteration and design by maximizing the number of hydrogen bond donors/acceptors in the defined binding 
        site. This means that determining composition will be dynamic, which may be an issue down the line...
        """
        # Dict of contact sources to use for selections
        fuzzball_motif_residues = dict()
        fuzzball_motif_residues['conformer'] = conformer

        # Decide how many residues to take for each group
        hbond_contact_count = 150
        packing_contact_count = 50

        # Apply banned residue filter to motif_residue_attributes_df
        banned_resnames = {'ALA', 'GLY', 'PRO', 'CYS'}
        working_resnames = list(set(self.resnames) - banned_resnames)
        self.motif_residue_attributes_df = self.motif_residue_attributes_df[self.motif_residue_attributes_df['resname'].isin(working_resnames)]

        # Get residues for current conformer
        current_conformer_motifs = self.motif_residue_attributes_df.groupby('conformer').get_group(conformer)
        current_conformer_motifs_by_contact = current_conformer_motifs.groupby('ligand_contact')
        motif_residues_aggr_list = list()

        # --- Hydrogen Bonding Residues --- #

        fuzzball_motif_residues['hbond_motifs'] = {}

        for contact_atom in self.hbond_atom_set:
            contacts_df = current_conformer_motifs_by_contact.get_group(contact_atom)
            hbond_contacts_df = contacts_df[(contacts_df['hbond_sc'] < 0) & (contacts_df['total_score'] < 0)]

            # Take best scoring contacts and add to fuzzball_motif_residues
            best_hbond_contacts = hbond_contacts_df.nsmallest(hbond_contact_count, 'total_score')
            fuzzball_motif_residues['hbond_motifs'][contact_atom] = best_hbond_contacts['contact_source'].values
            motif_residues_aggr_list += list(best_hbond_contacts['contact_source'].values)

        # --- Packing Residues --- #

        fuzzball_motif_residues['packing_motifs'] = {}

        # Packed atoms are just atoms in the ligands that are not in defined_hbonding_atoms
        packed_atoms = set(current_conformer_motifs['ligand_contact'].values) - set(self.hbond_atom_set)

        for packed_atom in packed_atoms:
            contacts_df = current_conformer_motifs_by_contact.get_group(packed_atom)
            good_contacts_df = contacts_df[(contacts_df['hbond_sc'] == 0) & contacts_df['total_score'] < 0]

            # Take best scoring contacts and add to fuzzball_motif_residues
            best_packing_contacts = good_contacts_df.nsmallest(packing_contact_count, 'total_score')
            fuzzball_motif_residues['packing_motifs'][packed_atom] = best_packing_contacts['contact_source'].values
            motif_residues_aggr_list += list(best_packing_contacts['contact_source'].values)

        # --- Motif Residue Records --- #
        fuzzball_motif_df = self.motif_residue_attributes_df[self.motif_residue_attributes_df['contact_source'].isin(motif_residues_aggr_list)]

        return fuzzball_motif_residues, fuzzball_motif_df

    def assemble_defined_fuzzball(self, conformer, limit=1200):
        """
        Yiss.

        Output both ag.npz and .pdb representations of fuzzball. PDB representation doesn't save residue source info.

        :param conformer: Name of current conformer e.g. 'MEH_0001'
        :param limit: number of residues to include in a fuzzball
        :return:
        """

        # --- User-Defined Residues --- #

        # Just iterate through .pdb files in {user_defined_dir}/Inputs/Defined_Interactions

        # --- Hydrogen Bonding and Packing Motif Residues --- #

        # Get list of unique motif residue IDs (organized by contact atom)
        motif_dict, motif_residue_summary_df = self.define_fuzzball_composition(conformer)
        unique_motif_id_set = set(motif_residue_summary_df['contact_source'].values)

        # Pull residues (deepcopy!) out of clusters and store in a list
        motif_residue_list_untransformed = list()

        for fragment in self.cluster_dict:
            for cluster in self.cluster_dict[fragment]:

                current_cluster = self.cluster_dict[fragment][cluster]
                residues_in_cluster = set(current_cluster.getData('contact_source'))
                residues_to_pull = set(residues_in_cluster) & unique_motif_id_set

                if len(residues_to_pull) > 0:
                    selection_string = ' '.join([f'`{residue_to_pull}`' for residue_to_pull in residues_to_pull])
                    pulled_residue_prody = current_cluster.select(f'contact_source {selection_string}').getHierView()
                    for residue in pulled_residue_prody.iterResidues():
                        motif_residue_list_untransformed.append(residue)

        # --- Transform all residues about target ligand --- #

        current_conformer = prody.parsePDB(os.path.join(self.rosetta_inputs, f'{conformer}.pdb'))
        current_conformer_fuzzball = current_conformer.copy()

        # --- Add motif residues to fuzzball --- #

        # Add user-defined residues first
        offset = user_defined_offset = 2  # Place user-defined residues in beginning of fuzzball

        print('Processing defined interactions...\n')
        for index, defined_interaction in enumerate(pdb_check(os.path.join(self.user_defined_dir, 'Inputs', 'Defined_Interactions')), start=user_defined_offset):

            print(f'Working on {os.path.basename(os.path.normpath(defined_interaction))}...\n')

            defined_interaction_prody= prody.parsePDB(defined_interaction)
            defined_motif_residue = defined_interaction_prody.select('protein').copy()
            defined_ligand = defined_interaction_prody.select(f'resname {self.chemical_component_identifier}').copy()

            # --- Superpose defined_interaction ligand onto current_conformer, then use current_conformer to get constraints --- #

            # Super jank, WILL NOT WORK with ligand possessing multiple conformations...
            # Should be (and tried) mapping atom names of defined onto current, but there are issues with symmetry

            """
            Okay here's how this is going to work:
            1.  Iterate through all available conformers and calculate RMSD (conformer should exist in my library, or
                else I'm messing up...)
            2.  Take conformer with lowest RMSD to defined ligand and map reference atom names directly onto defined.
                Probably going to iterate through defined atoms and adopt name of closest reference atom
                I tried using RDKit substructure search but there are issues with symmetry
            3.  Use defined interaction to get constraint atoms, save and apply for all conformers
            """

            # Identify lowest RMSD conformer in library
            def find_defined_conformer_in_library(defined_ligand):
                rmsd = 100  # lowest RMSD of a conformer to defined
                conformer_prody = None  # Path to conformer in self.rosetta_inputs

                for conf_path in pdb_check(self.rosetta_inputs):
                    conf_prody = prody.parsePDB(conf_path)
                    current_rmsd = prody.calcRMSD(defined_ligand, conf_prody)
                    if current_rmsd < rmsd:
                        rmsd = current_rmsd
                        conformer_prody = conf_prody

                return conformer_prody

            # Map conformer heavy atoms onto defined ligand
            def map_conformer_atom_names_onto_defined(defined, conf):
                defined_ligand_mapped = defined.copy()

                # Build distance matrix between conformer and defined
                distance_matrix = prody.buildDistMatrix(defined, conf)

                for index, atom in enumerate(defined):

                    if atom.getElement() != 'H':
                        # Get index of minimum distance of defined to conformer
                        min_index = np.argmin(distance_matrix[index])
                        print(f'Defined atom {atom.getName():<4} maps to conformer atom {conf[min_index].getName():<4}')

                        mapped_name = conf[min_index].getName()
                        defined_ligand_mapped[index].setName(mapped_name)

                return defined_ligand_mapped

            conformer_with_lowest_rmsd_to_defined = find_defined_conformer_in_library(defined_ligand)

            # Superpose defined onto conformer
            superposed_defined, superpose_transformation = prody.superpose(defined_ligand, conformer_with_lowest_rmsd_to_defined)

            # Map conformer atom names onto defined (defined_ligand {rows}, conformer_with_lowest_rmsd_to_defined {columns})
            defined_ligand_mapped = map_conformer_atom_names_onto_defined(superposed_defined, conformer_with_lowest_rmsd_to_defined)

            # --- Apply transformation to motif residue --- #

            # Residue needs to be in same reference frame as the superposed defined ligand before determining transformation
            superposed_motif = prody.applyTransformation(superpose_transformation, defined_motif_residue)

            motif_residue_transformed, constraint_atoms_dict = self.transform_residue_about_current_conformer(current_conformer, superposed_motif, reference_conformer=defined_ligand_mapped)

            # Tag defined motif residue with input filename
            defined_interaction_filename = os.path.basename(os.path.normpath(defined_interaction)).split('.')[0]
            motif_residue_transformed.setData('contact_source', [defined_interaction_filename] * len(motif_residue_transformed))
            motif_residue_transformed.setData('defined_contact', [True] * len(motif_residue_transformed))

            # Set resnum, occupany, coordset
            motif_residue_transformed.setResnums([index] * len(motif_residue_transformed))
            motif_residue_transformed.setChids(['A'] * len(motif_residue_transformed))
            motif_residue_transformed.setAltlocs([None] * len(motif_residue_transformed))
            motif_residue_transformed.setOccupancies([1] * len(motif_residue_transformed))

            # Add defined motif to fuzzball
            current_conformer_fuzzball += motif_residue_transformed
            offset += 1

        # Maintain dict of transformed residues (keys resname) to check for motif residue redundancy
        transformed_motifs = dict.fromkeys(self.resnames, [])

        # Add all other motif residues to fuzzball
        for index, motif_residue in enumerate(motif_residue_list_untransformed, start=offset):

            motif_residue_transformed, constraint_atoms_dict = self.transform_residue_about_current_conformer(current_conformer, motif_residue)

            # --- Redundant residue check --- #
            current_resname = motif_residue_transformed.getResnames()[0]
            current_motif_constraint_atoms = constraint_atoms_dict['residue']['atom_names']
            current_motif_ligand_contact = constraint_atoms_dict['ligand']['atom_names'][0]

            current_motif_coords = motif_residue_transformed.select(f'name {" ".join(current_motif_constraint_atoms)}').getCoords()
            current_motif_ca = motif_residue_transformed.select('name CA')

            if any([len(current_motif_coords) < 3, len(current_motif_ca) !=1]):
                print(f'Essential atom coordinates are missing from {motif_residue_transformed}')
                continue

            truthiness_list = list()

            for transformed_motif_tuple in transformed_motifs[current_resname]:
                constraint_atoms, ligand_contact_atom, motif_contact_coords, motif_ca = transformed_motif_tuple

                res_contact_atoms_match = current_motif_constraint_atoms == constraint_atoms
                ligand_contact_atom_match = current_motif_ligand_contact == ligand_contact_atom

                if res_contact_atoms_match:
                    contact_atoms_rmsd = prody.calcRMSD(motif_contact_coords, current_motif_coords) < 0.5
                    ca_distance = prody.calcDistance(motif_ca, current_motif_ca) < 1.5
                else:
                    contact_atoms_rmsd = False
                    ca_distance = False

                truthiness_list.append(all([res_contact_atoms_match, ligand_contact_atom_match, contact_atoms_rmsd, ca_distance]))

            if any(truthiness_list):
                # print(f'\n{motif_residue_transformed.getData("contact_source")[0]} was found to be redundant!\n')
                continue

            # print(f'\n{motif_residue_transformed.getData("contact_source")[0]} accepted!\n')
            # Set resnum, occupany, coordset
            motif_residue_transformed.setResnums([index] * len(motif_residue_transformed))
            motif_residue_transformed.setAltlocs([None] * len(motif_residue_transformed))
            motif_residue_transformed.setOccupancies([1] * len(motif_residue_transformed))

            # Add residue to fuzzball
            current_conformer_fuzzball += motif_residue_transformed

            # Remember transformed motifs for
            transformed_motifs[current_resname].append((current_motif_constraint_atoms, current_motif_ligand_contact, current_motif_coords, current_motif_ca))

        # --- Output fuzzballs --- #
        fuzzball_dir = os.path.join(self.user_defined_dir, 'Fuzzballs')
        os.makedirs(fuzzball_dir, exist_ok=True)
        current_conformer_fuzzball_path = os.path.join(fuzzball_dir, f'{conformer}-fuzzball')

        # .ag.npz representation (for iteration and contact sourcing)
        prody.saveAtoms(current_conformer_fuzzball, filename=current_conformer_fuzzball_path)

        # .pdb representation (for FeatureReporter)
        prody.writePDB(current_conformer_fuzzball_path, current_conformer_fuzzball)

        # Motif residue summary dataframe
        motif_residue_summary_df.to_csv(os.path.join(fuzzball_dir, f'{conformer}-fuzzball.csv'), index=False)

    # todo: actually use this somewhere...
    def validate_user_defined_fuzzball_residues(self):
        """
        Parse user-defined residue interactions with the ligand to include in the fuzzball. This function will expect a
        ligand-residue pair for each defined interaction. The ligand can be in any conformation as the atom names
        assigned by Rosetta will be mapped onto the ligand. The residue just needs to be there completely so that it can
        by imported by Rosetta (i.e. includes all relevant backbone atoms C, CA, N, O)

        :return: list of defined residue interactions with the ligand
        """

        defined_interaction_list = list()

        # Load reference ligand as RDKit Mol object
        # todo: define a project configuration file and pull values from there...
        ligand_reference_mol = Chem.MolFromPDBFile(os.path.join(self.rosetta_inputs, '38E_0001.pdb'), removeHs=True)

        for defined_interaction in pdb_check(self.defined_interactions):

            defined_interaction_prody = prody.parsePDB(defined_interaction)

            # --- Check that ligand is complete --- #

            # Pull ligand from defined_interaction and convert to Mol object
            defined_interaction_ligand = defined_interaction_prody.select(f'name {self.chemical_component_identifier}')
            defined_interaction_ligand_pdbstring = io.StringIO()
            prody.writePDBStream(defined_interaction_ligand_pdbstring, defined_interaction_ligand)
            ligand_defined_mol = Chem.MolFromPDBFile(defined_interaction_ligand_pdbstring.getvalue(), removeHs=True)

            # Map all heavy atoms from defined ligand onto reference
            ligand_mcs = rdFMCS.FindMCS([ligand_defined_mol, ligand_reference_mol], bondCompare=rdFMCS.BondCompare.CompareAny)
            ligand_substructure_mol = Chem.MolFromSmarts(ligand_mcs.smartsString)

            ligand_defined_match = ligand_defined_mol.GetSubstructMatch(ligand_substructure_mol)
            ligand_reference_match = ligand_reference_mol.GetSubstructMatch(ligand_substructure_mol)

            # Assert 1:1 mapping excluding hydrogens
            if len(ligand_defined_match) != len([a for a in ligand_reference_match.GetAtoms()]):
                print('\nFull ligand was not mapped!\n')
                continue

            # --- Check that residue is complete --- #

            defined_interaction_residue = defined_interaction_prody.select(f'not resname {self.chemical_component_identifier}')

            if len(defined_interaction_residue.numResidues() > 1):
                print('There can only be one residue in a defined interaciton!')
                continue
            elif len(defined_interaction_residue.numResidues() < 1):
                print('A residue is not defined!')
                continue

            # Assert that essential backbone atoms are present
            defined_residue_bb_atoms = [a.getName() for a in defined_interaction_residue.select('backbone')]
            if len(set(defined_residue_bb_atoms) & {'C', 'CA', 'N', 'O'}) != 4:
                print('Essential backbone atoms are missing from the defined residue!')

        return defined_interaction_list

# DEPRECIATED
class Generate_Motif_Residues(Generate_Constraints):  # todo: Generate_Motif_Residues SHOULD NOT be inheriting from Generate_Constraints
    """
    This is a class for generating and manipulating motif residues for hypothetical ligand binding sites
    
    For each ligand conformer during motif residue generation, I will calculate and output:
    1.  A sparse matrix of motif residues that clash
    2.  Residue-ligand interaction scores as calculated with Rosetta
    3.  Constraint file for each residue-ligand interaction
    """
    # todo: accommodate C-term residues... could just add 1 to everything and say either or... sloppy though
    # todo: replace with PyRosetta filter...
    expected_atoms = {'ALA': 5, 'CYS': 6, 'ASP': 8, 'GLU': 9, 'PHE': 11, 'GLY': 4, 'HIS': 10, 'ILE': 8,
                      'LYS': 9, 'LEU': 8, 'MET': 8, 'ASN': 8, 'PRO': 7, 'GLN': 9, 'ARG': 11, 'SER': 6,
                      'THR': 7, 'VAL': 7, 'TRP': 14, 'TYR': 12, 'MSE': 8, 'SEC': 6}

    def __init__(self, user_defined_dir, config_dict=None):
        """
        :param user_defined_dir: 
        :param motif_cluster_yaml: 
        :param generate_single_pose: Generate a single pose containing ligand and all motif residues?
        """
        Generate_Constraints.__init__(self, user_defined_dir)
        self.user_defined_dir = user_defined_dir
        self.motif_residue_dir = os.path.join(self.user_defined_dir, 'Motifs', 'Representative_Residue_Motifs')
        self.residue_ligand_interactions_dir = os.path.join(self.user_defined_dir, 'Motifs', 'Residue_Ligand_Interactions')
        self.score_interactions_list_path = os.path.join(self.residue_ligand_interactions_dir, 'PDBs_to_score.txt')
        self.rosetta_inputs_path = os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs')
        self.conformer_transformation_dict = None
        self.user_config = config_dict

    def _fix_mse_sec(self, representative_residue, resname):
        """
        Takes MSE/SEC and turns it into MET/CYS
        :param representative_residue: representative_residue with Se
        :return: representative_residue with S
        """
        # Find index of SE
        res_elements = representative_residue.getElements()
        seleno_index = [e for e in res_elements].index('SE')

        # Set SE to S
        res_elements[seleno_index] = 'S'

        # Set elements to MET
        representative_residue.setElements(res_elements)

        # Set resnames MSE>MET, SEC>CYS
        # Set seleno atom name to met/cys atom sulfur atom name
        if resname == 'MSE':
            representative_residue.setResnames(['MET'] * len(representative_residue))
            res_atom_names = representative_residue.getNames()
            res_atom_names[seleno_index] = 'SD'
            representative_residue.setNames(res_atom_names)
        elif resname == 'SEC':
            representative_residue.setResnames(['CYS'] * len(representative_residue))
            res_atom_names = representative_residue.getNames()
            res_atom_names[seleno_index] = 'SG'
            representative_residue.setNames(res_atom_names)
        else:
            raise Exception('Either MSE or SEC had to be set!')

        return representative_residue

    def generate_residue_ligand_pdbs(self, conformer, motif_residue_list,
                                     generate_residue_ligand_pairs=False, generate_single_pose=False, generate_reference_pose=False):
        """
        Generate residue-ligand PDBs used for scoring. Uses residues in <motif_residue_list> to populate a fuzzball based
        on whatever filters are in vogue...

        :param conformer: path to pdb of conformer
        :param motif_residue_list: list of (motif_filename, motif_prody) tuples
        :return:
        """
        # Write all PDBs to a text file so I can score them all with score_jd2 -in:file:l <text_file_with_list.txt>
        os.makedirs(self.residue_ligand_interactions_dir, exist_ok=True)

        ligand_ID = os.path.basename(os.path.normpath(conformer)).split('.')[0]
        ligand_prody = prody.parsePDB(conformer)

        # DEPRECIATED
        if generate_residue_ligand_pairs == True:
            conformer_residue_output_dir = os.path.join(self.residue_ligand_interactions_dir, ligand_ID)
            os.makedirs(conformer_residue_output_dir, exist_ok=True)

            # For each motif residue in Representative_Residue_Motifs directory, combine ligand and residue
            for motif_residue_name, residue_prody in motif_residue_list:
                residue_ligand_prody = residue_prody + ligand_prody
                motif_residue_index = motif_residue_name.split('-')[0]

                # Output ligand-residue pairs to a new Resdiue_Ligand_Interactions directory under Inputs
                pdb_output_path = os.path.join(conformer_residue_output_dir, '{}-{}.pdb'.format(ligand_ID, motif_residue_index))
                prody.writePDB(pdb_output_path, residue_ligand_prody)
                with open(self.score_interactions_list_path, 'a+') as pdb_list:
                    pdb_list.write('{}\n'.format(pdb_output_path))

        # Generate single pose of ligand and all motif residues if option is turned on (default)
        if generate_single_pose == True:

            # Make a convenient directories
            if generate_reference_pose:
                os.makedirs(self.reference_pose_dir, exist_ok=True)
            else:
                single_pose_dir = os.path.join(self.residue_ligand_interactions_dir, "Single_Poses")
                os.makedirs(single_pose_dir, exist_ok=True)

            # # Touch file for list of single_pose paths
            # open(os.path.join(single_pose_dir, 'single_pose_list.txt'), 'w').close()

            # Directory and file paths
            single_pose_dir = os.path.join(self.residue_ligand_interactions_dir, "Single_Poses")
            single_pose_list = os.path.join(single_pose_dir, 'single_pose_list.txt')

            # Make a residue explosion
            conformer_single_pose = prody.parsePDB(conformer)
            conformer_single_pose_clash_reference = prody.parsePDB(conformer)

            # Generate/append dict mapping residue indicies to source (will be exported/appended to .csv)
            residue_list_of_dicts = []

            # Count residues added (not enumerating b/c I already wrote most of this...)
            # todo: this will need to be fixed when we move to designing DNA/RNA interactions
            residue_count = 2 # Ligand is 1

            for motif_residue_name, residue_prody in motif_residue_list:

                if not generate_reference_pose:
                    # Don't add residues that Rosetta will throw out e.g. missing backbone atoms
                    essential_bb_atoms = ['C', 'CA', 'N']

                    # If a residue is too close to the ligand, continue
                    contact_distance, row_index, column_index = minimum_contact_distance(residue_prody, conformer_single_pose_clash_reference, return_indices=True)

                    if contact_distance <= 1.5 or not all([atom in residue_prody.getNames() for atom in essential_bb_atoms]):
                        continue

                # Renumber all residues and set all residues to the same chain
                residue_prody.setResnums([residue_count] * len(residue_prody))
                residue_prody.setChids(['X'] * len(residue_prody))

                # Add info to list of dicts
                residue_list_of_dicts.append({'struct_id': int(ligand_ID.split('_')[1]),
                                              'source_conformer': ligand_ID,
                                              'residue_index': residue_count,
                                              'source_pdb': motif_residue_name
                                              })

                # Add to residue explosion
                conformer_single_pose = conformer_single_pose + residue_prody
                residue_count += 1

            # Output residue explosion to it's own happy directory under Inputs
            if generate_reference_pose == True:
                pdb_output_path = os.path.join(self.reference_pose_dir, '{}-single_pose.pdb'.format(ligand_ID))
            else:
                pdb_output_path = os.path.join(single_pose_dir, '{}-single_pose.pdb'.format(ligand_ID))

            prody.writePDB(pdb_output_path, conformer_single_pose)

            # Add residue index to new/existing .csv
            if generate_reference_pose == True:
                # os.makedirs(self.reference_pose_dir, exist_ok=True)
                output_csv_path = os.path.join(self.reference_pose_dir, 'residue_index_mapping.csv')
            else:
                output_csv_path = os.path.join(single_pose_dir, 'residue_index_mapping.csv')

            if os.path.exists(output_csv_path):
                df_existing = pd.read_csv(output_csv_path, index_col=0)
                df_new = pd.DataFrame(residue_list_of_dicts)
                df_concat = df_existing.append(df_new, ignore_index=True)
                df_concat.to_csv(output_csv_path)
            else:
                df_new = pd.DataFrame(residue_list_of_dicts)
                df_new.to_csv(output_csv_path)

            if generate_reference_pose == False:
                # Add PDB path to single_pose list
                with open(single_pose_list, 'a') as muh_list:
                    muh_list.write(pdb_output_path + '\n')

    def generate_all_the_fuzzballs(self, filter_redundant_residues=True, pose_from_clusters=True, use_pyrosetta=True):
        """
        Generate fuzzballs for all ligand conformers using residues in reference pose
        Generate single pose containing all residues from selected clusters for each conformer

        :param filter_redundant_residues:
        :param pose_from_clusters:
        :param use_pyrosetta: use pyrosetta to determine fuzzball residues
        :return:
        """
        fragment_prody_dict = self.create_fragment_prody_dict(pose_from_clusters=pose_from_clusters)

        # Then for each conformer...
        for conformer in pdb_check(os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs'), conformer_check=True):
            self.generate_fuzzball(conformer, fragment_prody_dict, filter_redundant_residues=filter_redundant_residues)

    def create_fragment_prody_dict(self, pose_from_clusters=True):
        """
        Create dictionary of fragments and all associated clusters for fuzzballs
        :return:
        """

        # Assemble dict to store list of prody residue instances for each cluster
        fragment_prody_dict = collections.OrderedDict()

        # For each cluster in the cluster list
        for fragment in self.fragment_cluster_list:
            fragment_prody_dict[fragment] = collections.OrderedDict()
            fragment_cluster_path = os.path.join(self.user_defined_dir, 'Cluster_Results', fragment)

            # Make a list for each cluster for a given fragment
            if pose_from_clusters:
                for cluster_number in self.fragment_cluster_list[fragment]:
                    fragment_prody_dict[fragment][int(cluster_number)] = []
            else:
                number_of_fragment_clusters = len(set([int(re.split('-|_', motif)[1]) for motif in pdb_check(fragment_cluster_path)]))
                for cluster_number in range(1, number_of_fragment_clusters + 1):
                    fragment_prody_dict[fragment][cluster_number] = []

            # Check pdbs for given fragment and add to appropriate cluster list
            for fragment_pdb in pdb_check(fragment_cluster_path, base_only=True):
                cluster_number = int(re.split('[-_.]', fragment_pdb)[1])

                if cluster_number in self.fragment_cluster_list[fragment]:

                    fragment_prody = prody.parsePDB(os.path.join(fragment_cluster_path, fragment_pdb))
                    fragment_prody_dict[fragment][cluster_number].append(
                        [os.path.join(fragment, fragment_pdb), fragment_prody])  # fragment path, fragment prody

        return fragment_prody_dict

    def generate_fuzzball(self, conformer, fragment_prody_dict, filter_redundant_residues=True, use_pyrosetta=True):
        """
        Creates a fuzzball for the given conformer
        :return:
        """
        # Reference Conformer
        reference_conformer_name = f'{self.user_defined_dir[:3]}_0001'
        reference_conformer_prody = prody.parsePDB(os.path.join(self.rosetta_inputs_path, f'{reference_conformer_name}.pdb'))

        conformer_name = os.path.basename(os.path.normpath(conformer)).split('.')[0]
        conformer_prody = prody.parsePDB(conformer)

        # Make a list of all transformed prody motif residues
        motif_residue_list = []

        # For each user-defined fragment...
        for dict_fragment in fragment_prody_dict:

            # For each "quality" cluster...
            for cluster_index in fragment_prody_dict[dict_fragment]:

                print('\n\nEvaluating {}, cluster {}...\n\n'.format(dict_fragment, cluster_index))

                accepted_residues = []

                # For each residue where (cluster_residue_name, cluster_residue_prody)...
                for cluster_member_tuple in fragment_prody_dict[dict_fragment][cluster_index]:

                    # todo: decide whether this is necessary
                    # Pass ALA/GLY/CSY/PRO residues
                    if cluster_member_tuple[1].getResnames()[0] in ['GLY', 'ALA', 'CYS', 'PRO']:
                        continue

                    # Deep copy required... applyTransformation uses pointers to residue location
                    deepcopy_residue = copy.deepcopy(cluster_member_tuple[1])

                    # todo: fix this as well... so hacky.
                    if filter_redundant_residues:

                        # Known issue with messed up residues in the PDB where neighbormap cannot be formed
                        try:
                            contraint_dict = self.determine_constraint_atoms(deepcopy_residue, conformer_prody, fragment=dict_fragment, verbose=False)

                        except Exception as e:
                            print('\n\nThere was an error determining constraint atoms for the currrent residue\n{}\n\n'.format(e))
                            continue

                        relevant_residue_CA = deepcopy_residue.select('name CA')
                        relevant_residue_atoms = deepcopy_residue.select('name {}'.format(' '.join(contraint_dict['residue']['atom_names'])))
                        relevant_residue_atoms_coords = relevant_residue_atoms.getCoords()

                        if len(accepted_residues) != 0:

                            truthiness_list = [all(
                                (contraint_dict['residue']['atom_names'] == res_feature_tuple[0],
                                 prody.calcRMSD(res_feature_tuple[1], relevant_residue_atoms_coords) < 0.75,
                                 prody.calcDistance(res_feature_tuple[2], relevant_residue_CA)[0] < 1.5)
                            ) for res_feature_tuple in accepted_residues]

                            if any(truthiness_list):
                                print('\n{0} was found to be redundant!\n'.format(deepcopy_residue))
                                continue

                            accepted_residues.append((contraint_dict['residue']['atom_names'], relevant_residue_atoms_coords, relevant_residue_CA))

                        else:
                            accepted_residues.append((contraint_dict['residue']['atom_names'], relevant_residue_atoms_coords, relevant_residue_CA))

                    # Get fragment atoms for transformation of mobile onto reference fragment
                    constraint_atoms_dict = self.determine_constraint_atoms(deepcopy_residue, reference_conformer_prody, fragment=dict_fragment, verbose=False)
                    ligand_constraint_atoms = constraint_atoms_dict['ligand']['atom_names']

                    # Select atoms from reference ligand conformer (conformer #1 e.g. MEH_0001)
                    reference_conformer_atoms = reference_conformer_prody.select(
                        'name {}'.format(' '.join(ligand_constraint_atoms))).getCoords()

                    # Select atoms from current conformer
                    target_conformer_atoms = conformer_prody.select(
                        'name {}'.format(' '.join(ligand_constraint_atoms))).getCoords()

                    transformation_matrix = self.calculate_transformation_matrix(reference_conformer_atoms, target_conformer_atoms)
                    transformed_motif = prody.applyTransformation(transformation_matrix, deepcopy_residue)

                    if use_pyrosetta:
                        motif_residue_list.append((cluster_member_tuple[0], transformation_matrix))  # Path to residue, transformation matrix

                    else:
                        motif_residue_list.append((cluster_member_tuple[0], transformed_motif))  # Path to residue, transformed residue

        if use_pyrosetta:
            self.fuzzball_pyrosetta(conformer, motif_residue_list)
        else:
            self.generate_residue_ligand_pdbs(conformer, motif_residue_list, generate_single_pose=True)

    def calculate_transformation_matrix(self, reference_conformer_atoms, target_conformer_atoms):
        """
        Caluclate transformation matrix by hand...

        :param reference_conformer_atoms:
        :param target_conformer_atoms:
        :return:
        """
        # Calculate centroid for selected atoms
        reference_centroid = np.sum(reference_conformer_atoms, axis=0) / len(reference_conformer_atoms)
        target_centroid = np.sum(target_conformer_atoms, axis=0) / len(target_conformer_atoms)

        # Generate centered coordinate sets
        reference_centered = np.asarray(
            [atom - reference_centroid for atom in reference_conformer_atoms])
        target_centered = np.asarray([atom - target_centroid for atom in target_conformer_atoms])

        # Calculate Covariance Matrix
        covariance_matrix = np.dot(reference_centered.T, target_centered)

        # Singular Value Decomposition of covariance matrix
        V, s, Ut = np.linalg.svd(covariance_matrix)

        # Calculate determinant of Ut.T and V.T
        det = np.linalg.det(np.dot(Ut.T, V.T))
        det_matrix = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, det]])

        # Calculate Rotation Matrix
        rotation_matrix = np.dot(np.dot(Ut.T, det_matrix), V.T)

        # Translation vector
        translation_vector = np.dot(rotation_matrix, -(reference_centroid[np.newaxis].T)) + target_centroid[np.newaxis].T

        return prody.Transformation(rotation_matrix, np.ravel(translation_vector))

    def fuzzball_pyrosetta(self, conformer, motif_residue_list):
        """
        20180724 - It dawned on me that I should be scoring the motif residue interaction energies at the fuzzball
        generation step... Half the residues I pass to the feature reporter aren't used by Gurobi since they have shit
        interaction energies with the ligand. Feature reporter scales very poorly (time and memory) with more residues.
        Now that I have PyRosetta working on my laptop, I should score pairwise energies for all residues within 4A of
        each fragment and accept all residues that pass my score cutoff (-0.5 REU where only considering fa_atr+fa_rep+
        fa_sol+fa_elec+sc_hbond). This should save on a lot of superfluous shit...
        :param conformer: name of current conformer
        :param motif_residue_list: list of all motif residues to consider as prody objects
        :return:
        """

        conformer_name = os.path.basename(os.path.normpath(conformer)).split('.')[0]
        current_params_path = os.path.join(self.rosetta_inputs_path, f'{conformer_name}.params')

        pyrosetta.init(options="-extra_res_fa {0} -ex1 -ex2 -extrachi_cutoff 0".format(current_params_path))
        sfxn = rosetta.core.scoring.get_score_function()

        """    
        20181204 - Another idea: once fuzzballs can be generated based on REU, I can be a little more sophisticated in 
        how I populate fuzzballs.
        
        1. For each fragment, load all residues that pass clustering in a single pose
            > Keep track of source pdb for each residue!
            > Append by jump
        2. Get ligand-residue interaction energy for all residues and rank
            > Presumably another dataframe...
        3. Decide how residues will be added to the fuzzball
            > No more clustering, just take top ~1200 residues globally or weighted % based on fragments.
            > Partition residues into greasy/polar and adjust proportions basic on hbond donors/acceptors on ligand
            > Add residues after the fact based on user-defined design constraints
            
        The nice thing is that handling residues with PyRosetta here will prevent issues with the FeatureReporter...
        I think all this logic should be separate from the conventional cluster method, which relies on a lot of crap
        that isn't necessary for a PyRosetta implementation. 
        """

        # Load current ligand conformer
        current_fuzzball = rosetta.core.pose.Pose()
        rosetta.core.import_pose.pose_from_file(current_fuzzball, conformer)
        sfxn(current_fuzzball)

        # List of dicts for residue information
        fuzzball_dict_of_dicts = {}

        # --- Generate a reference pose with all residues PROPERLY TRANSFORMED --- #

        # Iterate through residues in motif_residue_list
        for res in motif_residue_list:

            # Load residue into its own pose
            path_to_cluster_residue = res[0]
            cluster_residue = rosetta.core.pose.Pose()

            try:
                rosetta.core.import_pose.pose_from_file(cluster_residue, os.path.join(self.user_defined_dir, 'Cluster_Results', path_to_cluster_residue))
            except Exception as e:
                print(f"There was an error loading this residue into Rosetta:\n{e}")
                continue

            # Append residue to current fuzzball and transform
            current_fuzzball.append_residue_by_jump(cluster_residue.residue(1).clone(), 1)
            current_resnum = current_fuzzball.size()

            rosetta_rot_matrix = rosetta.numeric.xyzMatrix_double_t.rows(*np.ndarray.flatten(res[1].getRotation()))
            rosetta_trans_vector = rosetta.numeric.xyzVector_double_t(*res[1].getTranslation())
            current_fuzzball.residue(current_resnum).apply_transform_Rx_plus_v(rosetta_rot_matrix, rosetta_trans_vector)

            # Add info to fuzzball_list_of_dicts
            source_pdb_split = re.split('[-/]', path_to_cluster_residue)
            temp_dict = {'fragment': source_pdb_split[0],
                         'cluster': source_pdb_split[1],
                         'resname': source_pdb_split[2].split('_')[0],
                         'source_pdb': source_pdb_split[3],
                         'resnum': current_resnum}
            fuzzball_dict_of_dicts[current_resnum] = temp_dict

        # todo: dump reference fuzzballs with all possible cluster residues...

        # Score entire fuzzball
        sfxn(current_fuzzball)
        e_edge = current_fuzzball.energies().energy_graph()

        for res in range(2, current_fuzzball.size() + 1):

            what = e_edge.find_energy_edge(res, 1).fill_energy_map()

            total_score = sum([what[rosetta.core.scoring.fa_rep] * 0.55,
                               what[rosetta.core.scoring.fa_atr],
                               what[rosetta.core.scoring.hbond_sc],
                               what[rosetta.core.scoring.fa_sol],
                               what[rosetta.core.scoring.fa_elec]
                               ]
                              )

            fuzzball_dict_of_dicts[res]['total_score'] = total_score
            fuzzball_dict_of_dicts[res]['hbond_sc'] = what[rosetta.core.scoring.hbond_sc]

        df = pd.DataFrame.from_dict(fuzzball_dict_of_dicts, orient='index')
        # todo: Figure out where to dump these...
        df.to_csv(f'{conformer_name}-fuzzball_scores.csv')