#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
import prody
from rdkit import Chem
import pprint
import yaml
import io
import re
import shutil
import copy
import itertools
import subprocess
import sqlite3
import MySQLdb
import collections
import json
from ast import literal_eval
from scipy.special import comb
from pathos.multiprocessing import ProcessingPool as Pool
from .alignments import Align_PDB
from .utils import *

class Generate_Constraints():

    three_to_one = {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K',
                    'LEU':'L', 'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V',
                    'TRP':'W', 'TYR':'Y'}

    def __init__(self, user_defined_dir):
        self.user_defined_dir = user_defined_dir
        self.fragment_cluster_list = yaml.load(open(os.path.join(self.user_defined_dir, 'Inputs', 'User_Inputs', 'Motif_Clusters.yml'), 'r'))
        self.single_pose_dir = os.path.join(self.user_defined_dir, 'Motifs', 'Residue_Ligand_Interactions', 'Single_Poses')
        self.reference_pose_dir = os.path.join(self.user_defined_dir, 'Motifs', 'Reference_Pose')
        self.res_idx_map_df = False

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

        constraints_dict = self._determine_constraint_atoms(residue_prody, current_conformer, fragment)

        constraints_dict['source_pdb'] = source_pdb
        constraints_dict['residue_index'] = residue_index
        constraints_dict['fragment'] = fragment
        constraints_dict['cluster'] = cluster
        constraints_dict['resname'] = resname

        return constraints_dict

    def _determine_constraint_atoms(self, residue_prody, current_conformer, fragment, verbose=False):
        """
        Determine the three atoms in the ligand and the residue that will be used to calculate ideal values for the DOFs
        required for matcher. 

        :param residue_prody: prody object of residue
        :param current_conformer: string of current ligand conformer (e.g. "MEH_0001")
        :param fragment: string of current fragment (e.g. "Fragment_1")
        :return: 
        """

        # I need:
        # prody obect of residue
        # prody object of ligand
        # All that source PDB/fragment/cluster/resname information

        ligand_code = os.path.basename(os.path.normpath(self.user_defined_dir))[:3]

        # Get conformer prody
        # todo: fix conformer_io indicies
        # Could have done a selection from conformer_residue_prody, but there's something off about the indexing
        # So I need to reset the indexing on the selection, the selection keeps indices from original object
        # fragment_source_conformer = conformer_residue_prody.select('resname {}'.format(ligand))
        fragment_source_conformer_path = os.path.join(self.user_defined_dir,
                                                      'Inputs',
                                                      'Rosetta_Inputs',
                                                      '{}.pdb'.format(
                                                          os.path.basename(os.path.normpath(current_conformer)))
                                                      )

        fragment_source_conformer = prody.parsePDB(fragment_source_conformer_path)

        residue_index_atom_map = {atom.getIndex(): atom.getName() for atom in residue_prody.select('not hydrogen')}
        residue_atom_index_map = {v: k for k, v in residue_index_atom_map.items()}
        ligand_index_atom_map = {atom.getIndex(): atom.getName() for atom in fragment_source_conformer.select('not hydrogen')}

        # Write residue_prody to StringIO. This is required to select neighbor atoms in _determine_next_residue_constraint_atom()
        residue_io = io.StringIO()
        prody.writePDBStream(residue_io, residue_prody)

        RD_residue = Chem.MolFromPDBBlock(residue_io.getvalue(), removeHs=True)
        RD_ligand = Chem.MolFromPDBBlock(open(fragment_source_conformer_path, 'r').read(), removeHs=False)  # Trouble? (if removeHs=True; EDIT: THIS NEEDS TO BE FALSE)

        # I need to determine the closest atom-atom contacts and two additional atoms for determining bond torsions and angles
        # NOTE: Contact distance and indicies are for residue and ligand with hydrogens stripped!

        # So it occurred to me that contraints would be more meaningful if the residue was constrained to atoms in
        # the ligand from its source fragment...

        # Select source fragment atoms from current conformer
        fragment_atoms_prody = prody.parsePDB(os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs', '{}.pdb'.format(fragment)))
        fragment_atom_names = fragment_atoms_prody.select('not hydrogen').getNames()

        conformer_fragment_atoms = fragment_source_conformer.select('resname {} and name {}'.format(ligand_code, ' '.join(fragment_atom_names)))
        contact_distance, residue_index_low, ligand_index_low = minimum_contact_distance(residue_prody, conformer_fragment_atoms, return_indices=True)

        # So I derped, residue_index_low and ligand_index_low are indices from the distance matrix and do not
        # necessarily correspond to residue/ligand atom indicies... just happened to work for residues since most
        # of them have hydrogens stripped already

        residue_atom_list = [atom for atom in residue_prody.select('not hydrogen')]
        ligand_atom_list = [atom for atom in conformer_fragment_atoms.select('not hydrogen')]

        residue_contact_atom = residue_atom_list[residue_index_low]
        ligand_contact_atom = ligand_atom_list[ligand_index_low]

        # Okay, now the hard part: selecting two additional atoms downstream from the contact atoms for meaningful constraints... programmatically...
        # Using RDKit to determine neighboring atoms, removes Hydrogens by default

        # RESIDUE CONSTRAINT ATOM INDICES
        # Special cases: O (O>C>CA), N (N>CA>C), C (C>CA>N), CA (CA>C>O)
        # todo: build checks to ensure these atoms actually exist...
        # todo: get C-Term and N-term special atoms...
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

        # debugging
        # pprint.pprint(ligand_atom_list)
        # pprint.pprint(ligand_contact_atom)
        # pprint.pprint([atom.GetIdx() for atom in RD_ligand.GetAtoms()])

        ligand_second_atom, ligand_third_atom = self._determine_ligand_constraint_atoms(ligand_contact_atom.getIndex(), RD_ligand, conformer_fragment_atoms)

        if verbose:

            # So I derped, residue_index_low and ligand_index_low are indices from the distance matrix and do not
            # necessarily correspond to residue/ligand atom indicies... just happened to work for residues since most
            # of them have hydrogens stripped already

            print('RESIDUE - FIRST ATOM: {} {}'.format(residue_contact_atom.getIndex(), residue_index_atom_map[residue_contact_atom.getIndex()]))
            print('RESIDUE - SECOND ATOM: {} {}'.format(residue_second_atom, residue_index_atom_map[residue_second_atom]))
            print('RESIDUE - THIRD ATOM: {} {}'.format(residue_third_atom, residue_index_atom_map[residue_third_atom]))

            # LIGAND CONSTRAINT ATOM INDICES
            # So as far as choosing ligand constraints go, I think it's safe to chose whatever atoms as long as the
            # terminal atom has >1 neighbor. The rationale being this will propogate atom selection toward less
            # flexible parts of the ligand... hopefully...

            # todo: may need to generate multiple ligand constraints depending on where rotatable bonds lie...
            # This shouldn't be a huge issue since I generate all of the transformed ligand-residue contacts
            # Just need to consider the conformer when I concat complete constraint files together

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

                    # Need to filter for residues with missing atoms...
                    if residue_type == constraint_atoms_dict['resname'] and cluster_residue_prody.numAtoms() == \
                            self.expected_atoms[residue_type]:
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

    def conventional_constraints_from_gurobi_solutions(self, gurobi_solutions_csv_dir, constraints_to_generate=1000, offset=0, angle_dihedral_tolerance=5, angle_dihedral_sample_number=1, iteration=False, greasy_sampling=False, json_output=False):
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

    # DEPRECIATED
    def assemble_cluster_dict(self, names_only=False):
        """
        Assembles an ordered dict containing prody or names of cluster residues organized by source fragment and cluster
        index number. 
        
        :param names_only: return pdb names only if true, else prody objects for each each
        :return: 
        """
        # Assemble dict to store list of prody residue instances for each cluster
        fragment_prody_dict = collections.OrderedDict()

        # For each cluster in the cluster list
        for fragment in self.fragment_cluster_list:
            fragment_prody_dict[fragment] = collections.OrderedDict()
            fragment_cluster_path = os.path.join(self.user_defined_dir, 'Cluster_Results', fragment)

            # Make a list for each cluster for a given fragment
            for cluster_number in self.fragment_cluster_list[fragment]:
                fragment_prody_dict[fragment][int(cluster_number)] = []

            # Check pdbs for given fragment and add to appropriate cluster list
            for fragment_pdb in pdb_check(fragment_cluster_path, base_only=True):
                cluster_number = int(re.split('-|_', fragment_pdb)[1])

                if cluster_number in self.fragment_cluster_list[fragment]:
                    if names_only:
                        fragment_prody_dict[fragment][cluster_number].append(fragment_pdb)
                    else:
                        fragment_prody_dict[fragment][cluster_number].append(
                            prody.parsePDB(os.path.join(fragment_cluster_path, fragment_pdb)).select('not hydrogen'))

        return fragment_prody_dict


class Generate_Motif_Residues(Generate_Constraints):
    """
    This is a class for generating and manipulating motif residues for hypothetical ligand binding sites
    
    For each ligand conformer during motif residue generation, I will calculate and output:
    1.  A sparse matrix of motif residues that clash
    2.  Residue-ligand interaction scores as calculated with Rosetta
    3.  Constraint file for each residue-ligand interaction
    """
    # todo: accommodate C-term residues... could just add 1 to everything and say either or... sloppy though
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

    def _select_respresentative_residue(self, cluster_residues):
        """
        Select representative residues from a given cluster 
        The method for determining representative residues depends on the contact type and residue ID
        * Work in progress *
        
        :param cluster_residues: list of residues in a cluster
        :return: representative residue prody
        """
        average_coordinates = np.mean([a.getCoords() for a in cluster_residues], axis=0)

        # Select residue closest to average
        rmsd_list = [prody.calcRMSD(average_coordinates, residue.getCoords()) for residue in cluster_residues]
        representative_residue = cluster_residues[rmsd_list.index(min(rmsd_list))]


        # Check for MSE and SEC
        # Change MSE into MET, SEC into CYS
        representative_residue_resname = representative_residue.getResnames()[0]
        if representative_residue_resname == 'MSE' or representative_residue_resname == 'SEC':
            representative_residue = self._fix_mse_sec(representative_residue, representative_residue_resname)

        return representative_residue

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

    def single_pose_cluster_residue_dump(self):
        """
        Generate single pose containing all residues from selected clusters for each conformer
        :return: 
        """
        self.generate_fuzzballs(filter_redundant_residues=True, pose_from_clusters=True)

    def generate_residue_ligand_pdbs(self, conformer, motif_residue_list, generate_residue_ligand_pairs=False, generate_single_pose=False, generate_reference_pose=False):
        """
        Generate residue-ligand PDBs used for scoring
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

    def generate_fuzzballs(self, filter_redundant_residues=True, pose_from_clusters=True):
        """
        Generate fuzzballs for all ligand conformers using residues in reference pose

        :param filter_redundant_residues:
        :param pose_from_clusters:
        :return:
        """
        # 1. Generate dict of fragments and all associated clusters

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
                number_of_fragment_clusters = len(
                    set([int(re.split('-|_', motif)[1]) for motif in pdb_check(fragment_cluster_path)]))
                for cluster_number in range(1, number_of_fragment_clusters + 1):
                    fragment_prody_dict[fragment][cluster_number] = []

            # Check pdbs for given fragment and add to appropriate cluster list
            for fragment_pdb in pdb_check(fragment_cluster_path, base_only=True):
                cluster_number = int(re.split('-|_', fragment_pdb)[1])

                if cluster_number in self.fragment_cluster_list[fragment]:
                    # Check for residues with missing BB atoms
                    fragment_prody = prody.parsePDB(os.path.join(fragment_cluster_path, fragment_pdb))

                    # if all([bb_atom in fragment_prody.getNames() for bb_atom in ['C', 'CA', 'N']]):
                    fragment_prody_dict[fragment][cluster_number].append(['-'.join([fragment, fragment_pdb]), fragment_prody])

        # Reference Conformer
        reference_conformer_name = f'{self.user_defined_dir[:3]}_0001'
        reference_conformer_prody = prody.parsePDB(os.path.join(self.rosetta_inputs_path, f'{reference_conformer_name}.pdb'))

        # Then for each conformer...
        for conformer in pdb_check(os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs'), conformer_check=True):

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

                        # Deep copy required... applyTransformation uses pointers to residue location
                        deepcopy_residue = copy.deepcopy(cluster_member_tuple[1])

                        # todo: fix this as well... so hacky.
                        if filter_redundant_residues:

                            # Known issue with messed up residues in the PDB where neighbormap cannot be formed
                            try:
                                contraint_dict = self._determine_constraint_atoms(deepcopy_residue, conformer_name, dict_fragment, verbose=False)

                            except Exception as e:
                                print(
                                    '\n\nThere was an error determining constraint atoms for the currrent residue\n{}\n\n'.format(
                                        e))
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
                        constraint_atoms_dict = self._determine_constraint_atoms(deepcopy_residue, reference_conformer_name, dict_fragment, verbose=False)
                        ligand_constraint_atoms = constraint_atoms_dict['ligand']['atom_names']

                        # Select atoms from reference ligand conformer (conformer #1 e.g. MEH_0001)
                        reference_conformer_atoms = reference_conformer_prody.select('name {}'.format(' '.join(ligand_constraint_atoms))).getCoords()

                        # Select atoms from current conformer
                        target_conformer_atoms = conformer_prody.select('name {}'.format(' '.join(ligand_constraint_atoms))).getCoords()

                        transformation_matrix = self.calculate_transformation_matrix(reference_conformer_atoms, target_conformer_atoms)
                        transformed_motif = prody.applyTransformation(transformation_matrix, deepcopy_residue)

                        motif_residue_list.append((cluster_member_tuple[0], transformed_motif))

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

    # Not necessary...
    def generate_reference_pose(self):
        """
        Dumps all cluster motif residues into a reference pose

        :return: 
        """
        motif_residue_list = []
        cluster_results = os.path.join(self.user_defined_dir, 'Cluster_Results')

        for fragment in cluster_results:

            for motif_residue in pdb_check(os.path.join(cluster_results, fragment), base_only=True):
                fragment_prody = prody.parsePDB(os.path.join(os.path.join(cluster_results, fragment, motif_residue)))

                if all([bb_atom in fragment_prody.getNames() for bb_atom in ['C', 'CA', 'N']]):
                    motif_residue_list.append(('-'.join([fragment, motif_residue]), fragment_prody))

        self.generate_residue_ligand_pdbs(f'{self.user_defined_dir[:3]}_0001.pdb', motif_residue_list, generate_single_pose=True, generate_reference_pose=True)

    # Not necessary...
    def generate_motif_transformation_dict(self, pose_from_clusters=False):
        """
        Generates prody transformation matricies for each motif residue in a single pose based on its matcher constraint
        atoms.
        :param pose_from_clusters: If true, generate transformation dict for all residues in selected clusters
        :return:
        """
        # Precalculate transformation matrices for each fragment for each conformer
        ligand_code = os.path.basename(os.path.normpath(self.user_defined_dir))[:3]

        # Generate single pose with either representative motif residues or cluster residues
        motif_residue_list = []

        # Representative Residues (NOT SUPPORTED ANYMORE, OLD AF)
        if os.path.exists(self.motif_residue_dir) and not pose_from_clusters:
            for representative_residue in pdb_check(self.motif_residue_dir, base_only=True):
                motif_residue_list.append((representative_residue,
                                          prody.parsePDB(os.path.join(self.motif_residue_dir, representative_residue))
                                           )
                                          )
            conformer = os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs', '{}_0001.pdb'.format(ligand_code))
            self.generate_residue_ligand_pdbs(conformer, motif_residue_list, generate_single_pose=True, generate_reference_pose=True)

        # Single pose from cluster residues
        else:
            reference_pose_dir = os.path.join(self.user_defined_dir, 'Motifs', 'Reference_Pose')
            if os.path.exists(reference_pose_dir):
                shutil.rmtree(reference_pose_dir)
            self.generate_fuzzballs(generate_reference_pose=True)

        # Import residue_index_map after generating the first single pose
        self.res_idx_map_df = pd.read_csv(os.path.join(self.reference_pose_dir, 'residue_index_mapping.csv'))

        # Use single pose to determine constraint atoms for each motif residue, add to dict
        # For each motif residue, the transformation required to align the ligand contact atoms in the reference ligand
        # onto the current conformer is the same transformation required to maintain the motif residue contact
        ligand_constraint_atoms_dict = {}
        for single_pose_path in pdb_check(self.reference_pose_dir):
            single_pose_prody = prody.parsePDB(single_pose_path)
            single_pose_hv = single_pose_prody.getHierView()
            ligand_prody = single_pose_prody.select('resname {}'.format(ligand_code))

            current_conformer = os.path.basename(os.path.normpath(single_pose_path)).split('-')[0]

            for motif_residue in single_pose_hv.iterResidues():
                if motif_residue.getResname() != ligand_code:
                    constraint_atoms_dict = self.determine_constraint_atoms_from_single_pose(single_pose_prody, current_conformer, motif_residue.getResnum(), verbose=False)
                    ligand_constraint_atoms_dict[motif_residue.getResnum()] = constraint_atoms_dict['ligand']['atom_names']
            break # Redundant

        # todo: convert all uses of this function into for loops so I can convert this into a generator to save on memory
        # Generate conformer_transformation_dict
        conformer_transformation_dict = {}

        for conformer in pdb_check(os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs'), conformer_check=True):
            conformer_name = os.path.basename(os.path.normpath(conformer)).split('.')[0]
            conformer_transformation_dict[conformer_name] = {}
            conformer_prody = prody.parsePDB(conformer)

            # todo: just use prody here, they fixed the reflection issue...
            for key, value in ligand_constraint_atoms_dict.items():

                # Select atoms from reference ligand conformer (conformer #1 e.g. MEH_0001)
                reference_conformer_atoms = ligand_prody.select('name {}'.format(' '.join(ligand_constraint_atoms_dict[key]))).getCoords()

                # Select atoms from current conformer
                target_conformer_atoms = conformer_prody.select('name {}'.format(' '.join(ligand_constraint_atoms_dict[key]))).getCoords()

                # Calculate centroid for selected atoms
                reference_centroid = np.sum(reference_conformer_atoms, axis=0)/len(reference_conformer_atoms)
                target_centroid = np.sum(target_conformer_atoms, axis=0)/len(target_conformer_atoms)

                # Generate centered coordinate sets
                reference_centered = np.asarray([atom -  reference_centroid for atom in reference_conformer_atoms])
                target_centered = np.asarray([atom -  target_centroid for atom in target_conformer_atoms])

                # Calculate Covariance Matrix
                covariance_matrix = np.dot(reference_centered.T, target_centered)

                # Singular Value Decomposition of covariance matrix
                V, s, Ut = np.linalg.svd(covariance_matrix)

                # Calculate determinant of Ut.T and V.T
                det = np.linalg.det(np.dot(Ut.T, V.T))
                det_matrix = np.matrix([[1,0,0],[0,1,0],[0,0,det]])

                # Calculate Rotation Matrix
                rotation_matrix = np.dot(np.dot(Ut.T, det_matrix), V.T)

                # Translation vector
                translation_vector = np.dot(rotation_matrix, -(reference_centroid[np.newaxis].T)) + target_centroid[np.newaxis].T

                conformer_transformation_dict[conformer_name][key] = prody.Transformation(rotation_matrix, np.ravel(translation_vector))

        return conformer_transformation_dict

    # DEPRECIATED
    def prepare_motifs_for_conformers(self):
        """
        Prepare things for conformer binding site generation
        :return:
        """
        os.makedirs(self.residue_ligand_interactions_dir, exist_ok=True)

        # Generate transformation dict
        self.conformer_transformation_dict = self.generate_motif_transformation_dict()

        # Generate poses for representative residues only
        self.generate_residue_residue_clash_matrix(actually_generate_matrix=True)
        self.import_res_idx_map()
        self.generate_residue_ligand_clash_list()
        self.score_residue_ligand_interactions()
        self.generate_residue_ligand_constraints()

    # DEPRECIATED
    def score_residue_ligand_interactions(self):
        """
        Score each unique residue-ligand interaction with Rosetta and export to .csv (fa_atr, hbond_sc, fa_elec)
        :return:
        """
        # Calculate scores for all PDBs in this new directory (batch process, don't forget .params)
        current_ligand = os.path.basename(os.path.normpath(self.user_defined_dir))
        scores_dir = os.path.join(self.residue_ligand_interactions_dir, '{}_scores_raw.txt'.format(current_ligand))

        # todo: make a config file where users can specify Rosetta path and other things...
        run_jd2_score = subprocess.Popen([os.path.join(self.user_config['Rosetta_path'], 'main/source/bin/score_jd2.{}'.format(self.user_config['Rosetta_compiler'])),
                                          '-l',
                                          self.score_interactions_list_path,
                                          '-extra_res_fa',
                                          os.path.join(self.rosetta_inputs_path, '{}.params'.format(current_ligand)),
                                          # '-scorefile_format',
                                          # 'json',
                                          '-out:file:scorefile',
                                          scores_dir,
                                          '-out:file:score_only',
                                          '-overwrite'
                                          ])
        run_jd2_score.wait()

        # Process score file into a .csv that can be imported by numpy as a matrix
        print('Raw score file outputted as {}!'.format(scores_dir))

        score_df = pd.read_csv(scores_dir,
                               delim_whitespace=True,
                               index_col=5,
                               usecols=['hbond_sc',
                                        'hbond_bb_sc',
                                        'fa_elec',
                                        'fa_atr',
                                        'fa_rep',
                                        'description'],
                               skiprows=0,
                               header=1)

        # Output only necessary scores to a .csv
        score_df.to_csv(os.path.join(self.residue_ligand_interactions_dir, '{}_scores_df.csv'.format(current_ligand)))

    # DEPRECIATED
    def generate_residue_residue_clash_matrix(self, clashing_cutoff=2, actually_generate_matrix=True):
        """
        Generates a sparse matrix for all representative motif residues that clash with each other. This is determined
        based on the distance of the closest atom-atom interaction between two residues.
        Matrix is output as a .csv that can be imported with numpy.

        :param motif_residue_dir: directory with representative motif residue PDBs
        :param clashing_cutoff: distance cutoff for clashing residues in angstroms. Default set to 2.
        :return:
        """

        residue_residue_clash_dict = {}

        # touch text file to keep track of paths for PDBs to score
        open(self.score_interactions_list_path, 'w').close()

        # Import residue_index_mapping for reference pose
        residue_index_mapping_df = pd.read_csv(
            os.path.join(self.user_defined_dir, 'Motifs', 'Reference_Pose', 'residue_index_mapping.csv'))

        # Okay so for the actually residue-residue clashing stuff for each conformer
        # For each conformer I want to determine motif clashes with...
        for conformer in pdb_check(os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs'),
                                   conformer_check=True):
            conformer_name = os.path.basename(os.path.normpath(conformer)).split('.')[0]

            # For each residue...
            # Make a list of all transformed prody motif residues, then pass to minimum_contact_distance()
            motif_residue_list = []
            for motif_residue in pdb_check(
                    os.path.join(self.user_defined_dir, 'Motifs', 'Representative_Residue_Motifs')):
                motif_prody = prody.parsePDB(motif_residue)

                # Get residue index from residue_index_mapping_df
                residue_index_row = residue_index_mapping_df.loc[
                    residue_index_mapping_df['source_pdb'] == os.path.basename(os.path.normpath(
                        motif_residue))]  # & residue_index_mapping_df['source_conformer'] == conformer_name
                residue_index = residue_index_row['residue_index'].values[0]

                print(motif_prody, residue_index)

                # Get translation and rotation for fragment onto conformer
                transformation_matrix = self.conformer_transformation_dict[conformer_name][residue_index]

                transformed_motif = prody.applyTransformation(transformation_matrix, motif_prody)
                motif_residue_list.append((os.path.basename(os.path.normpath(motif_residue)), transformed_motif))

            # Generate residue-conformer PDBs for scoring here so I don't have to regenerate the transformed residues...
            self.generate_residue_ligand_pdbs(conformer, motif_residue_list, generate_residue_ligand_pairs=True,
                                              generate_single_pose=True, generate_reference_pose=False)

            if actually_generate_matrix == True:
                residue_residue_clash_set = set()
                for outer_index, outer_motif_tuple in enumerate(motif_residue_list):
                    for inner_index, inner_motif_tuple in enumerate(motif_residue_list[outer_index + 1:]):
                        outer_motif_index = int(outer_motif_tuple[0].split('-')[0])
                        inner_motif_index = int(inner_motif_tuple[0].split('-')[0])
                        if minimum_contact_distance(outer_motif_tuple[1], inner_motif_tuple[1]) < clashing_cutoff:

                            # If clashing, append tuples in both orders... look up times? Whatever!
                            if outer_motif_index != inner_motif_index:
                                residue_residue_clash_set.add((outer_motif_index, inner_motif_index))
                                residue_residue_clash_set.add((inner_motif_index, outer_motif_index))

                residue_residue_clash_dict[conformer_name] = list(residue_residue_clash_set)

        if actually_generate_matrix == True:
            yaml.dump(residue_residue_clash_dict, open(
                os.path.join(self.user_defined_dir, 'Inputs', 'User_Inputs', 'Residue_Residue_Clash_COO.yml'), 'w'))

    # DEPRECIATED
    def transform_motif_residues_onto_fragments(self):
        """
        Exactly what is says!

        :return:
        """
        # Precalculate transformation matrices for each fragment for each conformer
        conformer_transformation_dict = {}

        # Open fragment prodys and stick them in a dicitonary
        fragment_prody_dict = {}

        # First need to do translation and rotation of all motif residues relative to fragments in new conformations
        # Select atom names used for matcher constraints...
        for conformer in pdb_check(os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs'),
                                   conformer_check=True):
            conformer_name = os.path.basename(os.path.normpath(conformer)).split('.')[0]
            conformer_transformation_dict[conformer_name] = {}

            for fragment in pdb_check(os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs')):
                # Map fragment onto conformer
                # Add fragment prody to fragment_prody_dict
                fragment_name = os.path.basename(os.path.normpath(fragment)).split('.')[0]
                fragment_prody_dict[fragment_name] = prody.parsePDB(fragment)

                # Determine transformation for fragment onto conformer
                align = Align_PDB(self.user_defined_dir)

                # Select fragment atoms from conformer
                conformer_prody = prody.parsePDB(conformer)
                fragment_prody = prody.parsePDB(fragment)
                fragment_prody_atoms_names = fragment_prody.getNames()
                corresponding_target_atoms = conformer_prody.select(
                    'name {}'.format(' '.join(fragment_prody_atoms_names)))

                frag_atom_coords = fragment_prody.getCoords()
                trgt_atom_coords = corresponding_target_atoms.getCoords()

                # Align only on rigid atoms (if they are defined in Rigid_Fragment_Atoms dir)
                frag_inputs_dir = os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs',
                                               'Rigid_Fragment_Atoms')

                if os.path.exists(frag_inputs_dir):
                    frag_rigid_pdb_name = '{}-rigid.pdb'.format(fragment_name)

                    if frag_rigid_pdb_name in os.listdir(frag_inputs_dir):
                        frag_atom_rigid, trgt_atom_rigid = align.return_rigid_atoms(fragment_name, frag_atom_coords,
                                                                                    trgt_atom_coords)

                        transformation = prody.calcTransformation(frag_atom_rigid, trgt_atom_rigid)

                    else:
                        transformation = prody.calcTransformation(frag_atom_coords, trgt_atom_coords)

                else:
                    transformation = prody.calcTransformation(frag_atom_coords, trgt_atom_coords)

                conformer_transformation_dict[conformer_name][fragment_name] = transformation

                # Debugging: Outputs PDBs of fragment transformed onto corresponding atoms on conformer
                # import io
                # os.makedirs(os.path.join(self.user_defined_dir, 'Conformer_Tests'), exist_ok=True)
                # transformed_fragment = prody.applyTransformation(transformation, prody.parsePDBStream(io.StringIO(fragment_pdb_string)))
                # prody.writePDB(os.path.join(self.user_defined_dir, 'Conformer_Tests', '{}-{}.pdb'.format(con_name, frag_name)), transformed_fragment)

        return conformer_transformation_dict

    # DEPRECIATED
    def generate_residue_ligand_clash_list(self, cutoff_distance=2):
        """
        Filters representative motif residues that are clashing with the ligand for each conformer in Inputs/Rosetta_Inputs
        Outputs a yaml file for each conformer with a list of motif residue indicies that clash
        Clashing is determined by a distance cutoff for the closest atom-atom distance between each residue and the ligand

        :param motif_residue_dir: directory with representative motif residue PDBs
        :param clashing_cutoff: distance cutoff for clashing residues in angstroms. Default set to 2.
        :return:
        """
        # Generate and output a yaml file for each conformer and list which residues clash
        # This way I can output all representative motif residues once and just ignore the ones that clash for each conformer
        # YAML will be a dict with conformer pdbs as keys and list of residue indices as values

        residue_ligand_clash_dict = {}

        target_molecule = os.path.normpath(self.user_defined_dir)
        for conformer in directory_check(self.residue_ligand_interactions_dir):
            clashing_residue_indices = []
            for motif_residue in pdb_check(conformer):
                interaction_prody = prody.parsePDB(motif_residue)
                residue_prody = interaction_prody.select('protein')
                ligand_prody = interaction_prody.select('hetero')

                if minimum_contact_distance(residue_prody, ligand_prody) <= cutoff_distance:
                    residue_index = re.split('-|\.', os.path.basename(motif_residue))[1]
                    clashing_residue_indices.append(residue_index)

            residue_ligand_clash_dict[os.path.basename(conformer)] = clashing_residue_indices

        yaml.dump(residue_ligand_clash_dict,
                  open(os.path.join(target_molecule, 'Inputs', 'User_Inputs', 'Residue_Ligand_Clash_List.yml'), 'w'))

    # DEPRECIATED
    def generate_motif_residues(self):
        """
        Define motif residues from clusters outlined in a user defined yaml file for a single ligand conformation.
        A yaml file will also be generated to keep track of which residues clash with each conformer. This way I can
        output all representative motif residues once and just ignore the ones that clash for each conformer
        * Inputs/User_Inputs/Motif_Clusters.yaml
        """
        # Assemble dict to store list of prody residue instances for each cluster
        fragment_prody_dict = self.assemble_cluster_dict(names_only=False)

        # Okay, cool. Now that I have all of my clusters neatly organized in a dict, I can go through and generate motif
        # residues from each of the clusters based on residue type

        fragment_output_dir = os.path.join(self.user_defined_dir, 'Motifs', 'Representative_Residue_Motifs')
        os.makedirs(fragment_output_dir, exist_ok=True)

        # Start going through fragments and clusters
        total_motif_count = 1
        for fragment in fragment_prody_dict:
            motif_count = 1

            for cluster in fragment_prody_dict[fragment]:
                # For each type of residue in cluster
                residue_types = sorted(list(set([residue.getResnames()[0] for residue in fragment_prody_dict[fragment][cluster]])))

                for res in residue_types:
                    # Calculate average residue coordinates
                    cluster_residues = [residue for residue in fragment_prody_dict[fragment][cluster]
                                        if all([residue.getResnames()[0] == res, len(residue) == self.expected_atoms[residue.getResnames()[0]]])]

                    if len(cluster_residues) > 0:
                        representative_residue = self._select_respresentative_residue(cluster_residues)
                        residue_type = representative_residue.getResnames()[0]

                        # Output residue
                        prody.writePDB(os.path.join(fragment_output_dir, '{4}-{5}-{0}-Cluster_{1}-Motif_{2}-{3}'.format(fragment, cluster, motif_count, residue_type, total_motif_count, len(fragment_prody_dict[fragment][cluster]))),
                                       representative_residue)
                        total_motif_count += 1

                motif_count += 1


class Generate_Binding_Sites():
    """
    This is a class for combining specified groups of binding motifs into hypothetical binding sites
    User input required to deal with residues that obviously do not get along with other motif residues or the ligand
    e.g. steric clashing, non-favorable interactions
    """

    ref2015_weights = {'fa_atr': 1,
                       'fa_elec': 1,
                       'hbond_sc': 1,
                       'hbond_bb_sc': 1,
                       'fa_rep': 0.55
                       }

    def __init__(self, user_defined_dir, residue_groups=None, hypothetical_binding_sites=None, config_dict=None):
        self.user_defined_dir = user_defined_dir
        self.residue_groups = residue_groups
        self.hypothetical_binding_sites = hypothetical_binding_sites

        self.constraints_path = os.path.join(self.user_defined_dir, 'Complete_Matcher_Constraints')
        self.binding_site_pdbs = os.path.join(self.constraints_path, 'Binding_Site_PDBs')
        self.complete_constraint_files = os.path.join(self.constraints_path, 'Constraint_Files')

        self.user_config = config_dict

    def calculate_energies_and_rank(self, use_mysql=True, motif_size=4):
        """
        Calculate total binding site interaction energies for all possible conformations of representative binding 
        motifs and rank.
        :return: 
        """
        print('\nGenerating binding motifs of size {}...\n'.format(motif_size))

        # Retrieve list of Residue-Residue and Residue-Ligand Clashes
        residue_ligand_clash_dict = yaml.load(open(os.path.join(self.user_defined_dir, 'Inputs', 'User_Inputs', 'Residue_Ligand_Clash_List.yml'), 'r'))
        residue_residue_clash_dict = yaml.load(open(os.path.join(self.user_defined_dir, 'Inputs', 'User_Inputs', 'Residue_Residue_Clash_COO.yml'), 'r'))

        # Import Rosetta Score df
        current_ligand = os.path.basename(os.path.normpath(self.user_defined_dir))
        score_df = pd.read_csv(os.path.join(self.user_defined_dir, 'Motifs', 'Residue_Ligand_Interactions', '{}_scores_df.csv'.format(current_ligand)), index_col='description')

        # Set up agg score dict (in the lamest fashion possible)
        score_agg_dict = {}
        for conformer in pdb_check(os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs'), base_only=True, conformer_check=True):
            score_agg_dict[conformer.split('.')[0]] = {}

        for index, row in score_df.iterrows():
            current_conformer = index.split('-')[0]
            motif_res_index = re.split('-|_', index)[2]

            score_agg_dict[current_conformer][motif_res_index] = sum([row['fa_atr'] * self.ref2015_weights['fa_atr'],
                                                                      row['fa_elec'] * self.ref2015_weights['fa_elec'],
                                                                      row['hbond_bb_sc'] * self.ref2015_weights['hbond_bb_sc'],
                                                                      row['hbond_sc'] * self.ref2015_weights['hbond_sc'],
                                                                      row['fa_rep'] * self.ref2015_weights['fa_rep']
                                                                      ])

        if use_mysql:
              self._use_mysql(residue_residue_clash_dict, residue_ligand_clash_dict, score_agg_dict, motif_size)
        else:
              self._use_sqlite(residue_residue_clash_dict, residue_ligand_clash_dict, score_agg_dict)

    def _generate_sqlite_db(self):
        """
        Generate SQLite database + table for storing binding site scores
        :return: 
        """
        db_path = os.path.join(self.user_defined_dir, 'Inputs', 'Scored_Binding_Motifs-sqlite.db')
        db_existed = os.path.exists(db_path)
        connection = sqlite3.connect(db_path)
        cur = connection.cursor()

        if not db_existed:
            cur.execute("create table binding_motif_scores (conformer TEXT NOT NULL, first INTEGER NOT NULL, second INTEGER NOT NULL, third INTEGER NOT NULL, fourth INTEGER NOT NULL, score REAL NOT NULL, PRIMARY KEY (conformer, first, second, third, fourth))")

        return connection, cur

    def _use_sqlite(self, residue_residue_clash_dict, residue_ligand_clash_dict, score_agg_dict):
        # Generate SQLite DB
        sqlite_connection, sqlite_cusor = self._generate_sqlite_db()

        # List of residue indicies
        rep_motif_path = os.path.join(self.user_defined_dir, 'Representative_Residue_Motifs')
        representative_motif_residue_indices = [motif.split('-')[0] for motif in
                                                pdb_check(rep_motif_path, base_only=True)]

        # For every conformer I've generated
        for conformer, residue_list in residue_residue_clash_dict.items():

            # Get rid of any residues that clash with the ligand straight away
            useable_residues_distance = list(
                set(representative_motif_residue_indices) - set(residue_ligand_clash_dict[conformer + '.pdb']))
            useable_residues = list(set(useable_residues_distance) - set(
                [key for key, value in score_agg_dict[conformer].items() if value > 0]))

            # Generate all XXX choose 4 binding site configurations
            list_of_residue_combinations = itertools.combinations(useable_residues, 4)

            # For each combination of four residues...
            for combo in list_of_residue_combinations:
                # For each pair or residues, check if tuple is in residue_residue_clash_dict
                fail_residue_combo = any([(pair in residue_list) for pair in itertools.combinations(combo, 2)])

                # If there are no clashes, calculate Rosetta score (sum(fa_atr, fa_elec, hbond_sc))
                if not fail_residue_combo:
                    total_score = sum([score_agg_dict[conformer][res] for res in combo])

                    # Push to SQLite DB table if score is < 0
                    if total_score < 0:
                        sorted_combo = sorted(combo)
                        # print('Inserting {:10} {:3} {:3} {:3} {:3} {:7.3f} into SQLite3 DB'.format(conformer, sorted_combo[0], sorted_combo[1], sorted_combo[2], sorted_combo[3], total_score))
                        sqlite_cusor.execute(
                            "INSERT OR IGNORE INTO binding_motif_scores (conformer, first, second, third, fourth, score) VALUES (?,?,?,?,?,?)",
                            (str(conformer), int(sorted_combo[0]), int(sorted_combo[1]), int(sorted_combo[2]),
                             int(sorted_combo[3]), float(total_score)))

            sqlite_connection.commit()

        table = pd.read_sql_query("SELECT * from binding_motif_scores", sqlite_connection)
        table.to_csv("ASDF" + '.csv', index_label='index')
        sqlite_connection.close()

    def _generate_mysql_db(self):
        """
        Generate a MySQL database + table for storing binding site scores
        Switching to MySQL to support multiprocessing
        :return: 
        """
        # Requires initial setup of MySQL database
        # Create ~/.my.cnf with [client] header (username, password)
        # Ensure MySQL permissions are set so that user can read, write, and execute
        connection = MySQLdb.connect(host='localhost', read_default_file="~/.my.cnf")
        connection.query('CREATE DATABASE IF NOT EXISTS scored_binding_motifs_{}'.format(os.path.basename(os.path.normpath(self.user_defined_dir))))
        connection.query('USE scored_binding_motifs_{}'.format(os.path.basename(os.path.normpath(self.user_defined_dir))))
        # connection.query('CREATE TABLE IF NOT EXISTS binding_motif_scores (conformer VARCHAR(10) NOT NULL, first INTEGER NOT NULL, second INTEGER NOT NULL, third INTEGER NOT NULL, fourth INTEGER NOT NULL, score REAL NOT NULL, PRIMARY KEY (conformer, first, second, third, fourth))')
        connection.query(
            '''
            CREATE TABLE IF NOT EXISTS binding_motif_score(
            ID int NOT NULL AUTO_INCREMENT,
            conformer VARCHAR(10) NOT NULL,
            score REAL NOT NULL,
            PRIMARY KEY(ID, conformer, score)
            )
            '''
        )
        connection.query(
            '''
            CREATE TABLE IF NOT EXISTS binding_motif_residue(
            ID int NOT NULL,
            residue_id int NOT NULL,
            PRIMARY KEY(ID, residue_id),
            CONSTRAINT FOREIGN KEY (ID) REFERENCES binding_motif_score (ID)
            )
            '''
        )
        cursor = connection.cursor()
        return connection, cursor

    def _use_mysql(self, residue_residue_clash_dict, residue_ligand_clash_dict, score_agg_dict, motif_size=4):
        """
        Evaluates viable binding sites for each conformer and pushes the binding motif + score to a database
        Multiprocessing supported for each unique conformer
        
        :param residue_residue_clash_dict: 
        :param residue_ligand_clash_dict: 
        :param score_agg_dict: 
        :return: 
        """
        # Generate SQLite DB
        mysql_connection, mysql_cursor = self._generate_mysql_db()

        # List of residue indicies
        rep_motif_path = os.path.join(self.user_defined_dir, 'Motifs', 'Representative_Residue_Motifs')
        representative_motif_residue_indices = [motif.split('-')[0] for motif in
                                                pdb_check(rep_motif_path, base_only=True)]

        def push_scores_to_db(residue_residue_clash_dict_tuple, filter=True, filter_percentage=0.8, score_cutoff=0.4):
            """
            Push cumulative weighted residue-ligand scores scores to MySQL database
            :param residue_residue_clash_dict_tuple: 
            :param filter: If true, only use motif residues with best weighted residue-ligand interaction energies
            :param filter_percentage: Upper percentage of motif residues to use
            :param score_cutoff: Only commit combinations with scores in the top score_cutoff percentile
            :return: 
            """
            connection_embed = MySQLdb.connect(host='localhost',
                                               db='scored_binding_motifs_{}'.format(os.path.basename(os.path.normpath(self.user_defined_dir))),
                                               read_default_file="~/.my.cnf")
            cursor_embed = connection_embed.cursor()

            # For every conformer I've generated...
            conformer = residue_residue_clash_dict_tuple[0]
            residue_list = residue_residue_clash_dict_tuple[1]

            # Get rid of any residues that clash with the ligand straight away
            useable_residues_distance = list(set(representative_motif_residue_indices) - set(residue_ligand_clash_dict[conformer]))

            if filter:
                # we only really care about residues that make good interactions with the ligand
                import operator
                residues_sorted_by_score = sorted(score_agg_dict[conformer].items(), key=operator.itemgetter(1))

                # Drop scores > 0
                residues_filtered_zero = [residue for residue in residues_sorted_by_score if residue[1] < 0]

                # Drop remaining residues based on filter_percentage
                best_scoring_residues = int(len(residues_filtered_zero) * filter_percentage)
                residues_filtered_top_scoring = residues_filtered_zero[:best_scoring_residues]
                useable_residues = list(set(useable_residues_distance) & set([residue[0] for residue in residues_filtered_top_scoring]))

            else:
                useable_residues = list(set(useable_residues_distance) - set([key for key, value in score_agg_dict[conformer].items() if value > 0]))

            # Generate all XXX choose X binding site configurations
            list_of_residue_combinations = itertools.combinations(useable_residues, motif_size)

            # Calculate number of rows to be calculated
            number_of_jobs = comb(len(useable_residues), motif_size, exact=True)
            print('\nWorking on {}: {} unique hypothetical binding motifs to be evaluated'.format(residue_residue_clash_dict_tuple[0], number_of_jobs))

            # Count completed binding motif combinations
            overall_count = 0
            committed_count = 0

            # Track top score and only keep scores within, say, 0.2 of the best. This should cut down on required
            # commits and runtime by reducing the number of required commits to DB. Commits are the limiting factor
            # since the DB is locked everytime I call LAST_INSERT_ID()
            best_score = 0

            for combo in list_of_residue_combinations:

                total_score = sum([score_agg_dict[conformer][res] for res in combo])

                if total_score < best_score:
                    best_score = total_score

                if total_score < best_score * score_cutoff:
                    cursor_embed.execute("""INSERT IGNORE INTO binding_motif_score (conformer, score) VALUES (%s,%s)""", (conformer, total_score,))
                    for residue_id in combo:
                        cursor_embed.execute("""INSERT IGNORE INTO binding_motif_residue (ID, residue_id) VALUES (LAST_INSERT_ID(), %s)""", (residue_id,))

                    connection_embed.commit()
                    committed_count += 1

                    # Trying small transactions to speed things up...
                    if committed_count % 100 == 0:
                        connection_embed.commit()

                overall_count += 1

                if overall_count % 10000 == 0:
                    print('Evaluated {} of {} for {}!'.format(overall_count, number_of_jobs, conformer))

            print('{} completed!!!'.format(conformer))

        process = Pool()
        process.map(push_scores_to_db, residue_residue_clash_dict.items())
        process.close()
        process.join()

        mysql_connection.close()

    def generate_binding_site_constraints(self, score_cutoff=-10, secondary_matching=False, use_mysql=True):
        """
        Generate constraint files (and optionally binding site PDBs) for binding sites that pass score filters
        :return: 
        """
        if use_mysql:
            mysql_connection = MySQLdb.connect(host='localhost',
                                               db='scored_binding_motifs_{}'.format(os.path.basename(os.path.normpath(self.user_defined_dir))),
                                               read_default_file="~/.my.cnf")
            mysql_cursor = mysql_connection.cursor()
            mysql_cursor.execute(
                """
                SELECT binding_motif_score.ID, binding_motif_score.conformer, binding_motif_residue.residue_id, 
                binding_motif_score.score FROM binding_motif_score INNER JOIN binding_motif_residue on 
                binding_motif_score.ID = binding_motif_residue.ID WHERE score < {}
                """.format(score_cutoff)
            )
            score_table_rows = pd.DataFrame([{'ID': row[0],
                                              'conformer': row[1],
                                              'residue_id': row[2],
                                              'score': row[3]
                                              } for row in mysql_cursor.fetchall()])

            motifs_to_consider = len(score_table_rows['ID'].unique())
            print('\nConsidering {} binding motifs...\n'.format(motifs_to_consider))

        # SQLITE3 is very out of date....
        else:
            sqlite_connection, sqlite_cusor = self._generate_sqlite_db()
            sqlite_cusor.execute("SELECT * FROM binding_motif_score WHERE score < {}".format(score_cutoff))
            score_table_rows = sqlite_cusor.fetchall()

        # Create directories and files
        os.makedirs(self.constraints_path, exist_ok=True)
        os.makedirs(self.binding_site_pdbs, exist_ok=True)
        os.makedirs(self.complete_constraint_files, exist_ok=True)

        # Touch
        open(os.path.join(self.binding_site_pdbs, 'Binding_sites_to_score.txt'), 'w').close()

        if len(score_table_rows) > 0:

            # Import residue-residue clash dict
            residue_residue_clash_dict = yaml.load(open(os.path.join(self.user_defined_dir, 'Inputs', 'User_Inputs', 'Residue_Residue_Clash_COO.yml'),'r'))

            # Count number of viable binding motif combinations
            motif_count = 0
            motifs_considered = 0

            for index, motif in score_table_rows.groupby(['ID']):
                row_conformer = motif['conformer'].values[0]
                row_motif_indicies = motif['residue_id'].values
                row_score = motif['score'].values[0]

                # Check that none of the residues clash
                residue_list = residue_residue_clash_dict[row_conformer]
                if all([(pair not in residue_list) for pair in itertools.combinations(row_motif_indicies, 2)]):

                    # Generate binding site PDB
                    self._generate_constraint_file_binding_site(row_conformer, row_motif_indicies)

                    # Get individual motif residue constraint blocks
                    motif_constraint_block_list = [open(os.path.join(self.user_defined_dir, 'Motifs', 'Single_Constraints', '{}-{}.cst'.format(row_conformer, index))).read() for index in row_motif_indicies]

                    # Write constraint blocks to file
                    with open(os.path.join(self.complete_constraint_files, '{}-{}.cst'.format(row_conformer, '-'.join([str(a) for a in row_motif_indicies]))), 'w') as complete_constraint_file:

                        # Options for secondary matching
                        for index, block in enumerate(motif_constraint_block_list):
                            complete_constraint_file.write('CST::BEGIN\n')
                            complete_constraint_file.write(motif_constraint_block_list[index])
                            complete_constraint_file.write('\n')
                            if secondary_matching and index != 0:
                                complete_constraint_file.write('  ALGORITHM_INFO:: match\n')
                                complete_constraint_file.write('    SECONDARY_MATCH: DOWNSTREAM\n')
                                complete_constraint_file.write('  ALGORITHM_INFO::END\n')
                            complete_constraint_file.write('CST::END\n')

                    motif_count += 1

                motifs_considered += 1

                if motifs_considered % 1000 == 0:
                    print('{}/{} motifs considered...'.format(motifs_considered, motifs_to_consider))

            # Report number of viable binding motifs
            print('\nFound {} viable binding motifs, proceeding to scoring...\n'.format(motif_count))

            # Score complete binding sites
            self._score_constraint_file_binding_site()
            print('Done!')

        else:
            print('There are no motifs with a score less than {}!'.format(score_cutoff))

    def _generate_constraint_file_binding_site(self, conformer, motif_indicies):
        """
        Generates a binding site PDB based on a constraint file
        :param conformer: 
        :param motif_indicies: 
        :return: 
        """
        #todo: update so that conformers constraints are supported
        conformer_path = os.path.join(self.user_defined_dir, 'Motifs', 'Residue_Ligand_Interactions', conformer)

        # Use residue-ligand pairs in Inputs/Residue_Ligand_Interactions
        relevant_residue_list = [open(os.path.join(conformer_path, '{}-{}.pdb'.format(conformer, str(index)))) for index in motif_indicies]
        conformer_pdb = open(os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs', '{}.pdb'.format(conformer))).read()

        #  Ligand HETATM + Residue ATOM records
        binding_site_pdb_path = os.path.join(self.binding_site_pdbs, '{}-{}.pdb'.format(conformer, '-'.join([str(a) for a in motif_indicies])))
        with open(binding_site_pdb_path, 'w') as binding_site_pdb:
            binding_site_pdb.write('# {} {}\n'.format(conformer, ' '.join([str(a) for a in motif_indicies])))
            for res in relevant_residue_list:
                for line in res:
                    if line.split()[0] == 'ATOM':
                        binding_site_pdb.write(line)
            binding_site_pdb.write(conformer_pdb)

        # Add path to .txt for scoring later
        with open(os.path.join(self.binding_site_pdbs, 'Binding_sites_to_score.txt'), 'a') as score_list:
            score_list.write('{}\n'.format(binding_site_pdb_path))

    def _score_constraint_file_binding_site(self):
        """
        Scores a binding site as defined by a constraint file generated using calculate_energies_and_rank()
        This is to screen for residue-residue clashes that were not caught with the distance filter
        :return: 
        """

        # Calculate scores for all PDBs in this new directory (batch process, don't forget .params)
        current_ligand = os.path.basename(os.path.normpath(self.user_defined_dir))
        scores_dir = os.path.join(self.constraints_path, '{}_scores_raw.txt'.format(current_ligand))
        score_sites_list_path = os.path.join(self.binding_site_pdbs, 'Binding_sites_to_score.txt')

        # todo: update paths with config inputs
        run_jd2_score = subprocess.Popen([os.path.join(self.user_config['Rosetta_path'], '/main/source/bin/score_jd2.{}'.format(self.user_config['Rosetta_compiler'])),
                                          '-l',
                                          score_sites_list_path,
                                          '-extra_res_fa',
                                          os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs', '{}.params'.format(current_ligand)),
                                          '-out:file:scorefile',
                                          scores_dir,
                                          '-out:file:score_only',
                                          '-overwrite'
                                          ])
        run_jd2_score.wait()

        # Process final score file with weighted scores
        print('Raw score file outputted as {}!'.format(scores_dir))

        score_df = pd.read_csv(scores_dir,
                               delim_whitespace=True,
                               index_col=5,
                               usecols=['hbond_sc',
                                        'hbond_bb_sc',
                                        'fa_elec',
                                        'fa_atr',
                                        'fa_rep',
                                        'description'],
                               skiprows=0,
                               header=1)

        for index, row in score_df.iterrows():
            score_df.loc[index:, 'total_score'] = sum([row['fa_atr'] * self.ref2015_weights['fa_atr'],
                                                       row['fa_elec'] * self.ref2015_weights['fa_elec'],
                                                       row['hbond_sc'] * self.ref2015_weights['hbond_sc'],
                                                       row['hbond_bb_sc'] * self.ref2015_weights['hbond_bb_sc'],
                                                       row['fa_rep'] * self.ref2015_weights['fa_rep']
                                                       ])

        # Output only necessary scores to a .csv
        score_df.to_csv(os.path.join(self.constraints_path, '{}_scores_df.csv'.format(current_ligand)))

    # DEPRECIATED
    def generate_binding_sites_by_hand(self):
        """
        This method takes the user defined residue groups and combines them as specified in the hypothetical_binding_sites
        yaml file. This method will output all hypothetical binding sites into a new directory and rank them based on the
        sum total of cluster members from which the representative residues were derived.
        
        {Rank}-{res_#}_{res_#}_{res_#}_{res_#}-{Description}.pdb
        
        :return: 
        """
        # For each hypothetical binding site defined in the yaml file
        for hbs in self.hypothetical_binding_sites:
            # Generate output path
            output_path = os.path.join(self.user_defined_dir, 'Hypothetical_Binding_Sites', hbs)
            os.makedirs(output_path, exist_ok=True)

            # Generate all possible residue combinations for all positions defined
            list_of_residue_combinations = itertools.product(*[self.residue_groups[group] for group in self.hypothetical_binding_sites[hbs]])

            # For each unique binding site
            for residue_combination in list_of_residue_combinations:

                binding_residue_list = []
                # Open all residue pdbs from representative binding residues directory as defined in list
                # Calculate sum of all cluster representitives
                cluster_sum = 0

                for residue in residue_combination:
                    for pdb in pdb_check(os.path.join(self.user_defined_dir, 'Motifs', 'Representative_Residue_Motifs'), base_only=True):
                        if pdb.split('-')[0] == str(residue):
                            binding_residue_list.append(prody.parsePDB(os.path.join(self.user_defined_dir, 'Motifs', 'Representative_Residue_Motifs', pdb)))
                            cluster_sum += int(pdb.split('-')[1])

                # Simple distance check for all residues to prevent clashing
                clashing = False
                cutoff_distance = 2

                for outer_index, a in enumerate(binding_residue_list):
                    for inner_index, b in enumerate(binding_residue_list[(outer_index+1):]):
                        if minimum_contact_distance(a, b) < cutoff_distance:
                            clashing = True
                            print('CLASHING!!!!~!!!!!!!!')
                            break
                    if clashing:
                        break

                # Proceed to generating and outputting binding site if non-clashing
                if not clashing:
                    # Combine all residues into the same file, including the ligand
                    complete_binding_site = prody.parsePDB(os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs', '{}_0001.pdb'.format(os.path.normpath(self.user_defined_dir))))
                    for binding_residue in binding_residue_list:
                        complete_binding_site = complete_binding_site + binding_residue

                    # Generate Constraints
                    binding_site_description = '{0}-{1}.pdb'.format(cluster_sum, '_'.join([str(a) for a in residue_combination]))
                    self.generate_binding_site_constraints(binding_site_description, complete_binding_site)

                    # Output PDB
                    prody.writePDB(os.path.join(output_path, binding_site_description), complete_binding_site)

    def generate_binding_sites_from_gurobi(self):
        """
        20170829 - So I was thinking about it... once I generate a single pose with all of my motif residues, I can pull
        residues directly from it to generate my binding motif PDBs. No need to keep track of cluster residue names and 
        all that except for traceback purposes.
        
        I keep track of residue indicies in the single poses. Once I get results from Gurobi, I can just use prody to
        select those residues from the single pose and output them... will also need to come up with a way to generate
        the constraint files in parallel... should be easy enough...
        :return: 
        """
        pass