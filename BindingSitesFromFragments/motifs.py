#!/usr/bin/env python3

import io
import re
import os
import copy
import json
import pickle
import shutil
import collections
from pprint import pprint
from ast import literal_eval

import yaml
import prody
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import rdFMCS

import pyrosetta
from pyrosetta import rosetta

from .utils import *

# --- Silence ProDy --- #
prody.confProDy(verbosity='none')

# todo: there is no reason for this to be a class anymore...
# todo: move Constraints to its own file
class Generate_Constraints(object):

    three_to_one = {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K',
                    'LEU':'L', 'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V',
                    'TRP':'W', 'TYR':'Y'}

    def __init__(self, chemical_component_identifier):
        self.chemical_component_identifier = chemical_component_identifier

    def determine_constraint_atoms(self, residue_prody_H, conformer_prody_H, fragment=None, verbose=False):
        """
        Determine the three atoms in the ligand and the residue that will be used to calculate ideal values for the DOFs
        required for matcher.

        :param residue_prody: prody object of residue
        :param conformer_prody: prody object of current conformer
        :param fragment: Path to PDB representation of current fragment (e.g. "Inputs/Fragment_Inputs/Fragment_1.pdb")
        :return: 
        """

        residue_prody = residue_prody_H.select('not hydrogen').copy()
        conformer_prody = conformer_prody_H.select('not hydrogen').copy()

        ligand_code = conformer_prody.getResnames()[0]

        residue_index_atom_map = {atom.getIndex(): atom.getName() for atom in residue_prody}
        residue_atom_index_map = {v: k for k, v in residue_index_atom_map.items()}
        ligand_index_atom_map = {atom.getIndex(): atom.getName() for atom in conformer_prody}

        # Load residue and ligand as RDKit Mol objects
        RD_residue = RDKit_Mol_from_ProDy(residue_prody, removeHs=True)
        RD_ligand = RDKit_Mol_from_ProDy(conformer_prody, removeHs=True)

        # Calculate closest atom-atom contacts and two additional atoms for determining bond torsions and angles
        # NOTE: Contact distance and indicies are for residue and ligand with hydrogens stripped!

        # Select source fragment atoms from current conformer is fragment is provided, else use entire ligand
        if fragment:
            # fragment_atoms_prody = prody.parsePDB(os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs', f'Fragment_{fragment}.pdb'))
            fragment_atoms_prody = prody.parsePDB(fragment).select('not hydrogen').copy()
            fragment_atom_names = fragment_atoms_prody.getNames()
            ligand_atoms_for_determining_contacts = conformer_prody.select(f'resname {ligand_code} and name {" ".join(fragment_atom_names)}')

        else:
            ligand_atoms_for_determining_contacts = conformer_prody

        contact_distance, residue_index_low, ligand_index_low = minimum_contact_distance(residue_prody, ligand_atoms_for_determining_contacts, return_indices=True)

        residue_contact_atom = residue_prody[residue_index_low]
        ligand_contact_atomname = ligand_atoms_for_determining_contacts[ligand_index_low].getName()
        ligand_contact_atom = conformer_prody.select(f'name {ligand_contact_atomname}')[0]


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

        # I literally don't even know, I think it's catching deuteriums...
        elif residue_contact_atom.getName()[0] not in ['C', 'N', 'O', 'S', 'P']:
            print('wtf...')
            prody.writePDB(f'{residue_prody_H.getData("contact_source")[0]}.pdb', residue_prody_H)
            return False

        else:
            residue_second_atom = self._determine_next_residue_constraint_atom(residue_contact_atom.getIndex(), RD_residue, residue_prody)
            if residue_second_atom == None: return False

            residue_third_atom = self._determine_next_residue_constraint_atom(residue_second_atom, RD_residue, residue_prody)
            if residue_third_atom == None: return False

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
        residue_index_atom_map = {atom.getIndex(): atom.getName() for atom in residue_prody}
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

            return None
            # raise Exception('Somehow you\'ve run out of neighbors for your contact atom...')

    def _determine_ligand_constraint_atoms(self, current_atom_index, RD_ligand, ligand_fragment_prody):
        """
        Grabs the next ligand atoms down the line for ligand constraint records
        :param current_atom_index: index of ligand contact atom in RD_ligand
        :param RD_ligand: RDKit molecule of FULL ligand
        :param ligand_fragment_prody: prody of ligand FRAGMENT from current conformer
        :return: 
        """
        ligand_index_atom_map = {atom.getIndex(): atom.getName() for atom in ligand_fragment_prody.select('not hydrogen')}
        ligand_contact_atom_neighbors = [atom for atom in RD_ligand.GetAtomWithIdx(int(current_atom_index)).GetNeighbors() if atom.GetIdx() in ligand_index_atom_map.keys()]

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

    def calculate_matcher_dofs(self, residue_prody, ligand_prody, fragment=None, return_constraint_dict=True):
        """
        Calcualtes DOFs used by the matcher to build a ligand out from a residue
        :return:
        """

        constraint_atoms_dict = self.determine_constraint_atoms(residue_prody, ligand_prody, fragment=fragment)

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

        dof_tuple = ideal_distance, ideal_angle_A, ideal_angle_B, ideal_torsion_A, ideal_torsion_AB, ideal_torsion_B

        if return_constraint_dict:
            return dof_tuple, constraint_atoms_dict
        else:
            return dof_tuple

    def generate_single_constraint_block(self, fragment_dir, single_pose_prody, residue_index, fragment=None,
                                         distance_tolerance=0.5, angle_A_tolerance=5, angle_B_tolerance=5,
                                         torsion_A_tolerance=5, torsion_AB_tolerance=5, torsion_B_tolerance=5,
                                         torsion_constraint_sample_number=1, angle_constraint_sample_number=1,
                                         distance_constraint_sample_number=0, greasy_sampling=False):
        """
        Wrapper to generate a single constraint block from a fuzzball interaction
        :param fragment_dir: Path to directory containing fragment PDBs (Fragment_Inputs)
        :param single_pose_prody:
        :param residue_index:
        :return:
        """

        # single_pose_prody should be loaded from .ag.npz
        ligand_prody = single_pose_prody.select('resname {}'.format(self.chemical_component_identifier)).copy()
        residue_prody = single_pose_prody.select('resnum {} and not hydrogen'.format(residue_index)).copy()

        # Get fragment path
        residue_fragment_data = residue_prody.getData("fragment_id")
        current_fragment = None if residue_fragment_data is None else f'Fragment_{residue_fragment_data[0]}'
        current_fragment_path = os.path.join(fragment_dir, f'{current_fragment}.pdb') if current_fragment else None

        return self.generate_single_constraint_block_base(residue_prody, ligand_prody, fragment=current_fragment_path,
                                                          distance_tolerance=distance_tolerance, angle_A_tolerance=angle_A_tolerance, angle_B_tolerance=angle_B_tolerance,
                                                          torsion_A_tolerance=torsion_A_tolerance, torsion_AB_tolerance=torsion_AB_tolerance, torsion_B_tolerance=torsion_B_tolerance,
                                                          torsion_constraint_sample_number=torsion_constraint_sample_number,
                                                          angle_constraint_sample_number=angle_constraint_sample_number,
                                                          distance_constraint_sample_number=distance_constraint_sample_number, greasy_sampling=greasy_sampling)

    def generate_single_constraint_block_base(self, residue_prody, ligand_prody, fragment=None,
                                         distance_tolerance=0.5, angle_A_tolerance=5, angle_B_tolerance=5,
                                         torsion_A_tolerance=5, torsion_AB_tolerance=5, torsion_B_tolerance=5,
                                         torsion_constraint_sample_number=1, angle_constraint_sample_number=1,
                                         distance_constraint_sample_number=0, greasy_sampling=False):
        """
        Generate a single constraint block for one residue-ligand interaction
        :return: 
        """

        dof_tuple, constraint_atoms_dict = self.calculate_matcher_dofs(residue_prody, ligand_prody, fragment=fragment)
        ideal_distance, ideal_angle_A, ideal_angle_B, ideal_torsion_A, ideal_torsion_AB, ideal_torsion_B = dof_tuple

        # Updated this for backbone contacts to allow for any residue identity
        if constraint_atoms_dict['residue']['atom_names'][0] in ['C', 'CA', 'N', 'O', 'OXT']:
            residue_resname = 'A'
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

    def conventional_constraints_from_gurobi_solutions(self, user_defined_dir, gurobi_solutions_csv_dir, fuzzball_dir, constraints_to_generate=100000,
                                                       offset=0, angle_dihedral_tolerance=5, angle_dihedral_sample_number=1,
                                                       iteration=False, greasy_sampling=False, json_output=False, consolidate_solutions=False):
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
                            fragment_inputs_path, conformer_fuzzball, residue_index,
                            angle_A_tolerance=angle_dihedral_tolerance, angle_B_tolerance=angle_dihedral_tolerance,
                            torsion_A_tolerance=angle_dihedral_tolerance,
                            torsion_AB_tolerance=angle_dihedral_tolerance,
                            torsion_B_tolerance=angle_dihedral_tolerance,
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
                            fragment_inputs_path, conformer_fuzzball, residue_index,
                            angle_A_tolerance=angle_dihedral_tolerance, angle_B_tolerance=angle_dihedral_tolerance,
                            torsion_A_tolerance=angle_dihedral_tolerance,
                            torsion_AB_tolerance=angle_dihedral_tolerance,
                            torsion_B_tolerance=angle_dihedral_tolerance,
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
        fragment_inputs_path = os.path.join(user_defined_dir, 'Inputs', 'Fragment_Inputs')
        gurobi_things_path = os.path.join(user_defined_dir, 'Gurobi_Constraints')
        gurobi_constraints_path = os.path.join(gurobi_things_path, 'Constraints')
        gurobi_motif_path = os.path.join(gurobi_things_path, 'Motif_PDBs')

        os.makedirs(gurobi_things_path, exist_ok=True)
        os.makedirs(gurobi_constraints_path, exist_ok=True)
        os.makedirs(gurobi_motif_path, exist_ok=True)

        if consolidate_solutions:
            consolidated_solutions_dir = os.path.join(gurobi_solutions_csv_dir, 'Consolidated_Solutions')
            os.makedirs(consolidated_solutions_dir, exist_ok=True)

        if json_output:
            json_dict = {}

        # --- Get unique fuzzball solution sets --- #

        fuzzball_identifiers = [tuple(file.split('-')[2:4]) for file in os.listdir(gurobi_solutions_csv_dir) if file.endswith('.csv')]
        unique_fuzzball_identifiers = set(fuzzball_identifiers)
        print(unique_fuzzball_identifiers)

        # Get solutions for unique solution set
        for unique_id in unique_fuzzball_identifiers:
            solution_set_list = [file for file in os.listdir(gurobi_solutions_csv_dir) if tuple(file.split('-')[2:4]) == unique_id and file.endswith('.csv')]

            print(unique_id)
            print(solution_set_list)

            # Import all solutions from csv
            gurobi_solutions = pd.DataFrame(columns=['Obj_score', 'Residue_indicies', 'Conformer'])

            for solution_set in solution_set_list:
                temp_solution_df = pd.read_csv(os.path.join(gurobi_solutions_csv_dir, solution_set), usecols=['Obj_score', 'Residue_indicies', 'Conformer'])
                gurobi_solutions = gurobi_solutions.append(temp_solution_df, ignore_index=True)

            # Remove duplicates
            gurobi_solutions = gurobi_solutions.drop_duplicates(subset=['Residue_indicies', 'Conformer'])

            # if not iteration:
            gurobi_solutions = gurobi_solutions.sort_values(by=['Obj_score'], ascending=True).head(n=constraints_to_generate + offset).tail(constraints_to_generate)

            # Group dataframe by conformer
            for current_conformer, conformer_df in gurobi_solutions.groupby(['Conformer']):

                if consolidate_solutions:
                    conformer_df.to_csv(os.path.join(consolidated_solutions_dir, f'{current_conformer}.csv'), index=False)

                # Open pdb only once!
                conformer_fuzzball = prody.loadAtoms(os.path.join(fuzzball_dir, '{}.ag.npz'.format(current_conformer)))

                # Keep track of previously calculated constraint blocks
                previously_calculated_constraint_blocks = {}

                # Generate constraint files or JSON
                _populate_constraint_json() if json_output else _generate_constraint_files()

        if json_output:
            json_output_path = os.path.join(gurobi_constraints_path, '{}-constraint_blocks.json'.format(user_defined_dir))
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
            load transformed residues and current conformer as Pose in PyRosetta (prody > PDBStream > PyRosetta)
            score pose
            save relevant metrics into a dataframe (contact REU, hbond energy, hbond contact atom on ligand)

        decide which residues to use based on metrics from dataframe
        pull residues into a single prody

        write prody to disk as PDB
    """

    resnames = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
                'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']

    def __init__(self, user_defined_dir, fuzzball_dir, config_dict=None):

        # --- Paths --- #
        self.user_defined_dir = user_defined_dir
        self.user_config = config_dict
        self.fuzzball_dir = fuzzball_dir
        self.rosetta_inputs = os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs')
        self.defined_icorenteractions = os.path.join(self.user_defined_dir, 'Inputs', 'Defined_Interactions')

        # --- Other Variables --- #
        self.chemical_component_identifier = os.path.normpath(os.path.basename(self.user_defined_dir))[:3]
        self.reference_conformer_prody = None  # Initialized in _init_load_reference_conformer_and_clusters()
        self.cluster_dict = None  # Initialized in _init_load_reference_conformer_and_clusters()
        self.hbond_atom_set = None  # Initialized in _init_get_ligand_hbond_atoms()

        self.motif_residue_attributes_df = None  # Either loaded in assemble_fuzzball() or score_motif_conformer_interactions()
        self.reference_complex_atomname_map = None  # Assigned if assembling fuzzball for existing complex
        
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
                cluster_name_split = re.split('[_.]', os.path.normpath(os.path.basename(cluster)))
                cluster_cat_index = int(cluster_name_split[1])
                cluster_contact_index = int(cluster_name_split[2])
                if cluster_cat_index not in temp_cluster_dict[fragment].keys():
                    temp_cluster_dict[fragment][cluster_cat_index] = dict()
                temp_cluster_dict[fragment][cluster_cat_index][cluster_contact_index] = prody.loadAtoms(cluster)

        self.cluster_dict = temp_cluster_dict

    def _init_get_ligand_hbond_atoms(self):
        """
        Use RDKit to get the names of ligand atoms that can mediate hydrogen bond interactions.

        :return:
        """
        ligand_mol = RDKit_Mol_from_ProDy(self.reference_conformer_prody)

        # todo: update to use new Rosetta HBond class
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

    def score_motif_conformer_interactions(self, conformer_dir=None, apply_mapping=False):
        """
        Performs the bulk of work for generating fuzzballs

        Conformer directory should contain PDB/params pairs for each ligand conformer
        :param conformer_dir: directory containing ligand conformers to score against
        :param apply_mapping: maps fuzzball atom names onto ligand atoms in an existing complex
        :return:
        """

        # todo: decide whether to save a bunch of dataframes or consolidate... dataframe is currently written only after all residues for all conformers are scored
        # List of dicts for residue information
        fuzzball_list_of_dicts = list()

        # Set conformer directory
        if not conformer_dir:
            conformer_dir = self.rosetta_inputs

        # Make sure conformer directory exists
        if not os.path.exists(conformer_dir):
            raise Exception(f'{conformer_dir} does not exist!')

        # For each conformer
        for conformer in pdb_check(conformer_dir):

            # Initialize PyRosetta with stuff
            conformer_name, _ = os.path.splitext(os.path.basename(os.path.normpath(conformer)))
            # current_params = conformer_name if conformer_dir == self.rosetta_inputs else self.chemical_component_identifier
            # current_params_path = os.path.join(conformer_dir, f'{current_params}.params')
            current_params_path = os.path.join(conformer_dir, f'{conformer_name}.params')

            # core.conformation.Conformation muted, OXT is missing from all motif residues
            my_options = [f"-extra_res_fa {current_params_path} -flip_HNQ -packing:linmem_ig 20",
                          "-mute core.conformation core.chemical"]
            pyrosetta.init(options=' '.join(my_options))

            # Load conformer as Prody object
            conformer_prody = prody.parsePDB(conformer)

            # Load conformer as PyRosetta Pose
            conformer_pose = rosetta.core.import_pose.pose_from_file(conformer)

            # Iterate through self.cluster_dict (for loops all the way down...)
            for fragment in self.cluster_dict:

                fragment_path = os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs', f'{fragment}.pdb')

                # Save transformed clusters in a dict()
                transformed_cluster_dict = dict()

                for cat_cluster in self.cluster_dict[fragment]:
                    for cluster in self.cluster_dict[fragment][cat_cluster]:
                        current_cluster = self.cluster_dict[fragment][cat_cluster][cluster]

                        # Make copy of conformer for current cluster
                        current_conformer_and_transformed_residues = conformer_prody.copy()

                        # --- Transform cluster residues onto target conformer --- #

                        for residue in current_cluster.iterResidues():

                            # todo: condense this
                            # Skip residues if essential backbone atoms are missing
                            defined_residue_bb_atoms = [a.getName() for a in residue.select('backbone')]
                            if len(set(defined_residue_bb_atoms) & {'C', 'CA', 'N', 'O'}) != 4:
                                continue

                            # Skip residues if they possess multiple occupancies (4bn4 ARG 345)
                            if len(set(residue.getAltlocs())) > 1:
                                print(f'{residue.getData("contact_source")[0]} possesses alternate locations and will not be considered.')
                                continue

                            # Skip residue if more than 50% of atoms have zero occupancy (cutoff is completely arbitrary...)
                            if len([a for a in residue.getOccupancies() if a == 0]) > (len(residue) / 2):
                                print(f'More than 50% of atoms in {residue.getData("contact_source")[0]} have 0 occupancy and will not be considered.')
                                continue

                            transformed_motif, constraint_atoms_dict = self.transform_residue_about_current_conformer(conformer_prody, residue, fragment=fragment_path, apply_mapping=apply_mapping)

                            # Add residue to current_conformer_and_transformed_residues collection for scoring
                            if transformed_motif is not None:
                                current_conformer_and_transformed_residues += transformed_motif

                        # Save all transformed cluster residues to transformed_cluster_dict
                        if current_conformer_and_transformed_residues.numResidues() > 1:
                            transformed_cluster_dict[f'{cat_cluster}-{cluster}'] = current_conformer_and_transformed_residues
                        else:
                            print(f'\nAll motif residues from {fragment} Cluster {cat_cluster}-{cluster} were rejected!\n')

                # --- Score all clusters for current fragment using PyRosetta and record scores/attributes --- #

                bsff_sfxn = rosetta.core.scoring.ScoreFunction()
                bsff_sfxn.set_weight(rosetta.core.scoring.fa_atr, 1)
                bsff_sfxn.set_weight(rosetta.core.scoring.fa_rep, 0.55)
                bsff_sfxn.set_weight(rosetta.core.scoring.fa_sol, 1)
                bsff_sfxn.set_weight(rosetta.core.scoring.fa_elec, 1)
                bsff_sfxn.set_weight(rosetta.core.scoring.hbond_sc, 1)
                bsff_sfxn.set_weight(rosetta.core.scoring.hbond_bb_sc, 1)

                sfxn_weights = bsff_sfxn.weights()

                for transformed_cluster in transformed_cluster_dict:

                    transformed_cluster_prody = transformed_cluster_dict[transformed_cluster]

                    transformed_residues_only = transformed_cluster_prody.select(f'not resname {self.chemical_component_identifier}')
                    if transformed_residues_only is None:
                        continue
                    transformed_residues_only_hv = transformed_residues_only.getHierView()
                    ligand_prody = transformed_cluster_prody.select(f'resname {self.chemical_component_identifier} and not hydrogen')

                    # Convert prody of transformed cluster into string
                    # import_pose.pose_from_pdb_stream() uses import_pose.pose_from_pdbstring() anyways and doesn't require BS options
                    transformed_cluster_stream = io.StringIO()
                    prody.writePDBStream(transformed_cluster_stream, transformed_cluster_prody)

                    # Load transformed cluster into PyRosetta
                    transformed_cluster_pose = rosetta.core.pose.Pose()
                    rosetta.core.import_pose.pose_from_pdbstring(transformed_cluster_pose, transformed_cluster_stream.getvalue())
                    bsff_sfxn(transformed_cluster_pose)

                    # Assert only the ligand was removed from selection...
                    # Need to catch residues with zero occupancy... rosetta throws them out while prody is fine

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

                        # OptH and rescore for all polar residues (except ARG) in independent pose
                        if residue_prody.getResname() in ['ASP', 'GLU', 'HIS', 'LYS', 'ASN', 'GLN', 'SER', 'THR', 'TYR']:
                            conformer_residue_opth = conformer_pose.clone()
                            residue_residue = transformed_cluster_pose.residue(index).clone()
                            conformer_residue_opth.append_residue_by_jump(residue_residue, 1)

                            rosetta.core.pack.optimizeH(conformer_residue_opth, bsff_sfxn)
                            bsff_sfxn(conformer_residue_opth)
                            e_edge_opth = conformer_residue_opth.energies().energy_graph()
                            what = e_edge_opth.find_energy_edge(1, 2)

                        else:
                            # Energy-related attributes
                            what = e_edge.find_energy_edge(index, 1)

                        if what is None:
                            continue

                        residue_attribute_dict = dict()

                        # Identification attributes
                        residue_attribute_dict['conformer'] = conformer_name
                        residue_attribute_dict['fragment'] = residue_prody.copy().getData('fragment_id')[0]
                        residue_attribute_dict['cluster'] = residue_prody.copy().getData('cluster')[0]
                        residue_attribute_dict['contact_source'] = residue_prody.copy().getData('contact_source')[0]
                        residue_attribute_dict['resname'] = residue_prody.getResname()

                        current_edge = what.fill_energy_map()
                        total_score = current_edge.dot(sfxn_weights)

                        residue_attribute_dict['total_score'] = total_score

                        hbond_sc_score = what[rosetta.core.scoring.hbond_sc]
                        hbond_bb_sc_score = what[rosetta.core.scoring.hbond_bb_sc]
                        residue_attribute_dict['hbond_sc'] = hbond_sc_score
                        residue_attribute_dict['hbond_bb_sc'] = hbond_bb_sc_score

                        # Report contact atom in self.hbond_atom_set if hbondsc < 0, else None
                        ligand_prody_stripped = ligand_prody.select('not hydrogen')
                        residue_prody_stripped = residue_prody.select('not hydrogen')

                        ligand_residue_distance_matrix = prody.buildDistMatrix(ligand_prody_stripped.getCoords(),
                                                                               residue_prody_stripped.getCoords())

                        if hbond_sc_score < 0 or hbond_bb_sc_score < 0:

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
        temp_df.drop_duplicates(subset=['conformer', 'contact_source', 'hbond_sc', 'hbond_bb_sc', 'ligand_contact', 'resname', 'total_score'], keep='first', inplace=True)
        self.motif_residue_attributes_df = temp_df

        # todo: should really be dumping into SQLite3 database...
        design_dir, conformer_file = os.path.split(conformer_dir)
        dump_dir = design_dir if conformer_dir != self.rosetta_inputs else self.user_defined_dir
        self.motif_residue_attributes_df.to_csv(os.path.join(dump_dir, 'motif_residue_attributes.csv'), index=False)

    def transform_residue_about_current_conformer(self, conformer_prody, residue, fragment=None, reference_conformer=None, apply_mapping=False):
        """
        Transforms a given motif residue about a given conformer to maintain contact geometries with reference fragment
        as observed in PDB. (How to English?)

        :param conformer_prody: prody object of target conformer
        :param residue: prody object of motif residue to transform about conformer_prody
        :param fragment: Path to PDB of current fragment. Uses entire ligand to determine constraint atoms if None.
        :param reference_conformer: ligand conformer to use for determining transformation reference atoms
        :param apply_mapping: maps fuzzball atom names onto ligand atoms in an existing complex
        :return:
        """
        # --- RIPPED FORM GENERATE_MOTIF_RESIDUES.GENERATE_FUZZBALL() --- #

        deepcopy_residue = residue.copy()

        # Check whether current residue can be parsed using RDKit, skip if returns None
        if RDKit_Mol_from_ProDy(deepcopy_residue) is None:
            print(f'RDKit was unable to parse the current residue: {deepcopy_residue.getData("contact_source")[0]}')
            return None, None

        if not reference_conformer:
            reference_conformer = self.reference_conformer_prody

        # Get fragment atoms for transformation of mobile onto reference fragment
        # Use reference conformer to get constraint atoms! Residues were pulled from PDB using fragments in the same
        # coordinate frame as the reference conformer.
        constraint_atoms_dict = self.generate_constraints.determine_constraint_atoms(deepcopy_residue,
                                                                                     reference_conformer,
                                                                                     fragment=fragment, verbose=False)

        # todo: FIX _determine_next_residue_constraint_atoms()!!!!
        if constraint_atoms_dict is False:
            print(f'There was an issue determining residue contact atoms... sorry :|')
            return None, None

        ligand_constraint_atoms = constraint_atoms_dict['ligand']['atom_names']

        # Select atoms from reference ligand conformer (conformer #1 e.g. MEH_0001)
        reference_conformer_coord_list = [reference_conformer.select(f'name {atom}').getCoords()[0] for atom in ligand_constraint_atoms]
        reference_conformer_atoms = np.asarray(reference_conformer_coord_list)

        # Select atoms from current conformer
        if apply_mapping:
            target_conformer_coord_list = [conformer_prody.select(f'name {self.reference_complex_atomname_map[atom]}').getCoords()[0] for atom in ligand_constraint_atoms]
            constraint_atoms_dict['ligand']['atom_names'] = [self.reference_complex_atomname_map[atom] for atom in ligand_constraint_atoms]
        else:
            target_conformer_coord_list = [conformer_prody.select(f'name {atom}').getCoords()[0] for atom in ligand_constraint_atoms]
        target_conformer_atoms = np.asarray(target_conformer_coord_list)

        transformation_matrix = Generate_Constraints.calculate_transformation_matrix(reference_conformer_atoms, target_conformer_atoms)
        return prody.applyTransformation(transformation_matrix, deepcopy_residue), constraint_atoms_dict

    def minimum_fuzzball_definition(self, current_conformer, current_conformer_fuzzball):
        """
        Produces a minimum fuzzball definition with the current conformer and any user-defined residues.

        :param current_conformer: prody of current conformer
        :param current_conformer_fuzzball: current conformer prody to which defined residues will be added
        :return:
        """

        for index, defined_interaction in enumerate(pdb_check(os.path.join(self.user_defined_dir, 'Inputs', 'Defined_Interactions')), start=2):

            print(f'Working on {os.path.basename(os.path.normpath(defined_interaction))}...\n')

            defined_interaction_prody = prody.parsePDB(defined_interaction)
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
            print(f'\n{conformer_with_lowest_rmsd_to_defined} successfully mapped onto defined ligand in {os.path.basename(os.path.normpath(defined_interaction))}!\n')

            # --- Apply transformation to motif residue --- #

            # Residue needs to be in same reference frame as the superposed defined ligand before determining transformation
            superposed_motif = prody.applyTransformation(superpose_transformation, defined_motif_residue)

            motif_residue_transformed, constraint_atoms_dict = self.transform_residue_about_current_conformer(current_conformer, superposed_motif, reference_conformer=defined_ligand_mapped)

            # Tag defined motif residue with input filename
            defined_interaction_filename = os.path.basename(os.path.normpath(defined_interaction)).split('.')[0]
            motif_residue_transformed.setData('contact_source', [defined_interaction_filename] * len(motif_residue_transformed))
            motif_residue_transformed.setData('defined', [1] * len(motif_residue_transformed))

            # Set resnum, occupany, coordset
            motif_residue_transformed.setResnums([index] * len(motif_residue_transformed))
            motif_residue_transformed.setChids(['A'] * len(motif_residue_transformed))
            motif_residue_transformed.setAltlocs([' '] * len(motif_residue_transformed))
            motif_residue_transformed.setOccupancies([1] * len(motif_residue_transformed))

            # Add defined motif to fuzzball
            current_conformer_fuzzball += motif_residue_transformed

        return current_conformer_fuzzball

    def minimum_fuzzball_definition_from_iteration(self, current_conformer_fuzzball, iteration):
        """
        Generate minimum fuzzball definition from a previous fuzzball iteration. I can just copy residues over directly
        since the  motif residues have already been properly transformed.
        :param current_conformer_fuzzball: prody ligand to which defined residues will be added
        :param iteration: prody atomgroup of minimum fuzzball definition pulled from a match
        :return:
        """

        if iteration is None:
            return current_conformer_fuzzball

        iteration_hv = iteration.getHierView()
        for residue in iteration_hv.iterResidues():

            # Tag defined motif residue with input filename
            residue.setData('defined', [1] * len(residue))

            # Set resnum, occupany, coordset
            residue.setResnums([current_conformer_fuzzball.numResidues() + 1] * len(residue))
            residue.setChids(['A'] * len(residue))
            residue.setAltlocs([' '] * len(residue))
            residue.setOccupancies([1] * len(residue))

            current_conformer_fuzzball += residue.copy()

        return current_conformer_fuzzball

    def assemble_fuzzball_for_existing_complex(self, existing_complex_path, ligand_params, ligand_reference,
                                               complex_ligand_id=None, force_limit=False, skip_clean=False):
        """
        Map ligand atom names from a Rosetta params file onto the ligand atoms in an existing protein-ligand complex.
        This is required since the method uses the standardized param file atoms names to keep track of unique
        protein-ligand contacts.

        I think I'm going to write this function as a wrapper for assemble_fuzzball()...

        Fuzzball creation for an existing complex will only happen (as far as I can tell) if I am attempting to design
        around a ligand that has been placed through docking (which we won't be doing) or if I want to benchmark binding
        site recovery for an already existing binding site (which we definitely will be doing). As a result, all of the
        files that I generate for a complex will be limited to that complex and can be stored in a single directory.
        All of the files for performing design will be self contained to this directory when calling
        fuzzball_composition_design() or any complementary RotamerSets benchmark.

        The things I'll need are:
          1. The input protein-ligand complex
          2. A .params file
          3. Fuzzball generated by this wrapper function

        The only complicated bit is the .params file. Fragments used for collecting contact observations were created
        using a reference conformer produced by molfile2params. I'm not sure whether molfile2params is consistent with
        atom name assignments... it's vital for correctly mapping fragments onto the current ligand conformer and
        maintaining contact geometries of observed side-chain interactions. What I'm going to do is copy the .params
        file for the reference conformer (XXX_0001) and use the PDB_CONFORMER flag with an accompanying PDB of the
        complex ligand with the correct atom names mapped.

        :param existing_complex_path: path to PDB of the protein-ligand complex
        :param complex_ligand_id: use this identifier to select the ligand in the provided PDB. Defaults to use
        self.chemical_component_identifier
        :param ligand_params: params file generated for ligand in the existing complex by molfile_to_params.py
        :param ligand_reference: pdb file generated for ligand in existing complex by molfile_to_params.py
        """

        # --- Create paths to dump things --- #
        new_design_name = os.path.basename(os.path.normpath(existing_complex_path)).split('.')[0]
        new_design_path = os.path.join(self.user_defined_dir, 'Design', new_design_name)
        os.makedirs(new_design_path, exist_ok=True)

        conformer_dir = os.path.join(new_design_path, 'Conformers')
        os.makedirs(conformer_dir, exist_ok=True)
        shutil.copy2(ligand_params, conformer_dir)
        shutil.copy2(ligand_reference, conformer_dir)

        # --- Process input complex --- #

        if skip_clean:
            shutil.copy(existing_complex_path, new_design_path)
            clean_pdb_path = existing_complex_path
            clean_prody = prody.parsePDB(clean_pdb_path)

        else:
            existing_complex_prody = prody.parsePDB(existing_complex_path)
            protein_only_prody = existing_complex_prody.select('protein').toAtomGroup()

            # Clean protein-only portion of complex
            clean_prody, removed_residues = clean_pdb(protein_only_prody)

            ligand_selection_string = f'resname {self.chemical_component_identifier if not complex_ligand_id else complex_ligand_id}'
            existing_complex_ligands = existing_complex_prody.select(ligand_selection_string)

            try:
                existing_complex_ligands_hv = existing_complex_ligands.getHierView()
            except FailedAssemblyError:
                raise FailedAssemblyError(f'No ligands with identifier "{complex_ligand_id}" exist in {existing_complex_path}')

            reference_ligand_prody = prody.parsePDB(ligand_reference)
            reference_ligand_mol = RDKit_Mol_from_ProDy(reference_ligand_prody, removeHs=False)

            assert existing_complex_ligands_hv.numResidues() == 1

            # Use RDKit to map idealized ligand onto existing ligand atoms
            for complex_ligand in existing_complex_ligands_hv.iterResidues():

                complex_ligand_mol = RDKit_Mol_from_ProDy(complex_ligand, removeHs=False)
                ligand_mcs = rdFMCS.FindMCS([reference_ligand_mol, complex_ligand_mol], bondCompare=rdFMCS.BondCompare.CompareAny)
                ligand_mcs_mol = Chem.MolFromSmarts(ligand_mcs.smartsString)

                if complex_ligand_mol.GetNumHeavyAtoms() == ligand_mcs_mol.GetNumHeavyAtoms():

                    # Align ligand_ref onto existing complex_ligand
                    ligand_ref_no_H = reference_ligand_prody.select('not hydrogen').copy()
                    ligand_complex_no_H = complex_ligand.select('not hydrogen').copy()
                    # transformation_matrix = prody.calcTransformation(ligand_ref_no_H, ligand_complex_no_H)
                    transformation_matrix = self.generate_constraints.calculate_transformation_matrix(ligand_ref_no_H.getCoords(), ligand_complex_no_H.getCoords())
                    prody.applyTransformation(transformation_matrix, reference_ligand_prody)

                    # Remove old ligand from complex and replace with ligand generated by molfile_to_params.py
                    append_ligand_resnum = clean_prody.numResidues() + 1
                    reference_ligand_prody.setResnums([append_ligand_resnum] * len(reference_ligand_prody))
                    reference_ligand_prody.setChids(['X'] * len(reference_ligand_prody))
                    clean_prody += reference_ligand_prody

                else:
                    raise Exception('Unable to map reference ligand onto complex ligand!')

            # Dump clean scaffold complex
            clean_pdb_path = os.path.join(new_design_path, f'{os.path.basename(os.path.normpath(existing_complex_path)).split(".")[0]}-clean.pdb')
            prody.writePDB(clean_pdb_path, clean_prody)

        # Copy reference params file, write rotamers to file, and append PDB_CONFORMER flag
        # source_params = os.path.join(self.rosetta_inputs, f'{self.chemical_component_identifier}_0001.params')
        # dest_params = os.path.join(conformer_dir, f'{self.chemical_component_identifier}.params')
        # shutil.copy2(source_params, dest_params)

        # --- Perform assembly --- #

        # Make self-contained directory for fuzzball-related stuff
        actual_fuzzball_dir = os.path.join(self.fuzzball_dir, 'Fuzzballs')
        os.makedirs(actual_fuzzball_dir, exist_ok=True)

        for ligand in clean_prody.select(f'resname {self.chemical_component_identifier}').getHierView().iterResidues():
            iteration_dict = {'minimum_definition': None, 'match_path': clean_pdb_path}
            local_ligand_reference = os.path.join(conformer_dir, os.path.basename(ligand_reference))
            fuzzball_path = self.assemble_fuzzball(local_ligand_reference, iteration=iteration_dict, existing=True, existing_params=ligand_params, force_limit=force_limit)

            # Move fuzzball-related stuff into a self-contained directory
            fuzzball_identifier = os.path.basename(fuzzball_path)[:-4]
            fuzzball_identifier_split = fuzzball_identifier.split('-')
            shutil.move(fuzzball_path, actual_fuzzball_dir)
            shutil.move(os.path.join(self.fuzzball_dir, f'{fuzzball_identifier}.ag.npz'), actual_fuzzball_dir)
            shutil.move(os.path.join(self.fuzzball_dir, f'{fuzzball_identifier}.csv'), actual_fuzzball_dir)
            shutil.move(os.path.join(self.fuzzball_dir, f'{fuzzball_identifier_split[1]}-{fuzzball_identifier_split[2]}-match_residue_map.pickle'), actual_fuzzball_dir)

    # todo: fuzzball assembly for an existing complex uses options from iteration...
    def assemble_fuzzball(self, conformer_path, add_user_defined_motifs=False, fuzzball_limit=5000, hbond_limit=500,
                          iteration=None, iteration_name=None, iteration_index=0, fuzzball_index=0, existing=False,
                          existing_params=None, force_limit=False):
        """
        Yiss.

        Output both ag.npz and .pdb representations of fuzzball. PDB representation doesn't save residue source info.

        :param conformer_path: Path to a conformer PDB for which a fuzzball will be assembled
        :param add_user_defined_motifs:
        :param fuzzball_limit:
        :param hbond_limit:
        :param iteration:
        :param iteration_name:
        :param iteration_index:
        :param fuzzball_index:
        :param existing: building a fuzzball for an existing complex
        :param existing_params: path to params file for ligand in the existing complex
        :param force_limit: force iterations to obey motif limits
        :return:
        """

        #####################################################
        # Initialize Conformer and Available Motif Residues #
        #####################################################

        design_dir, current_conformer_path = os.path.split(conformer_path)
        current_conformer_name = current_conformer_path.split('.')[0]

        # Load current conformer and motif_residue_attributes_csv
        current_conformer = prody.parsePDB(conformer_path)
        current_conformer.setData('defined', [0] * len(current_conformer))
        current_conformer_fuzzball = current_conformer.copy()

        # Import previously generated motif_residue_attributes_df, else generate it
        motif_residue_attributes_csv = os.path.join(self.user_defined_dir, 'motif_residue_attributes.csv')
        design_csv = os.path.join(design_dir, '..', 'motif_residue_attributes.csv')

        # If assembling for an existing complex, map complex conformer onto Rosetta_Inputs conformers
        # molfile_to_params does not assign atom names in a consistent manner...
        if existing:
            complex_ligand_mol = RDKit_Mol_from_ProDy(current_conformer, removeHs=False)
            reference_ligand_mol = RDKit_Mol_from_ProDy(self.reference_conformer_prody, removeHs=False)

            ligand_mcs = rdFMCS.FindMCS([reference_ligand_mol, complex_ligand_mol], bondCompare=rdFMCS.BondCompare.CompareAny)
            ligand_mcs_mol = Chem.MolFromSmarts(ligand_mcs.smartsString)

            complex_mol_match = complex_ligand_mol.GetSubstructMatch(ligand_mcs_mol)
            reference_mol_match = reference_ligand_mol.GetSubstructMatch(ligand_mcs_mol)

            reference_complex_atomname_map = dict()
            for complex_idx, reference_idx in zip(complex_mol_match, reference_mol_match):
                complex_atomname = complex_ligand_mol.GetAtomWithIdx(complex_idx).GetMonomerInfo().GetName()
                reference_atomname = reference_ligand_mol.GetAtomWithIdx(reference_idx).GetMonomerInfo().GetName()
                reference_complex_atomname_map[reference_atomname.strip()] = complex_atomname.strip()
                
            self.reference_complex_atomname_map = reference_complex_atomname_map

        # Prepare motif score dataframe
        if existing and os.path.exists(design_csv):
            print('motif_residue_attributes.csv for design found...\n')
            self.motif_residue_attributes_df = pd.read_csv(design_csv)

        elif not existing and os.path.exists(motif_residue_attributes_csv):
            print('motif_residue_attributes.csv found...\n')
            self.motif_residue_attributes_df = pd.read_csv(motif_residue_attributes_csv)
        else:
            print('motif_residue_attributes.csv not found, generating it now...\n')
            conformer_dir = design_dir if existing else None
            self.score_motif_conformer_interactions(conformer_dir=conformer_dir, apply_mapping=existing)

        ###############################################
        # Add user-defined motif residues to fuzzball #
        ###############################################

        print('Processing defined interactions...\n')

        defined_interactions_dir = os.path.join(self.user_defined_dir, 'Inputs', 'Defined_Interactions')
        if iteration:
            current_conformer_fuzzball = self.minimum_fuzzball_definition_from_iteration(current_conformer_fuzzball, iteration['minimum_definition'])
            match_pdb = iteration['match_path']

        elif add_user_defined_motifs and os.path.exists(defined_interactions_dir):
            current_conformer_fuzzball = self.minimum_fuzzball_definition(current_conformer, current_conformer_fuzzball)
            match_pdb = None
        # todo: streamline this
        else:
            match_pdb = None

        # Set params path
        params_path = existing_params if existing else os.path.join(self.rosetta_inputs, f'{current_conformer_name}.params')

        my_options = [f"-extra_res_fa {params_path}",
                      "-mute core.conformation core.chemical core.pack",
                      '-extrachi_cutoff 0']
        pyrosetta.init(options=' '.join(my_options))

        # Import minimum fuzzball prody (ligand + user-defined residues) as pose
        minimum_fuzzball_pose = rosetta.core.pose.Pose()
        minimum_fuzzball_prody_string = io.StringIO()
        prody.writePDBStream(minimum_fuzzball_prody_string, current_conformer_fuzzball)
        rosetta.core.import_pose.pose_from_pdbstring(minimum_fuzzball_pose, minimum_fuzzball_prody_string.getvalue())

        # Number of residues in minimum fuzzball
        residues_in_minimum_fuzzball = current_conformer_fuzzball.numResidues()
        resnum_count = residues_in_minimum_fuzzball

        # sfxn = rosetta.core.scoring.get_score_function()
        # sfxn(minimum_fuzzball_pose)

        # Custom score function using select 2b energies, consider hbond_bb_sc now that binding site is nucleated
        bsff_sfxn = rosetta.core.scoring.ScoreFunction()
        bsff_sfxn.set_weight(rosetta.core.scoring.fa_atr, 1)
        bsff_sfxn.set_weight(rosetta.core.scoring.fa_rep, 0.55)
        bsff_sfxn.set_weight(rosetta.core.scoring.fa_sol, 1)
        bsff_sfxn.set_weight(rosetta.core.scoring.fa_elec, 1)
        bsff_sfxn.set_weight(rosetta.core.scoring.hbond_sc, 1)
        bsff_sfxn(minimum_fuzzball_pose)

        sfxn_weights = bsff_sfxn.weights()

        # Import match_pdb if it exists and align onto fuzzball
        if match_pdb:

            # Rosetta pose of current match
            match_pose = rosetta.core.pose.Pose()
            rosetta.core.import_pose.pose_from_file(match_pose, match_pdb)

            # Get ligand from match, always last residue
            match_pose_size = match_pose.size()
            match_ligand = match_pose.residue(match_pose_size)
            if existing:
                motif_resnums = list()
            else:
                motif_resnums = find_constraint_resnums(match_pdb)

            # Keep track of match positions and compatible residue identites
            match_residue_map = {position: dict() for position in range(1, match_pose.size())}  # Assumes one ligand appended to end of sequence

            # Get ligand from fuzzball, always first residue
            fuzzball_ligand = minimum_fuzzball_pose.residue(1)

            # Calculate rotation/translation by hand using first three atoms of ligand
            mobile_match = rosetta.numeric.xyzTransform_double_t(match_ligand.xyz(1), match_ligand.xyz(2), match_ligand.xyz(3))
            mobile_match_inverse = mobile_match.inverse()
            target_fuzzball = rosetta.numeric.xyzTransform_double_t(fuzzball_ligand.xyz(1), fuzzball_ligand.xyz(2), fuzzball_ligand.xyz(3))

            ligand_rotation = target_fuzzball.R * mobile_match_inverse.R
            ligand_translation = target_fuzzball.R * mobile_match_inverse.t + target_fuzzball.t

            # Apply transformation
            match_pose.apply_transform_Rx_plus_v(ligand_rotation, ligand_translation)

            # --- Generate backbone-independent rotamer library --- #

            residue_conformers = dict()
            canonical_residues_pose = pyrosetta.pose_from_sequence('ACDEFGHIKLMNPQRSTVWY')
            for residue_index in range(1, canonical_residues_pose.size() + 1):
                current_residue = canonical_residues_pose.residue(residue_index)
                current_residue_type = current_residue.type()
                residue_conformers[current_residue.name3()] = rosetta.core.pack.rotamer_set.bb_independent_rotamers(current_residue_type, True)

        ############################################
        # Add all other motif residues to fuzzball #
        ############################################

        # todo: break this function up...
        def add_motif_residues_to_current_fuzzball(current_conformer_fuzzball, motif_residue_list_untransformed, resnum_count, iteration=False, limit=0):
            """
            Add motif residues to the current fuzzball. This function checks for redundant residues and potential clashes
            with the ligand and any user-defined motif residues before adding a motif residue to the fuzzball.

            :param current_conformer_fuzzball: prody object of current fuzzball
            :param motif_residue_list_untransformed: list of motif residues as prody objects to be considered for addition
            to current_conformer_fuzzball
            :return:
            """

            # Motif residues added to fuzzball
            added = 0

            for motif_residue in self.pull_relevant_motif_residues_from_clusters(motif_residue_list_untransformed):

                if len(set(motif_residue.getAltlocs())) > 1:
                    # print(motif_residue.getData("contact_source")[0])
                    # print(motif_residue.getAltlocs())
                    continue

                motif_residue_transformed, constraint_atoms_dict = self.transform_residue_about_current_conformer(current_conformer, motif_residue, apply_mapping=existing)

                if motif_residue_transformed is None:
                    continue

                # --- Redundant residue check --- #

                current_resname = motif_residue_transformed.getResnames()[0]
                current_motif_constraint_atoms = constraint_atoms_dict['residue']['atom_names']
                current_motif_ligand_contact = constraint_atoms_dict['ligand']['atom_names'][0]

                current_motif_coord_list = [motif_residue_transformed.select(f'name {atom}').getCoords()[0] for atom in current_motif_constraint_atoms]
                current_motif_coords = np.asarray(current_motif_coord_list)

                current_motif_ca = motif_residue_transformed.select('name CA')

                if any([len(current_motif_coords) < 3, len(current_motif_ca) != 1]):
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

                    truthiness_list.append(
                        all([res_contact_atoms_match, ligand_contact_atom_match, contact_atoms_rmsd, ca_distance]))

                if any(truthiness_list):
                    print(f'\n{motif_residue_transformed.getData("contact_source")[0]} was found to be redundant!\n')
                    continue

                # --- Minimum fuzzball compatibility check --- #

                minimum_fuzzball_plus_current_motif = minimum_fuzzball_pose.clone()

                # Convert motif residue into PyRosetta Pose
                motif_pose = rosetta.core.pose.Pose()
                motif_prody_string = io.StringIO()
                prody.writePDBStream(motif_prody_string, motif_residue_transformed)
                rosetta.core.import_pose.pose_from_pdbstring(motif_pose, motif_prody_string.getvalue())

                # Add current motif residue to pose
                minimum_fuzzball_plus_current_motif.append_pose_by_jump(motif_pose, residues_in_minimum_fuzzball)

                # Score
                bsff_sfxn(minimum_fuzzball_plus_current_motif)
                e_edge = minimum_fuzzball_plus_current_motif.energies().energy_graph()
                interaction_energies = list()

                for i in range(1, residues_in_minimum_fuzzball + 1):
                    get_energy_edge = e_edge.find_energy_edge(i, residues_in_minimum_fuzzball + 1)

                    if get_energy_edge is not None:
                        current_edge = get_energy_edge.fill_energy_map()
                        interaction_energy = current_edge.dot(sfxn_weights)
                        interaction_energies.append(interaction_energy)

                # Only accept residues that do not clash with the minimum motif
                # Only score motif-ligand edge during iteration since inverse rotamers may be compatible
                if iteration:
                    if interaction_energies[0] > 25:
                        print(f'{motif_residue_transformed.getData("contact_source")[0]} clashes with the ligand!')
                        continue
                else:
                    if any([e > 25 for e in interaction_energies]):
                        print(f'{motif_residue_transformed.getData("contact_source")[0]} clashes with a defined residue!\n')
                        continue

                # --- Inverse Rotamer Match Compatability Filter --- #

                if iteration:

                    # Keep track of match positions that can accommodate current residue identity
                    compatible_match_positions = set()

                    # Get coords from matcher contact atoms
                    pyrosetta_motif_residue = motif_pose.residue(1).clone()
                    atm_1, atm_2, atm_3 = tuple(current_motif_constraint_atoms)

                    # todo: figure out how atom coords can become missing here...
                    try:
                        target_motif_coords = rosetta.numeric.xyzTransform_double_t(pyrosetta_motif_residue.xyz(atm_1),
                                                                                    pyrosetta_motif_residue.xyz(atm_2),
                                                                                    pyrosetta_motif_residue.xyz(atm_3),
                                                                                    )
                    except Exception as e:
                        print(e)
                        continue

                    motif_constraint_indicies = [pyrosetta_motif_residue.atom_index(atom) for atom in current_motif_constraint_atoms]
                    assert len(motif_constraint_indicies) == 3

                    # --- Try fuzzball conformer --- #

                    # Get all nearby backbone positions in match
                    match_neighbors = [position for position in range(1, match_pose_size) if pyrosetta_motif_residue.xyz('CA').distance(match_pose.residue(position).xyz('CA')) < 10]

                    # Test current motif rotamer against match neighbors
                    motif_bb_atoms = pyrosetta_motif_residue.mainchain_atoms()
                    motif_mainchain_coords = np.asarray([list(pyrosetta_motif_residue.xyz(atom)) for atom in motif_bb_atoms])

                    for match_position in match_neighbors:
                        match_neighbor_residue = match_pose.residue(match_position)
                        match_neighbor_residue_bb_atoms = match_neighbor_residue.mainchain_atoms()
                        match_mainchain_coords = np.asarray([list(match_neighbor_residue.xyz(atom)) for atom in match_neighbor_residue_bb_atoms])
                        mainchain_rmsd = prody.calcRMSD(match_mainchain_coords, motif_mainchain_coords)

                        if mainchain_rmsd < 2:
                            compatible_match_positions.add(match_position)

                    # --- Try other conformers --- #

                    # Get conformers for current residue type
                    current_rotamer_library = residue_conformers[pyrosetta_motif_residue.name3()]
                    sample_motif_rotamers = pyrosetta_motif_residue.clone()

                    # Iterate over rotamers and check backbone RMSD
                    for current_rotamer in current_rotamer_library:

                        # Apply torsions from rotamer to test motif resdiue
                        current_rotamer_chis = current_rotamer.chi()
                        sample_motif_rotamers.set_all_chi(current_rotamer_chis)

                        # Calculate rotation/translation by hand using first three residues of ligand
                        # todo: figure out how atom coords can become missing here...
                        try:
                            mobile_motif_coords = rosetta.numeric.xyzTransform_double_t(
                                sample_motif_rotamers.xyz(atm_1),
                                sample_motif_rotamers.xyz(atm_2),
                                sample_motif_rotamers.xyz(atm_3),
                                )
                        except Exception as e:
                            print(e)
                            continue

                        mobile_motif_inverse = mobile_motif_coords.inverse()

                        motif_rotation = target_motif_coords.R * mobile_motif_inverse.R
                        motif_translation = target_motif_coords.R * mobile_motif_inverse.t + target_motif_coords.t

                        # Apply transformation
                        sample_motif_rotamers.apply_transform_Rx_plus_v(motif_rotation, motif_translation)

                        # --- Repeat code, but whatever --- #

                        match_neighbors = [position for position in range(1, match_pose_size) if sample_motif_rotamers.xyz('CA').distance(match_pose.residue(position).xyz('CA')) < 4]

                        # Test current motif rotamer against match neighbors
                        motif_bb_atoms = sample_motif_rotamers.mainchain_atoms()
                        motif_mainchain_coords = np.asarray([list(sample_motif_rotamers.xyz(atom)) for atom in motif_bb_atoms])

                        for match_position in match_neighbors:
                            match_neighbor_residue = match_pose.residue(match_position)
                            match_neighbor_residue_bb_atoms = match_neighbor_residue.mainchain_atoms()
                            match_mainchain_coords = np.asarray([list(match_neighbor_residue.xyz(atom)) for atom in match_neighbor_residue_bb_atoms])
                            mainchain_rmsd = prody.calcRMSD(match_mainchain_coords, motif_mainchain_coords)

                            if mainchain_rmsd < 2:
                                compatible_match_positions.add(match_position)

                    if len(compatible_match_positions) == 0:
                        print(f'{motif_residue_transformed.getData("contact_source")[0]} cannot be accommodated by the current match!\n')
                        continue

                    # Use matcher 6-DOFs to determine accepted rotamers for rotamerset generation
                    dof_tuple, constraint_atoms_dict = self.generate_constraints.calculate_matcher_dofs(motif_residue_transformed, current_conformer)
                    contact_info = [current_motif_coord_list, [float(a) for a in dof_tuple],
                                    constraint_atoms_dict['residue']['atom_names'],
                                    constraint_atoms_dict['ligand']['atom_names']]

                    # Add residue identity to residue_match_mapping
                    for position in compatible_match_positions:

                        # Don't mess with CPG residues
                        new_resname = pyrosetta_motif_residue.name1()
                        if new_resname in ['C', 'P', 'G']:
                            continue

                        if new_resname in match_residue_map[position].keys():
                            match_residue_map[position][new_resname].append(contact_info)
                        else:
                            match_residue_map[position][new_resname] = [contact_info]

                # --- Add residue to fuzzball if all filters are passed --- #

                # Set resnum, occupany, coordset
                resnum_count += 1
                motif_residue_transformed.setResnums([resnum_count] * len(motif_residue_transformed))
                motif_residue_transformed.setAltlocs([' '] * len(motif_residue_transformed))
                motif_residue_transformed.setOccupancies([1] * len(motif_residue_transformed))
                motif_residue_transformed.setData('defined', [0] * len(motif_residue_transformed))

                # Add residue to fuzzball
                current_conformer_fuzzball += motif_residue_transformed
                print(f'\n{motif_residue_transformed.getData("contact_source")[0]} accepted as residue {current_conformer_fuzzball.numResidues()}!\n')

                # Remember transformed motifs
                transformed_motifs[current_resname].append((current_motif_constraint_atoms,
                                                            current_motif_ligand_contact,
                                                            current_motif_coords,
                                                            current_motif_ca))
                added += 1

                if (force_limit or not iteration) and added >= limit:
                    break

            return current_conformer_fuzzball, resnum_count

        # Apply banned residue filter to motif_residue_attributes_df
        banned_resnames = {'ALA', 'GLY', 'PRO', 'CYS'}
        working_resnames = list(set(self.resnames) - banned_resnames)
        self.motif_residue_attributes_df = self.motif_residue_attributes_df[self.motif_residue_attributes_df['resname'].isin(working_resnames)]

        # Get residues for current conformer
        current_conformer_motifs = self.motif_residue_attributes_df.groupby('conformer').get_group(current_conformer_name)
        current_conformer_motifs_by_contact = current_conformer_motifs.groupby('ligand_contact')

        # Maintain dict of transformed residues (keys resname) to check for motif residue redundancy
        transformed_motifs = dict.fromkeys(self.resnames, [])

        # --- Add hydrogen bonding motif residues to fuzzball --- #

        print('\nAdding hydrogen bonding motif residues to fuzzball...\n')

        for contact_atom in self.hbond_atom_set:
            contacts_df = current_conformer_motifs_by_contact.get_group(contact_atom)
            hbond_contacts_df = contacts_df[((contacts_df['hbond_sc'] < 0) | (contacts_df['hbond_bb_sc'] < 0)) & (contacts_df['total_score'] < 25)].sort_values(by=['total_score'])
            motif_prody_list = list(hbond_contacts_df['contact_source'].values)
            print(f'Considering {len(motif_prody_list)} residues for contact atom {contact_atom}...\n')

            current_conformer_fuzzball, resnum_count = add_motif_residues_to_current_fuzzball(current_conformer_fuzzball, motif_prody_list, resnum_count, iteration=iteration, limit=hbond_limit)

        # --- Add packing motif residues to fuzzball --- #

        print('Adding packing motif residues to fuzzball...\n')

        # Packed atoms are just atoms in the ligands that are not in defined_hbonding_atoms
        packed_atoms = set(current_conformer_motifs['ligand_contact'].values) - set(self.hbond_atom_set)

        # Calculate how many packing residues to add per atom
        calculated_limit = int((fuzzball_limit - current_conformer_fuzzball.numResidues() - 1) / len(packed_atoms))

        for contact_atom in packed_atoms:
            contacts_df = current_conformer_motifs_by_contact.get_group(contact_atom)
            packing_contacts_df = contacts_df[(contacts_df['hbond_sc'] == 0) & (contacts_df['total_score'] < 25)].sort_values(by=['total_score'])
            motif_prody_list = list(packing_contacts_df['contact_source'].values)

            current_conformer_fuzzball, resnum_count = add_motif_residues_to_current_fuzzball(current_conformer_fuzzball, motif_prody_list, resnum_count, iteration=iteration, limit=calculated_limit)

        # --- Output fuzzballs --- #

        if match_pdb:
            match_residue_map_curated = {int(key): value for key, value in match_residue_map.items() if len(value.keys()) > 0 and key not in motif_resnums}
            with open(os.path.join(self.fuzzball_dir, f'iter_{iteration_index}-fuzz_{fuzzball_index}-match_residue_map.pickle'), 'wb') as match_residue_pickle:
                pickle.dump(match_residue_map_curated, match_residue_pickle)

        fuzzball_name = iteration_name if iteration_name else current_conformer_name
        print(f'\n{current_conformer_fuzzball.numResidues() - 1} motif residues found for {current_conformer_name}!\n')

        fuzzball_identifier = f'{fuzzball_name}-iter_{iteration_index}-fuzz_{fuzzball_index}'
        current_conformer_fuzzball_path = os.path.join(self.fuzzball_dir, fuzzball_identifier)

        # .ag.npz representation (for iteration and contact sourcing)
        prody.saveAtoms(current_conformer_fuzzball, filename=current_conformer_fuzzball_path)

        # todo: Remove ProDy REMARKS from output
        # .pdb representation (for FeatureReporter)
        prody.writePDB(current_conformer_fuzzball_path, current_conformer_fuzzball)

        # Motif residue summary dataframe
        fuzzball_id_list = [motif.getData('contact_source')[0] for motif in current_conformer_fuzzball.iterResidues()]
        motif_residue_summary_df = self.motif_residue_attributes_df[(self.motif_residue_attributes_df['contact_source'].isin(fuzzball_id_list)) & (self.motif_residue_attributes_df['conformer'] == current_conformer_name)]

        # Left join to add resnums
        fuzzball_resnums = [{'resnum': res.getResnum(), 'contact_source': res.getData('contact_source')[0]} for res in
                            current_conformer_fuzzball.iterResidues() if res.getResname() is not self.chemical_component_identifier]
        fuzzball_resnums_df = pd.DataFrame(fuzzball_resnums)
        motif_residue_summary_df = motif_residue_summary_df.merge(fuzzball_resnums_df, on='contact_source', how='left')
        motif_residue_summary_df.to_csv(os.path.join(self.fuzzball_dir, f'{fuzzball_identifier}.csv'), index=False)

        return os.path.join(self.fuzzball_dir, f'{fuzzball_identifier}.pdb')

    def pull_relevant_motif_residues_from_clusters(self, unique_motif_id_set):
        """
        Pull motif residues from clusters

        :return:
        """
        # Pull residues (deepcopy!) out of clusters and store in a list
        motif_residue_list_untransformed = list()

        # Get (fragment, cluster) set for residues I need
        fragment_cluster_set = set()
        for id_string in unique_motif_id_set:
            id_string_split = re.split('[_-]', id_string)
            fragment_cluster_set.add((int(id_string_split[1]), int(id_string_split[3]), int(id_string_split[4])))

        for fragment, cat_cluster, cluster in fragment_cluster_set:

            current_cluster = self.cluster_dict[f'Fragment_{fragment}'][cat_cluster][cluster]
            residues_in_cluster = set(current_cluster.getData('contact_source'))
            residues_to_pull = residues_in_cluster & set(unique_motif_id_set)

            if len(residues_to_pull) > 0:

                selection_string = ' '.join([f'`{residue_to_pull}`' for residue_to_pull in residues_to_pull])
                pulled_residue_prody = current_cluster.select(f'contact_source {selection_string}').getHierView()

                for residue in pulled_residue_prody.iterResidues():
                    yield residue.copy()
                    # motif_residue_list_untransformed.append(residue.copy())

        # return motif_residue_list_untransformed

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

class FailedAssemblyError(Exception):
    pass
