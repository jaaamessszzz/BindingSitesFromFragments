#!/usr/bin/env python3

"""
Utility functions!
"""
import os
import io
import prody
import numpy as np
import re

from rdkit import Chem
from pyrosetta import rosetta

# --- Path/file traversal --- #

def directory_check(dir, base_only=False):
    for subdir in os.listdir(dir):
        path = os.path.join(dir, subdir)
        if os.path.isdir(path):
            if base_only:
                yield subdir
            else:
                yield path


def pdb_check(dir, base_only=False, conformer_check=False, extension='.pdb'):
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if path.endswith(extension):
            if not conformer_check or (conformer_check and 'conformers' not in re.split('\.|_|-', file)):
                if base_only:
                    yield file
                else:
                    yield path


def processed_check(processed_dir, pdb, rejected_list):
    pdbid = pdb.split('.')[0].upper()

    in_reject_list = any([pdb == reject for reject in rejected_list])
    # already_processed = os.path.exists(os.path.join(processed_dir, '{}_processed.pdb'.format(pdbid)))
    already_processed = any([pdbid == file[:4].upper() for file in os.listdir(processed_dir)])

    return any([in_reject_list, already_processed])

# --- Common manipulations --- #

def minimum_contact_distance(a_H, b_H, return_indices=False, strip_H=True):
    """
    Calculates the minimum distance between two sets of coordinates
    :param a_H: prody object of first set (rows of dist matrix)
    :param b_H: prody object of second set (columns of dist matrix)
    :param return_indices: boolean, whether or not to return row and column indicies of atoms with min distance in matrix
    :return: minimum distance in angstroms
    """

    if strip_H:
        a = a_H.select('not hydrogen').getCoords()
        b = b_H.select('not hydrogen').getCoords()

    else:
        a = a_H.getCoords()
        b = b_H.getCoords()

    ligand_residue_distance_matrix = prody.buildDistMatrix(a, b)

    # Find minimum score in matrix
    row_min_indicies = np.amin(ligand_residue_distance_matrix, axis=0)
    ligand_index = np.argmin(row_min_indicies, axis=0)
    residue_index = np.argmin(ligand_residue_distance_matrix, axis=0)

    column_index_low = ligand_index
    row_index_low = residue_index[column_index_low]

    # Contact distance
    if return_indices:
        return ligand_residue_distance_matrix.item(row_index_low, column_index_low), row_index_low, column_index_low
    else:
        return ligand_residue_distance_matrix.item(row_index_low, column_index_low)


def find_conformer_and_constraint_resnums(pdb_name):
    """
    Generates name of ideal binding motif from match PDB name
    :param pdb_name:
    :return:
    """
    pdb_split = re.split('_|-|\.', pdb_name)

    ligand_name_index = 5
    conformer_id_index = 6

    if len(pdb_split[4]) == 4:
        ligand_name_index += 1
        conformer_id_index += 1

    ligand_name = pdb_split[ligand_name_index]
    conformer_id = pdb_split[conformer_id_index]
    conformer_name = '{}_{}'.format(ligand_name, conformer_id)

    constraint_resnum_block = re.split('-|\.', pdb_name)[1]
    constraint_resnums = [int(a) for a in constraint_resnum_block.split('_') if a != ''][1:]

    return conformer_name, constraint_resnums


def determine_matched_residue_positions(match_pdb_path):
    """
    Parse the filename of the match PDB to determine IDs and positions of match residues
    :return:
    """
    positions_block = os.path.basename(os.path.normpath(match_pdb_path)).split('_')[2]
    resnames = [a for a in re.split("[0-9]*", positions_block) if a]
    resnums = [int(a) for a in re.split("[a-zA-Z]*", positions_block) if a]

    return [(a, b) for a, b in zip(resnames, resnums)]


def clean_pdb(input_pdb, ligand_code=None, resnum=None):
    """
    Clean PDB. Roland: MSE to MET records, CSE to CYS records, discarding alternative conformations, setting all atom
    occupancies to 1.0, discarding all residues with any missing main chain atom(s), removing all ligands but the
    target, and ensuring internal Rosetta numbering

    :param input_pdb: prody of PDB to clean
    :param ligand_code: chemical component identifier for ligand to keep
    :param resnum: position of ligand in input_pdb to keep
    :return:
    """

    canon_residues = ['ALA', 'CYS', 'SEC', 'ASP', 'GLU',
                      'PHE', 'GLY', 'HIS', 'ILE', 'LYS',
                      'LEU', 'MET', 'MSE', 'ASN', 'PRO',
                      'GLN', 'ARG', 'SER', 'THR', 'VAL',
                      'TRP', 'TYR']

    # Necessary due to prody things...
    def _add_to_output(output_atoms, residue):
        if len(output_atoms) != 0:
            output_atoms = output_atoms + residue.copy()
        else:
            output_atoms = residue.copy()
        return output_atoms

    # Pirated from Generate_Motif_Residues
    def _fix_mse_sec(representative_residue, resname):
        """
        Takes MSE/SEC and turns it into MET/CYS
        :param representative_residue: representative_residue with Se
        :return: representative_residue with S
        """
        # Find index of SE
        res_elements = representative_residue.getElements()

        # If SE is missing due to missing atoms, just return the residue as is
        if 'SE' not in res_elements:
            return representative_residue

        # Set SE to S
        seleno_index = [e for e in res_elements].index('SE')
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

    # This selects Altloc A if there are alternate locations at all... makes things easy
    cleanish_pdb = input_pdb.select('(protein or hetero or nucleic) and not water').copy()
    hv = cleanish_pdb.getHierView()

    output_atoms = prody.atomgroup.AtomGroup('Output')
    res_count = 1
    removed_residues = []

    for chain in hv:
        for residue in chain:
            if residue.getResnum() == resnum and residue.getResname() == ligand_code:
                residue.setResnum(res_count)
                res_count += 1
                output_atoms = _add_to_output(output_atoms, residue)

            # Check Backbone atoms, else don't add to output atomgroup
            elif all(atom in residue.getNames() for atom in ['N', 'C', 'CA']) and residue.getResname() in canon_residues:
                residue.setResnum(res_count)
                res_count += 1

                if residue.getResname() in ['MSE', 'SEC']:
                    residue = _fix_mse_sec(residue, residue.getResname())

                residue.setOccupancies([1] * len(residue))
                output_atoms = _add_to_output(output_atoms, residue)

            elif residue.getResname() in ['A', 'T', 'C', 'G', 'U']:
                residue.setResnum(res_count)
                res_count += 1
                output_atoms = _add_to_output(output_atoms, residue)

            else:
                print('Removed {}'.format(residue))
                removed_residues.append('{0}{1}'.format(residue.getResname(), residue.getResnum()))

    return output_atoms, removed_residues

# --- Type conversion --- #

def RDKit_Mol_from_ProDy(prody_instance, removeHs=True):
    """
    Creates an RDKit Mol object from a ProDy AtomGroup instance
    :return:
    """
    residue_io = io.StringIO()
    prody.writePDBStream(residue_io, prody_instance)

    return Chem.MolFromPDBBlock(residue_io.getvalue(), removeHs=removeHs)


# --- PyRosetta --- #

def find_binding_motif_clashes(pose, sfxn, residue_index_list):
    """
    Returns the largest fa_rep value for any of the motif residues and ligand where all other residues in the match
    scaffold are mutated to ALA
    :param pose:
    :param residue_index_list:
    :return:
    """
    ligand_index = pose.size()
    motif_and_ligand_idx = residue_index_list + [ligand_index]

    sfxn(pose)  # Redundant
    edges = pose.energies().energy_graph()

    # --- Calculate max fa_rep for each motif residue --- #
    fa_rep_list = []
    fa_sol_list = []

    for motif_res in motif_and_ligand_idx:
        for res in range(1, ligand_index):
            if res not in motif_and_ligand_idx:
                current_edge = edges.find_energy_edge(motif_res, res)
                if current_edge is not None:
                    current_edge.fill_energy_map()
                    fa_rep_list.append(current_edge[rosetta.core.scoring.fa_rep])
                    fa_sol_list.append(current_edge[rosetta.core.scoring.fa_sol])

    return max(fa_rep_list), max(fa_sol_list)

def get_rotation_and_translation(mobile, target):
    """

    :param mobile: rosetta.numeric.xyzTransform_double_t
    :param target: rosetta.numeric.xyzTransform_double_t
    :return:
    """

    # Calculate rotation/translation by hand using first three residues of ligand
    mobile_inverse = mobile.inverse()

    mobile_rotation = target.R * mobile_inverse.R
    mobile_translation = target.R * mobile_inverse.t + target.t

    # Apply transformation
    return mobile_rotation, mobile_translation


