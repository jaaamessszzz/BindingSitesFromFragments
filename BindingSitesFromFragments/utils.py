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

def generate_new_project():
    """
    Generate a new project for a specific ligand. Sets up all necessary directory and template files
    :return: 
    """
    pass


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

    constraint_resnum_block = re.split('-|\.', pdb_name)[1][:-2]
    constraint_resnums = [int(a) for a in constraint_resnum_block.split('_') if a != ''][1:]

    return conformer_name, constraint_resnums


def RDKit_Mol_from_ProDy(prody_instance, removeHs=True):
    """
    Creates an RDKit Mol object from a ProDy AtomGroup instance
    :return:
    """
    residue_io = io.StringIO()
    prody.writePDBStream(residue_io, prody_instance)

    return Chem.MolFromPDBBlock(residue_io.getvalue(), removeHs=removeHs)