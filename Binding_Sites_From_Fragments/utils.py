#!/usr/bin/env python3

"""
Utility functions!
"""
import os
import prody
import numpy as np

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


def pdb_check(dir, base_only=False):
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if path.endswith('.pdb'):
            if base_only:
                yield file
            else:
                yield path


def processed_check(processed_dir, pdb, rejected_list):
    pdbid = pdb.split('.')[0].upper()

    in_reject_list = any([pdb == reject for reject in rejected_list])
    # already_processed = os.path.exists(os.path.join(processed_dir, '{}_processed.pdb'.format(pdbid)))
    already_processed = any([pdbid == file[:4] for file in os.listdir(processed_dir)])

    return any([in_reject_list, already_processed])


def minimum_contact_distance(a, b, return_indices=False):
    """
    Calculates the minimum distance between two sets of coordinates
    :param a: coordinates of first set
    :param b: coordinates of second set
    :param return_indices: boolean, whether or not to return row and column indicies of atoms with min distance in matrix
    :return: minimum distance in angstroms
    """
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