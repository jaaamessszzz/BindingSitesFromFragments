#!/usr/bin/env python3

"""
Utility functions!
"""
import os

def directory_check(dir):
    for subdir in os.listdir(dir):
        path = os.path.join(dir, subdir)
        if os.path.isdir(path):
            yield path


def pdb_check(dir):
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if path.endswith('.pdb'):
            yield path

def processed_check(processed_dir, pdb, rejected_list):
    pdbid = pdb.split('.')[0].upper()

    in_reject_list = any([pdb == reject for reject in rejected_list])
    # already_processed = os.path.exists(os.path.join(processed_dir, '{}_processed.pdb'.format(pdbid)))
    already_processed = any([pdbid == file[:4] for file in os.listdir(processed_dir)])

    return any([in_reject_list, already_processed])