#!/usr/bin/env python3

"""
Utility functions!
"""
import os

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