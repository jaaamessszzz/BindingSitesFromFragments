#!/usr/bin/env python3

import os
import sys
import numpy as np
import prody
import pprint
import re
from .utils import *

class Generate_Motif_Residues():
    """
    This is a class for generating and manipulating motif residues for hypothetical ligand binding sites
    """
    def __init__(self, cluster_path, motif_cluster_yaml):
        self.cluster_path = cluster_path
        self.fragment_cluster_list = motif_cluster_yaml
        self.expected_atoms = {'ALA': 5, 'CYS': 6, 'ASP': 8, 'GLU': 9, 'PHE': 11, 'GLY': 4, 'HIS': 10, 'ILE': 8,
                               'LYS': 9, 'LEU': 8, 'MET': 8, 'ASN': 8, 'PRO': 7, 'GLN': 9, 'ARG': 11, 'SER': 6,
                               'THR': 7, 'VAL': 7, 'TRP': 14, 'TYR': 12}

    def generate_motif_residues(self):
        """
        Define motif residues from clusters from clusters outlined in a yaml file
        :return: 
        """

        # Assemble dict to store list of prody residue instances for each cluster
        fragment_prody_dict = {}

        # For each cluster in the cluster list
        for fragment in self.fragment_cluster_list:
            fragment_prody_dict[fragment] = {}
            fragment_cluster_path = os.path.join(self.cluster_path, fragment)

            # Make a list for each cluster for a given fragment
            for cluster_number in self.fragment_cluster_list[fragment]:
                fragment_prody_dict[fragment][int(cluster_number)] = []

            # Check pdbs for given fragment and add to appropriate cluster list
            for fragment_pdb in pdb_check(fragment_cluster_path, base_only=True):
                cluster_number = int(re.split('-|_', fragment_pdb)[1])

                if cluster_number in self.fragment_cluster_list[fragment]:
                    fragment_prody_dict[fragment][cluster_number].append(prody.parsePDB(os.path.join(fragment_cluster_path, fragment_pdb)).select('not hydrogen'))

        pprint.pprint(fragment_prody_dict)

        # Okay, cool. Now that I have all of my clusters neatly organized in a dict, I can go through and generate motif
        # residues from each of the clusters based on residue type

        fragment_output_dir = os.path.join(os.path.split(self.cluster_path)[0], 'Representative_Residue_Motifs')
        os.makedirs(fragment_output_dir, exist_ok=True)

        # Start going through fragments and clusters
        for fragment in fragment_prody_dict:
            motif_count = 1

            for cluster in fragment_prody_dict[fragment]:
                # For each type of residue in cluster
                residue_types = list(set([residue.getResnames()[0] for residue in fragment_prody_dict[fragment][cluster]]))

                for res in residue_types:
                    # Calculate average residue coordinates
                    cluster_residues = [residue for residue in fragment_prody_dict[fragment][cluster]
                                        if all([residue.getResnames()[0] == res, len(residue) == self.expected_atoms[residue.getResnames()[0]]])]
                    average_coordinates = np.mean([a.getCoords() for a in cluster_residues], axis=0)

                    # Select residue closest to average
                    rmsd_list = [prody.calcRMSD(average_coordinates, residue.getCoords()) for residue in cluster_residues]
                    representative_residue = cluster_residues[rmsd_list.index(min(rmsd_list))]

                    # Debugging
                    print(rmsd_list)
                    print(min(rmsd_list))
                    print(average_coordinates)
                    print(representative_residue.getCoords())
                    print(prody.calcRMSD(average_coordinates, representative_residue.getCoords()))

                    # Output residue
                    prody.writePDB(os.path.join(fragment_output_dir,
                                                '{0}-Cluster_{1}-Motif-{2}-{3}.pdb'.format(fragment, cluster, motif_count, res)),
                                   representative_residue)

                motif_count += 1

    def define_second_shell_contants(self):
        """
        Identify second shell contacts with motif residues
        :return: 
        """
        pass

class Generate_Binding_Sites():
    """
    This is a class for combining specified groups of binding motifs into hypothetical binding sites
    User input required to deal with residues that obviously do not get along with other motif residues or the ligand
    e.g. steric clashing, non-favorable interactions
    """
    def __init__(self, residue_groups):
        self.residue_groups = residue_groups

    def generate_binding_sites(self):
        """
        sup.
        :return: 
        """
        pass

    def generate_constraints(self):
        """
        Generate matcher constraints for a given hypothetical binding site.
        :return: 
        """
        pass