#!/usr/bin/env python3

import os
import sys
import numpy as np
import prody
import pprint
import re
import itertools
from .utils import *

class Generate_Motif_Residues():
    """
    This is a class for generating and manipulating motif residues for hypothetical ligand binding sites
    
    For each ligand conformer during motif residue generation, I will calculate and output:
    1.  A sparse matrix of motif residues that clash
    2.  Residue-ligand interaction scores as calculated with Rosetta
    3.  Constraint file for each residue-ligand interaction
    """
    def __init__(self, cluster_path, motif_cluster_yaml):
        self.cluster_path = cluster_path
        self.fragment_cluster_list = motif_cluster_yaml
        # todo: accommodate C-term residues...
        self.expected_atoms = {'ALA': 5, 'CYS': 6, 'ASP': 8, 'GLU': 9, 'PHE': 11, 'GLY': 4, 'HIS': 10, 'ILE': 8,
                               'LYS': 9, 'LEU': 8, 'MET': 8, 'ASN': 8, 'PRO': 7, 'GLN': 9, 'ARG': 11, 'SER': 6,
                               'THR': 7, 'VAL': 7, 'TRP': 14, 'TYR': 12}


    def generate_motif_residues(self):
        """
        Define motif residues from clusters outlined in a user defined yaml file for a single ligand conformation.
        A yaml file will also be generated to keep track of which residues clash with each conformer. This way I can 
        output all representative motif residues once and just ignore the ones that clash for each conformer
        * Inputs/User_Inputs/Motif_Clusters.yaml
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

        # Okay, cool. Now that I have all of my clusters neatly organized in a dict, I can go through and generate motif
        # residues from each of the clusters based on residue type

        fragment_output_dir = os.path.join(os.path.split(self.cluster_path)[0], 'Representative_Residue_Motifs', 'Unfiltered_Residues')
        os.makedirs(fragment_output_dir, exist_ok=True)

        # Start going through fragments and clusters
        total_motif_count = 1
        for fragment in fragment_prody_dict:
            motif_count = 1

            for cluster in fragment_prody_dict[fragment]:
                # For each type of residue in cluster
                residue_types = list(set([residue.getResnames()[0] for residue in fragment_prody_dict[fragment][cluster]]))

                for res in residue_types:
                    # Calculate average residue coordinates
                    cluster_residues = [residue for residue in fragment_prody_dict[fragment][cluster]
                                        if all([residue.getResnames()[0] == res, len(residue) == self.expected_atoms[residue.getResnames()[0]]])]

                    if len(cluster_residues) > 0:
                        average_coordinates = np.mean([a.getCoords() for a in cluster_residues], axis=0)

                        # Select residue closest to average
                        rmsd_list = [prody.calcRMSD(average_coordinates, residue.getCoords()) for residue in cluster_residues]
                        representative_residue = cluster_residues[rmsd_list.index(min(rmsd_list))]

                        # Output residue
                        prody.writePDB(os.path.join(fragment_output_dir, '{4}-{5}-{0}-Cluster_{1}-Motif-{2}-{3}'.format(fragment, cluster, motif_count, res, total_motif_count, len(fragment_prody_dict[fragment][cluster]))),
                                       representative_residue)
                        total_motif_count += 1

                motif_count += 1

        self.filter_representative_residues(fragment_output_dir)


    def filter_representative_residues(self, unfiltered_dir, cutoff_distance=2):
        """
        Filters out residues that are too close to the ligand molecule and outputs the rest as PDBs
        :return: 
        """
        target_molecule = os.path.split(self.cluster_path)[0]
        # todo: update importing target molecule prody to accommodate conformer libraries
        target_molecule_prody = prody.parsePDB(os.path.join(target_molecule, 'Inputs', '{}.pdb'.format(target_molecule))).select('not hydrogen')

        # Generate and output a yaml file for each conformer and list which residues clash
        # This way I can output all representative motif residues once and just ignore the ones that clash for each conformer
        for residue in pdb_check(unfiltered_dir):
            residue_prody = prody.parsePDB(residue).select('not hydrogen')
            if minimum_contact_distance(residue_prody.getCoords(), target_molecule_prody.getCoords()) >= cutoff_distance:
                prody.writePDB(os.path.join(os.path.split(unfiltered_dir)[0], '{}'.format(os.path.basename(os.path.normpath(residue)))), residue_prody)


    def generate_clash_matrix(self, clashing_cutoff=2):
        """
        Generates a sparse matrix for all representative motif residues that clash with each other. This is determined
        based on the distance of the closest atom-atom interaction between two residues.
        Matrix is output as a .csv that can be imported with numpy.
        
        :param clashing_cutoff: distance cutoff for clashing residues in angstroms. Default set to 2.
        :return: 
        """
        pass


    def score_residue_ligand_interactions(self, ligand_pdb):
        """
        Score each unique residue-ligand interaction with Rosetta and export to .csv (fa_atr, hbond_sc, fa_elec) 
        
        :param ligand_pdb: PDB of the current ligand conformer to use for scoring.
        :return: 
        """
        for motif_residue in pdb_check(os.path.join(os.path.split(self.cluster_path)[0], 'Representative_Residue_Motifs')):
            residue_prody = prody.parsePDB(motif_residue)
        # Import ligand pdb as generated by molfile_to_params.py (for loop if conformer library has been generated)
        # For each motif residue in Representative_Residue_Motifs directory, combine ligand and residue
        # Output ligand-residue pairs to a new directory in Hypothetical_Binding_Sites directory
        # Write all PDBs to a text file so I can score them all with score_jd2 -in:file:l <text_file_with_list.txt>
        # Calculate scores for all PDBs in this new directory (batch process, don't forget .params)
        # Process score file into a .csv that can be imported by numpy as a matrix


    def generate_constraint(self, binding_site_description, prody_binding_site):
        """
        Generate matcher constraint for a given residue-ligand interaction using information from clusters

        :return: 
        """
        # For residue in binding site
        for residue in prody_binding_site:
            pass

            # I need to determine the closest atom-atom contacts and two additional atoms for determining bond torsions and angles
            # Get ideal distance, angle, and torsion values from prody residue
            # Get tolerance from clusters; I'm going to try making the tolerance +/- 1 SD of cluster values
            # Set penalty to whatever since it isn't used by the matcher
            # Set distance to 0, otherwise periodicity (per) to 360
            # Set sample number so that samples are made every 5 degrees or every 0.1A


    def define_second_shell_contraints(self):
        """
        Identify second shell contacts with motif residues
        :return: 
        """
        # This is most likely just going to take the PDB for a given representative motif residue and pull out a 10 A
        # shell or something like that... this is really going to be useful if I can use it on super tight clusters
        # so I can pull out second shell interaction patterns
        pass


class Generate_Binding_Sites():
    """
    This is a class for combining specified groups of binding motifs into hypothetical binding sites
    User input required to deal with residues that obviously do not get along with other motif residues or the ligand
    e.g. steric clashing, non-favorable interactions
    """
    def __init__(self, user_defined_dir, residue_groups, hypothetical_binding_sites):
        self.user_defined_dir = user_defined_dir
        self.residue_groups = residue_groups
        self.hypothetical_binding_sites = hypothetical_binding_sites


    def calculate_energies_and_rank(self):
        """
        Calculate total binding site interaction energies for all possible conformations of representative binding 
        motifs and rank.
        
        :return: 
        """
        pass


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
                    for pdb in pdb_check(os.path.join(self.user_defined_dir, 'Representative_Residue_Motifs'), base_only=True):
                        if pdb.split('-')[0] == str(residue):
                            binding_residue_list.append(prody.parsePDB(os.path.join(self.user_defined_dir, 'Representative_Residue_Motifs', pdb)))
                            cluster_sum += int(pdb.split('-')[1])

                # Simple distance check for all residues to prevent clashing
                clashing = False
                cutoff_distance = 2

                for outer_index, a in enumerate(binding_residue_list):
                    for inner_index, b in enumerate(binding_residue_list[(outer_index+1):]):
                        if minimum_contact_distance(a.getCoords(), b.getCoords()) < cutoff_distance:
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
                    self.generate_constraints(binding_site_description, complete_binding_site)

                    # Output PDB
                    prody.writePDB(os.path.join(output_path, binding_site_description), complete_binding_site)

    def generate_constraints(self):
        """
        Generate a constraint file for a complete hypothetical binding site by concatinating previously generated
        residue-ligand constraints
        :return: 
        """
        pass