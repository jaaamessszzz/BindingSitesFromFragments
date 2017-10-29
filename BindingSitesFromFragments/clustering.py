#!/usr/bin/env python3

import os
import sys
import numpy as np
import prody
import pandas as pd
import pprint
from .utils import *

class Cluster():
    """
    This class is responsible for taking aligned fragment PDBs and identifying representative contacts. 
    """

    def __init__(self, processed_PDBs_dir, distance_cutoff, weights):
        self.processed_PDBs_dir = processed_PDBs_dir
        self.weights = weights
        self.distance_cutoff = distance_cutoff
        self.pdb_object_list = self._import_pdbs()
        self.clusters = None
        self.df = None

    def _import_pdbs(self):
        """
        For each fragment ensemble, converts each residue in all processed PDBs into objects with representative 
        matrices and other relevant information
        
        :return: 
        """
        processsed_residue_list = []
        for pdb in pdb_check(self.processed_PDBs_dir):
            pdb_info = os.path.basename(os.path.normpath(pdb))
            prody_protein = prody.parsePDB(pdb).select('protein and not hetatm')
            # Check that residues exist within cutoff distance provided in alignments, otherwise pass
            if prody_protein == None:
                continue
            else:
                prody_protein_hv = prody_protein.getHierView()

            prody_ligand = prody.parsePDB(pdb).select('hetatm and resname {}'.format(pdb_info.split('_')[1]))

            # Iterate over residues in contacts and generate representative vector with weights applied
            for chain in prody_protein_hv:
                processsed_residue_list += [fragment_PDB(residue, pdb_info, prody_ligand, self.distance_cutoff, self.weights) for residue in chain]

        processsed_residue_list_cleaned = [residue for residue in processsed_residue_list if residue.vector is not None]

        return processsed_residue_list_cleaned

    def cluster_scipy(self, display_dendrogram=False):
        """
        Using Hierarchical clustering through scipy to identify representative contacts. 
        :return: 
        self.clusters: list of indicies corresponding to cluster numbers for each element in self.pdb_object_list
        """
        from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
        from matplotlib import pyplot as plt

        vector_array = np.asarray([thingy.vector for thingy in self.pdb_object_list])
        Z = linkage(vector_array, 'ward')

        if display_dendrogram == True:
            # Generate dendrogram
            # todo: CHANGE THIS. Taken straight from https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial
            fig = plt.figure(figsize=(25, 10))
            plt.title('Hierarchical Clustering Dendrogram')
            plt.xlabel('sample index')
            plt.ylabel('distance')
            ax = fig.gca()
            ax.set_ylim(5)
            dendrogram(
                Z,
                leaf_rotation=90.,  # rotates the x axis labels
                leaf_font_size=8.,  # font size for the x axis labels
            )
            plt.show()

        self.distance_cutoff = 2**0.5
        self.clusters = fcluster(Z, self.distance_cutoff, criterion='distance')

    def generate_output_directories(self, base_dir, fragment_path):
        residue_cluster_pairs = [(cluster, residue)for cluster, residue in zip(self.clusters, self.pdb_object_list)]
        fragment = os.path.basename(os.path.normpath(fragment_path))

        dict_list = []

        for cluster in set(self.clusters):
            # Generate list of fragment_PDB objects containing all relevant information on residue-ligand contacts in a cluster
            residue_list = [pair for pair in residue_cluster_pairs if pair[0] == cluster]
            fragment_cluster_path = os.path.join(base_dir,
                                     'Cluster_Results',
                                     fragment
                                     )

            os.makedirs(fragment_cluster_path, exist_ok=True)

            # Generate report on cluster qualities
            dict_list.append(self.generate_report_row(residue_list, cluster))

            # Output residue PDBs for each cluster
            count=0
            for index, residue in enumerate(residue_list):
                if residue[1].prody_residue.getResnames()[0] in ['ALA', 'CYS', 'SEC', 'ASP', 'GLU', 'PHE',
                                                                 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET',
                                                                 'MSE', 'ASN', 'PRO', 'GLN', 'ARG', 'SER',
                                                                 'THR', 'VAL', 'TRP', 'TYR']:
                    prody.writePDB(os.path.join(fragment_cluster_path, 'Cluster_{0}-Residue_{1}-{2}'.format(cluster, str(index), residue[1].pdb_info)), residue[1].prody_residue)

        # Export .csv
        self.df = pd.DataFrame(dict_list)
        self.df.to_csv(os.path.join(fragment_cluster_path,'{}_report.csv'.format(fragment)))

    def generate_report_row(self, residue_list, cluster_number):
        """
        Generates reports on each cluster
        I guess the first step would be to basically output the information contained in each of the representative vectors?
        :param residue_list: list of residue objects for a given cluster index
        :param cluster_number: cluster index, a bit redundant since I can get this from residue_list[0][0]
        :return: 
        """
        # Determine spread in min ligand-residue atom contact unit vectors
        atom_atom_centroid_vector = np.mean([residue[1].contact_unit_vector for residue in residue_list], axis=0)
        atom_atom_angle_deviations = [np.arccos(np.dot(atom_atom_centroid_vector, residue[1].contact_unit_vector)) for residue in residue_list]
        atom_atom_mean = np.degrees(np.mean(atom_atom_angle_deviations))
        atom_atom_SD = np.degrees(np.std(atom_atom_angle_deviations))

        # 20170901 - I've determined that the atom-centroid mean and SD aren't useful metrics...
        # Determine spread in unit vectors connecting min ligand contact atom and residue centroids
        # atom_centroid_vectors = [(residue[1].residue_center - residue[1].ligand_contact_atom.getCoords()[0]) /
        #                          np.linalg.norm(residue[1].residue_center - residue[1].ligand_contact_atom.getCoords()[0])
        #                          for residue in residue_list]
        # atom_centroid_centroid_vector = np.mean(atom_centroid_vectors, axis=0)
        #
        # atom_centroid_angle_deviations = [np.arccos(np.dot(atom_centroid_centroid_vector, atom_centroid_vector)) for atom_centroid_vector in atom_centroid_vectors]
        # atom_centroid_mean = np.degrees(np.mean(atom_centroid_angle_deviations))
        # atom_centroid_SD = np.degrees(np.std(atom_centroid_angle_deviations))

        return {'cluster_index': cluster_number,
                'cluster_members': len(residue_list),
                'cutoff': self.distance_cutoff,
                'atom-atom_mean': atom_atom_mean,
                'atom-atom_SD': atom_atom_SD
                }

    def automated_cluster_selection(self):
        """
        Automatically selects clusters based on number of members in each clusters. Clusters >1 SD above mean will be 
        selected for motif residue pool.
        :return: 
        """
        # Bruh. Paths.
        motif_yaml_path = os.path.join(self.processed_PDBs_dir.split('/')[0], 'Inputs', 'User_Inputs', 'Motif_Clusters.yml')
        current_fragment = os.path.basename(os.path.normpath(self.processed_PDBs_dir))

        # Select worthy clusters
        sorted_df = self.df.sort_values('cluster_members',ascending=False)

        cluster_size_mean = self.df.agg({'cluster_members': 'mean'})[0]
        cluster_size_SD = self.df.agg({'cluster_members': 'std'})[0]

        print('Cluster size mean:', cluster_size_mean)
        print('Cluster size SD:', cluster_size_SD)

        selected_cluster_rows = sorted_df[sorted_df.cluster_members > (cluster_size_mean + cluster_size_SD)]

        print('Selected clusters:')
        pprint.pprint(selected_cluster_rows)

        with open(motif_yaml_path, 'a') as motif_yaml:
            motif_yaml.write('{}:\n'.format(current_fragment))
            for index, row in selected_cluster_rows.iterrows():
                motif_yaml.write('- {}\n'.format(int(row['cluster_index'])))

class fragment_PDB():
    """
    Class for holding all information related to a processed fragment PDB
    
    :param pdb_info:
    :param prody_pdb:
    :param prody_residue:
    :param pdbid:
    :param vector:
    """
    def __init__(self, prody_residue, pdb_info, prody_ligand, distance_cutoff, weights):
        self.prody_residue = prody_residue
        self.prody_ligand = prody_ligand
        self.pdb_info = pdb_info
        self.residue_center = prody.calcCenter(prody_residue)
        self.ligand_center = prody.calcCenter(prody_ligand)
        self.distance_cutoff = distance_cutoff
        self.weights = weights
        self.vector = self.process_residue_into_vector()
        # For cluster evaluation, assigned in process_residue_into_vector()
        # Not a great practice, but meh...
        # self.residue_contact_atom
        # self.ligand_contact_atom
        # self.min_contact_distance
        # self.contact_unit_vector

    def process_residue_into_vector(self):
        """
        Converting each residue into a representative vector

        Metrics I care about for each residue:
        * Contact distance
        * Residue characteristics 
            * Amino acid identity or degenerate amino acid groups?
            * Amino acid chemical characteristics e.g. [Sigma Amino Acid Reference Chart]
              (http://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html)
        * Position of residue relative to fragment
            * Vector from fragment centroid to {residue centroid | closest residue atom }
        * Backbone or side chain
        
         {Angstroms} - X component, unit vector from fragment centroid to closest residue atom 
         {Angstroms} - Y component, unit vector from fragment centroid to closest residue atom 
         {Angstroms} - Z component, unit vector from fragment centroid to closest residue atom
         
         {Angstroms} - Normalized contact distance  
         { 0 | 1 }   - Backbone contact OR Sidechain contact
         { 0 | 1 }   - Ligand Polar Contact OR Ligand Non-polar Contact
         
         { 0 | 1 }   - Side chain has hydrogen bond donor/acceptor (DEHKNQRSTY)
         
         { 0 | 1 }   - Hydrophobic, aliphatic (AILV)
         { 0 | 1 }   - Hydrophobic, aromatic (FWY)
         { 0 | 1 }   - Polar (NCQMST)
         { 0 | 1 }   - Charged, Acidic (DE)
         { 0 | 1 }   - Charged, Basic (HKR)
         { 0 | 1 }   - Glycine
         { 0 | 1 }   - Proline
         { 0 | 1 }   - Backbone carbonyl
         { 0 | 1 }   - Backbone amino
         { 0 | 1 }   - Backbone C/CA
        
        :return: 
        """

        min_contact_distance, row_index_low, column_index_low = minimum_contact_distance(self.prody_residue, self.prody_ligand, return_indices=True)

        if min_contact_distance > self.distance_cutoff:
            return None

        else:
            # todo: use copies of prody selections so I don't do this base index BS
            residue_base_index = self.prody_residue.getIndices()[0]
            ligand_base_index = self.prody_ligand.getIndices()[0]

            self.residue_contact_atom = self.prody_residue.select('index {}'.format(residue_base_index + row_index_low))
            self.ligand_contact_atom = self.prody_ligand.select('index {}'.format(ligand_base_index + column_index_low))

            # Save min contact residue and ligand atom indicies for evaluating cluster quality later
            self.min_contact_distance = min_contact_distance

            # Residue Contact Type
            residue_contact_type = 0 if self.residue_contact_atom.getNames()[0] in ['C', 'CA', 'N', 'O'] else 1

            # Ligand Contact Type
            ligand_contact_type = 1 if self.ligand_contact_atom.getNames()[0][0] in ['C'] else 0

            # Vector from fragment centroid to closest residue atom
            contact_vector = (self.residue_contact_atom.getCoords() - self.ligand_center)[0]
            self.contact_unit_vector = contact_vector / np.linalg.norm(contact_vector)

            # Side chain has hydrogen bond donor/acceptor (DEHKNQRSTY)
            h_bond_donor_acceptor = 1 if self.residue_contact_atom.getResnames()[0] in ['ASP', 'GLU', 'HIS', 'LYS', 'ASN', 'GLN', 'ARG', 'SER', 'THR', 'TYR'] else 0

            # Residue characteristics
            # todo: UPDATE so that only one of the below can hold value of 1 at any given time
            # {0 | 1} - Hydrophobic, aliphatic(AILV)
            greasy_ali = 1 if self.residue_contact_atom.getResnames()[0] in ['ALA', 'ILE', 'LEU', 'VAL'] else 0
            # {0 | 1} - Hydrophobic, aromatic(FWY)
            greasy_aro = 1 if self.residue_contact_atom.getResnames()[0] in ['PHE', 'TYR', 'TRP'] else 0
            # {0 | 1} - Polar(NCQMST)
            polar = 1 if self.residue_contact_atom.getResnames()[0] in ['ASN', 'CYS', 'GLN', 'MET', 'SER', 'THR'] else 0
            # {0 | 1} - Charged, Acidic(DE)
            charged_acid = 1 if self.residue_contact_atom.getResnames()[0] in ['ASP', 'GLU'] else 0
            # {0 | 1} - Charged, Basic(HKR)
            charged_basic = 1 if self.residue_contact_atom.getResnames()[0] in ['HIS', 'LYS', 'ARG'] else 0
            # {0 | 1} - Glycine
            glycine = 1 if self.residue_contact_atom.getResnames()[0] in ['GLY'] else 0
            # {0 | 1} - Proline
            proline = 1 if self.residue_contact_atom.getResnames()[0] in ['PRO'] else 0
            # {0 | 1} - Backbone carbonyl
            bb_carbonyl = 1 if self.residue_contact_atom.getNames()[0] in ['O'] else 0
            # {0 | 1} - Backbone amino
            bb_amino = 1 if self.residue_contact_atom.getNames()[0] in ['N'] else 0
            # {0 | 1} - Backbone C / CA
            bb_c_ca = 1 if self.residue_contact_atom.getNames()[0] in ['C', 'CA'] else 0

            residue_vector = [
                self.contact_unit_vector[0] * self.weights[0],
                self.contact_unit_vector[1] * self.weights[0],
                self.contact_unit_vector[2] * self.weights[0],
                (self.min_contact_distance / self.distance_cutoff) * self.weights[0],
                residue_contact_type * self.weights[1],
                ligand_contact_type * self.weights[1],
                h_bond_donor_acceptor * self.weights[1],
                greasy_ali * self.weights[2],
                greasy_aro * self.weights[2],
                polar * self.weights[2],
                charged_acid * self.weights[2],
                charged_basic * self.weights[2],
                glycine * self.weights[2],
                proline * self.weights[2],
                bb_carbonyl * self.weights[3],
                bb_amino * self.weights[3],
                bb_c_ca * self.weights[3]
            ]

            return np.asarray(residue_vector)