#!/usr/bin/env python3

import os
import sys
import numpy as np
import prody
import sklearn
import pprint
from .utils import *

class Cluster():
    """
    This class is responsible for taking aligned fragment PDBs and identifying representative contacts. 
    """

    def __init__(self, processed_PDBs_dir, distance_cutoff, number_of_clusters, weights):
        self.processed_PDBs_dir = processed_PDBs_dir
        self.weights = weights
        self.number_of_clusters = number_of_clusters
        self.distance_cutoff = distance_cutoff
        self.pdb_object_list = self._import_pdbs()
        self.clusters = None

    def _import_pdbs(self):
        """
        For each fragment ensemble, converts each residue in all processed PDBs into objects with representative 
        matrices and other relevant information
        
        :return: 
        """
        processsed_residue_list = []
        for pdb in pdb_check(self.processed_PDBs_dir):
            pdb_info = os.path.basename(os.path.normpath(pdb))
            prody_protein = prody.parsePDB(pdb).select('protein')
            # Check that residues exist within cutoff distance provided in alignments, otherwise pass
            if prody_protein == None:
                continue
            else:
                prody_protein_hv = prody_protein.getHierView()

            prody_ligand = prody.parsePDB(pdb).select('resname {}'.format(pdb_info.split('_')[1]))

            # Iterate over residues in contacts and generate representative vector with weights applied
            for chain in prody_protein_hv:

                processsed_residue_list += [fragment_PDB(residue, prody_ligand, self.distance_cutoff, self.weights) for residue in chain]

        processsed_residue_list_cleaned = [residue for residue in processsed_residue_list if residue.vector is not None]

        return processsed_residue_list_cleaned

    def cluster_scipy(self):
        """
        Using Hierarchical clustering through scipy to identify representative contacts. 
        :return: 
        self.clusters: list of indicies corresponding to cluster numbers for each element in self.pdb_object_list
        """
        from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
        from matplotlib import pyplot as plt

        vector_array = np.asarray([thingy.vector for thingy in self.pdb_object_list])
        Z = linkage(vector_array, 'ward')

        # Generate dendrogram
        # todo: CHANGE THIS. Taken straight from https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial
        plt.figure(figsize=(25, 10))
        plt.title('Hierarchical Clustering Dendrogram')
        plt.xlabel('sample index')
        plt.ylabel('distance')
        dendrogram(
            Z,
            leaf_rotation=90.,  # rotates the x axis labels
            leaf_font_size=8.,  # font size for the x axis labels
        )
        plt.show()

        distance_cutoff = input('Enter desired distance cutoff: ').strip()
        self.clusters = fcluster(Z, distance_cutoff, criterion='distance')

    def cluster_sklearn_agglomerative(self):
        """
        Using Agglomerative Hierarchical clustering algorithm through sk-learn to identify representitive contacts.
        :return: 
        self.clusters: list of indicies corresponding to cluster numbers for each element in self.pdb_object_list
        """
        from sklearn.cluster import AgglomerativeClustering
        # todo: figure out a way to scale number of clusters based on the number of residues in the list...
        self.number_of_clusters = len(self.pdb_object_list)//20
        if self.number_of_clusters > 2:
            self.clusters = AgglomerativeClustering(n_clusters=self.number_of_clusters).fit_predict([thingy.vector for thingy in self.pdb_object_list])

    def generate_output_directories(self, base_dir, fragment_path):
        residue_cluster_pairs = [(cluster, residue)for cluster, residue in zip(self.clusters, self.pdb_object_list)]
        fragment = os.path.basename(os.path.normpath(fragment_path))

        for cluster in set(self.clusters):
            residue_list = [pair for pair in residue_cluster_pairs if pair[0] == cluster]
            fragment_cluster_path = os.path.join(base_dir,
                                     'Cluster_Results',
                                     fragment
                                     )
            os.makedirs(fragment_cluster_path,
                        exist_ok=True)
            count=0
            for residue in residue_list:
                # todo: add all residues to a single object and export that... the current onslaught of individual PDBs is way too intense
                prody.writePDB(os.path.join(fragment_cluster_path, 'Cluster_{0}_Resdiue_{1}.pdb'.format(cluster, str(count))), residue[1].prody_residue)
                count += 1

    def generate_reports(self, residue_list):
        """
        Generates reports on each cluster
        I guess the first step would be to basically output the information contained in each of the representative vectors?
        :return: 
        """
        # Create pandas dataframe, where each row has data on a specific residue
        # Aggregate and export .csv
        pass

class fragment_PDB():
    """
    Class for holding all information related to a processed fragment PDB
    
    :param pdb_info:
    :param prody_pdb:
    :param prody_residue:
    :param pdbid:
    :param vector:
    """
    def __init__(self, prody_residue, prody_ligand, distance_cutoff, weights):
        self.prody_residue = prody_residue
        self.prody_ligand = prody_ligand
        self.ligand_center = prody.calcCenter(prody_ligand)
        self.distance_cutoff = distance_cutoff
        self.weights = weights
        self.vector = self.process_residue_into_vector()

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
        ligand_residue_distance_matrix = prody.buildDistMatrix(self.prody_residue, self.prody_ligand)

        # Find minimum score in matrix
        row_min_indicies = np.amin(ligand_residue_distance_matrix, axis=0)
        ligand_index = np.argmin(row_min_indicies, axis=0)
        residue_index = np.argmin(ligand_residue_distance_matrix, axis=0)

        column_index_low = ligand_index
        row_index_low = residue_index[column_index_low]

        # Contact distance
        min_contact_distance = ligand_residue_distance_matrix.item(row_index_low, column_index_low)

        if min_contact_distance > self.distance_cutoff:
            return None

        else:
            residue_base_index = self.prody_residue.getIndices()[0]
            ligand_base_index = self.prody_ligand.getIndices()[0]

            residue_contact_atom = self.prody_residue.select('index {}'.format(residue_base_index + row_index_low))
            ligand_contact_atom = self.prody_ligand.select('index {}'.format(ligand_base_index + column_index_low))

            # Residue Contact Type
            residue_contact_type = 0 if residue_contact_atom.getNames()[0] in ['C', 'CA', 'N', 'O'] else 1

            # Ligand Contact Type
            ligand_contact_type = 1 if ligand_contact_atom.getNames()[0][0] in ['C'] else 0

            # Vector from fragment centroid to closest residue atom
            vector = ((residue_contact_atom.getCoords() - self.ligand_center)[0]) / min_contact_distance

            # Side chain has hydrogen bond donor/acceptor (DEHKNQRSTY)
            h_bond_donor_acceptor = 1 if residue_contact_atom.getResnames()[0] in ['ASP', 'GLU', 'HIS', 'LYS', 'ASN', 'GLN', 'ARG', 'SER', 'THR', 'TYR'] else 0

            # Residue characteristics
            # {0 | 1} - Hydrophobic, aliphatic(AILV)
            greasy_ali = 1 if residue_contact_atom.getResnames()[0] in ['ALA', 'ILE', 'LEU', 'VAL'] else 0
            # {0 | 1} - Hydrophobic, aromatic(FWY)
            greasy_aro = 1 if residue_contact_atom.getResnames()[0] in ['PHE', 'TYR', 'TRP'] else 0
            # {0 | 1} - Polar(NCQMST)
            polar = 1 if residue_contact_atom.getResnames()[0] in ['ASN', 'CYS', 'GLN', 'MET', 'SER', 'THR'] else 0
            # {0 | 1} - Charged, Acidic(DE)
            charged_acid = 1 if residue_contact_atom.getResnames()[0] in ['ASP', 'GLU'] else 0
            # {0 | 1} - Charged, Basic(HKR)
            charged_basic = 1 if residue_contact_atom.getResnames()[0] in ['HIS', 'LYS', 'ARG'] else 0
            # {0 | 1} - Glycine
            glycine = 1 if residue_contact_atom.getResnames()[0] in ['GLY'] else 0
            # {0 | 1} - Proline
            proline = 1 if residue_contact_atom.getResnames()[0] in ['PRO'] else 0
            # {0 | 1} - Backbone carbonyl
            bb_carbonyl = 1 if residue_contact_atom.getNames()[0] in ['O'] else 0
            # {0 | 1} - Backbone amino
            bb_amino = 1 if residue_contact_atom.getNames()[0] in ['N'] else 0
            # {0 | 1} - Backbone C / CA
            bb_c_ca = 1 if residue_contact_atom.getNames()[0] in ['C', 'CA'] else 0

            residue_vector = [
                vector[0] * self.weights[0],
                vector[1] * self.weights[0],
                vector[2] * self.weights[0],
                (min_contact_distance / self.distance_cutoff) * self.weights[0],
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