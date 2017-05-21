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

    def __init__(self, processed_PDBs_dir, weights):
        self.processed_PDBs_dir = processed_PDBs_dir
        self.weights = weights
        self.pdb_object_list = self._import_pdbs()

    def _import_pdbs(self):
        """
        For each fragment ensemble, converts each residue in all processed PDBs into objects with representative 
        matrices and other relevant information
        
        :return: 
        """
        processsed_residue_list = []
        for pdb in pdb_check(self.processed_PDBs_dir):
            pdb_info = os.path.basename(os.path.normpath(pdb))
            prody_protein_hv = prody.parsePDB(pdb).select('protein').getHierView()
            prody_ligand = prody.parsePDB(pdb).select('resname {}'.format(pdb_info.split('_')[1]))

            # Iterate over residues in contacts and generate representative vector with weights applied
            for chain in prody_protein_hv:
                processsed_residue_list += [fragment_PDB(residue, prody_ligand, self.weights) for residue in chain]

        return processsed_residue_list

    def cluster(self):
        """
        Using DBSCAN clustering algorithm to identify representitive contacts. I am using this algorithm because it
        takes density into consideration when forming clusters. The user can define what "dense" means
        
        :return: 
        """
        from sklearn.cluster import DBSCAN
        from sklearn.cluster import AgglomerativeClustering
        from sklearn import metrics

        # asdf = DBSCAN().fit_predict([thingy.vector for thingy in self.pdb_object_list])
        asdf = AgglomerativeClustering(n_clusters=6).fit_predict([thingy.vector for thingy in self.pdb_object_list])

        print(asdf)

class fragment_PDB():
    """
    Class for holding all information related to a processed fragment PDB
    
    :param pdb_info:
    :param prody_pdb:
    :param prody_residue:
    :param pdbid:
    :param vector:
    """
    def __init__(self, prody_residue, prody_ligand, weights):
        self.prody_residue = prody_residue
        self.prody_ligand = prody_ligand
        self.ligand_center = prody.calcCenter(prody_ligand)
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
        
         {Angstroms} - X component, Vector from fragment centroid to closest residue atom 
         {Angstroms} - Y component, Vector from fragment centroid to closest residue atom 
         {Angstroms} - Z component, Vector from fragment centroid to closest residue atom
         
         {Angstroms} - Contact distance  
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

        residue_base_index = self.prody_residue.getIndices()[0]
        ligand_base_index = self.prody_ligand.getIndices()[0]

        residue_contact_atom = self.prody_residue.select('index {}'.format(residue_base_index + row_index_low))
        ligand_contact_atom = self.prody_ligand.select('index {}'.format(ligand_base_index + column_index_low))

        # Residue Contact Type
        residue_contact_type = 0 if residue_contact_atom.getNames()[0] in ['C', 'CA', 'N', 'O'] else 1

        # Ligand Contact Type
        ligand_contact_type = 1 if ligand_contact_atom.getNames()[0][0] in ['C'] else 0

        # Vector from fragment centroid to closest residue atom
        vector = (residue_contact_atom.getCoords() - self.ligand_center)[0]

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
            min_contact_distance * self.weights[0],
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