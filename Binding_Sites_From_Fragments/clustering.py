#!/usr/bin/env python3

import os
import numpy as np
import prody
import sklearn
import pprint
from .utils import *

class Cluster():
    """
    This class is responsible for taking aligned fragment PDBs and identifying representative contacts. 
    
    """

    def __init__(self, user_defined_dir, current_fragment, weights):
        self.processed_PDBs_dir = os.path.join(user_defined_dir, 'Transformed_Aligned_PDBs', current_fragment)
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
            prody_pdb = prody.parsePDB(pdb)
            processsed_residue_list += [fragment_PDB(pdb_info, self.weights, prody_pdb, residue) for residue in prody_pdb]
        return processsed_residue_list

    def cluster(self):
        """
        Using DBSCAN clustering algorithm to identify representitive contacts. I am using this algorithm because it
        takes density into consideration when forming clusters. The user can define what "dense" means
        
        :return: 
        """

        # for fragment in os.listdir(self.processed_PDBs_dir):
        #     frag_dir = os.path.join(self.processed_PDBs_dir, fragment)
        #     for pdb in pdb_check(frag_dir):
        #         self.process_structure_into_matrix(pdb)

        # first = np.asmatrix(
        #     [
        #         [1,0,0],
        #         [1,0,0],
        #         [1,0,0],
        #         [1,0,0]
        #     ]
        # )
        # second = np.asmatrix(
        #     [
        #         [0,0,1],
        #         [0,0,1],
        #         [1,0,0],
        #         [1,0,0]
        #                       ]
        #                      )
        #
        # print(np.linalg.norm(first - second))
        pprint.pprint(self.pdb_object_list)

class fragment_PDB():
    """
    Class for holding all information related to a processed fragment PDB
    
    :param pdb_info:
    :param prody_pdb:
    :param prody_residue:
    :param pdbid:
    :param vector:
    """
    def __init__(self, pdb_info, weights, prody_pdb, prody_residue):
        self.prody_pdb = prody_pdb
        self.prody_residue = prody_residue
        self.pdbid = pdb_info.split('_')[0]
        self.ligand_resname = pdb_info.split('_')[1]
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
         { 0 | 1 }   - Backbone contact  
         { 0 | 1 }   - Sidechain contact 

         { 0 | 1 }   - Ligand Polar Contact
         { 0 | 1 }   - Ligand Non-polar Contact
         
         { 0 | 1 }   - Hydrophobic, aliphatic { 0 | 1 }
         { 0 | 1 }   - Hydrophobic, aromatic  { 0 | 1 }
         { 0 | 1 }   - Polar                  { 0 | 1 }
         { 0 | 1 }   - Charged, Acidic        { 0 | 1 }
         { 0 | 1 }   - Charged, Basic         { 0 | 1 }
         { 0 | 1 }   - Glycine                { 0 | 1 }
         { 0 | 1 }   - Proline                { 0 | 1 }
         { 0 | 1 }   - Backbone carbonyl      { 0 | 1 }
         { 0 | 1 }   - Backbone amino         { 0 | 1 }
         { 0 | 1 }   - Backbone C/CA          { 0 | 1 }
        
        :return: 
        """
        # Contact distance
        prody_ligand = self.prody_pdb.select('resname {}'.format(self.ligand_resname))
        ligand_residue_distance_matrix = prody.buildDistMatrix(self.prody_residue, prody_ligand)

        # Find minimum score in matrix
        row_min_indicies = np.amin(ligand_residue_distance_matrix, axis=0)
        residue_index = np.argmin(row_min_indicies, axis=0)
        ligand_index = np.argmin(ligand_residue_distance_matrix, axis=0)

        column_index_low = residue_index
        row_index_low = ligand_index[column_index_low]

        print("Residue atom index: {}".format(column_index_low))
        print("Ligand atom index: {}".format(row_index_low))
        print("Contact distance: {}".format(ligand_residue_distance_matrix.item(row_index_low, column_index_low)))

        # Residue characteristics


        # Amino acid identity or degenerate amino acid groups?
        # Amino acid chemical characteristics e.g. [Sigma Amino Acid Reference Chart]
        #       (http://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html)


        # Vector from fragment centroid to {residue centroid | closest residue atom }


        # Backbone or side chain

        return [1]