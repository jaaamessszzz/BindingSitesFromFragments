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

    def __init__(self, user_defined_dir, current_fragment):
        self.processed_PDBs_dir = os.path.join(self.user_defined_dir, 'Transformed_Aligned_PDBs', current_fragment)
        self.pdb_object_list = self.import_pdbs()

    def import_pdbs(self):
        """
        For each fragment ensemble, converts each residue in all processed PDBs into objects with representative 
        matrices and other relevant information
        
        :return: 
        """
        processsed_residue_list = []
        for pdb in pdb_check(self.processed_PDBs_dir):
            processsed_residue_list += [fragment_PDB(pdb, residue) for residue in prody.parsePDB(pdb)]
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


class fragment_PDB():
    """
    Class for holding all information related to a processed fragment PDB
    """
    def __init__(self, prody_pdb, prody_residue):
        self.prody_pdb = prody_pdb
        self.prody_residue = prody_residue
        self.pdbid = None
        self.matrix = self.process_residue_into_matrix()

    def process_residue_into_matrix(self):
        """
        Converting each residue into a representative matrix

        Metrics I care about for each residue:
        * Contact distance
        * Residue characteristics 
            * Amino acid identity or degenerate amino acid groups?
            * Amino acid chemical characteristics e.g. [Sigma Amino Acid Reference Chart]
              (http://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html)
        * Position of residue relative to fragment
            * Vector from fragment centroid to {residue centroid | closest residue atom }
        * Backbone or side chain

        :param pdb: path to processed PDB
        :return: 
        """
        # Contact distance


        # Residue characteristics


        # Amino acid identity or degenerate amino acid groups?
        # Amino acid chemical characteristics e.g. [Sigma Amino Acid Reference Chart]
        #       (http://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html)


        # Vector from fragment centroid to {residue centroid | closest residue atom }


        # Backbone or side chain
