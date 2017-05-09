#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pypdb
import pubchempy
import xmltodict
import yaml
import pprint
import urllib
import os

class Alignments():
    """
    This class is responsible for aligning all fragment-containing small molecule structures to the corresponding 
    fragments.
    
    """
    def __init__(self, user_defined_dir):
        self.user_defined_dir = user_defined_dir



    def fragment_target_mapping(self):
        """
        I'm going to utilize NetworkX (shoutout to Kale) for mapping fragment atoms to the corresponding atoms in each
        of the fragment-containing small molecules. [http://networkx.readthedocs.io]
        
        To do this, I will need atom and connectivity information for each small molecule to generate the nodes and 
        edges in the graphs. 
        
        :return: 
        """
        pass


    def extract_atoms_and_connectivities(self, ligand, pdb_file):
        """
        Extract information on atoms and connectivities for a fragment-containing small molecule.
        
        :param ligand: Three letter code for fragment-containing ligand
        :param pdb_file: Path to the PDB file containing the fragment-containing ligand
        :return: 
        """

        print(ligand)
        pdb = open(pdb_file)
        for line in pdb:
            print(line)
            break


    def evaluate_and_clean_pdbs(self):
        """
        Determine whether PDBs are usable, and if so, clean them up for fragment alignments.
        
        :return: 
        """
        pass


    def determine_rotation_and_translation(self):
        """
        Implementing the Kabsch algorithm for aligning all fragment-containing small molecules to the target ligand
        on the mapped atoms as determined by fragment_target_mapping()
        
        :return: 
        """
        pass