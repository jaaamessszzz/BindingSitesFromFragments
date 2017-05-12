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
import networkx
import prody
import biopandas

class Alignments():
    """
    This class is responsible for aligning all fragment-containing small molecule structures to the corresponding 
    fragments.
    
    Chris suggests looking into a package called RDKit for manipulating SMILES strings and structures
    
    """
    def __init__(self, user_defined_dir):
        self.user_defined_dir = user_defined_dir



    def fragment_target_mapping(self):
        """
        I'm going to utilize NetworkX (shoutout to Kale) for mapping fragment atoms to the corresponding atoms in each
        of the fragment-containing small molecules. [http://networkx.readthedocs.io]
        
        To do this, I will need atom and connectivity information for each small molecule to generate the nodes and 
        edges in the graphs. 
        
        OR. Or. I talked to Chris today and he showed me this chemoinformatics package called RDKit. This may be a much
        more streamlined method of doing this without needing to resort to subgraphs...
        > https://sourceforge.net/p/rdkit/mailman/message/34549645/
        
        :return: 
        """
        pass

    def generate_graph_from_pdb(self, pdb_file):
        """
        Generates a graph representation of a molecule using NetworkX
        
        :param pdb_file: Iterable containing information on molecule connectivity in PDB format. Iterable should either
        be a list or an opened instance of a .pdb file
        :return: 
        """
        pass


    def extract_atoms_and_connectivities(self, ligand, pdb_file):
        """
        Extract information on atoms and connectivities for a fragment-containing small molecule.
        
        Biopython is a pain in the ass
        BioPandas doesn't handle CONECT records
        Prody doesn't handle CONECT records
        
        :param ligand: Three letter code for fragment-containing ligand
        :param pdb_file: Path to the PDB file containing the fragment-containing ligand
        :return: list HETAM and CONECT records from the input pdb for the given ligand
        """

        # Determine which chains have my target ligand
        ag = prody.parsePDB(pdb_file)
        ag_ligands = ag.select('resname {}'.format(ligand))

        # Pick a chain and extract the target ligand... arbitrarily picks the first chain for now
        chain = ag_ligands.getChids()[0]

        # I only need CONECT and HETATM records for my ligand
        # This function assumes that all CONECT records will be at the very end of the PDB file
        pdb = open(pdb_file)
        line_list = []

        # Store HETAM lines for the ligand
        HETAM_num_list = []

        for line in pdb:
            split_line = line.split()
            if split_line[0] == 'HETATM' and split_line[3] == ligand and split_line[4] == chain:
                line_list.append(line)

                # Keep track of HETAM numbers so I know which CONET lines to pull
                HETAM_num_list.append(line.split()[1])

            # Save CONECT records if they reference one of the HETAMs in my ligand
            if line.split()[0] == 'CONECT':
                for element in line.split()[1:]:
                    if element in HETAM_num_list:
                        line_list.append(line)
                        break

        return line_list


    def clean_pdbs(self):
        """
        Clean up PDBs for fragment alignments.
        * Extract chain with target ligand bound
        
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