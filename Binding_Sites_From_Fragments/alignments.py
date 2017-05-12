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
from rdkit import Chem
from rdkit.Chem import rdFMCS
import io


class Alignments():
    """
    This class is responsible for aligning all fragment-containing small molecule structures to the corresponding 
    fragments.
    
    Chris suggests looking into a package called RDKit for manipulating SMILES strings and structures
    
    """
    def __init__(self, user_defined_dir):
        self.user_defined_dir = user_defined_dir


    def fragment_target_mapping(self, fragment_pdb, target_string):
        """
        I'm going to utilize NetworkX (shoutout to Kale) for mapping fragment atoms to the corresponding atoms in each
        of the fragment-containing small molecules. [http://networkx.readthedocs.io]
        
        To do this, I will need atom and connectivity information for each small molecule to generate the nodes and 
        edges in the graphs. 
        
        OR. Or. I talked to Chris today and he showed me this chemoinformatics package called RDKit. This may be a much
        more streamlined method of doing this without needing to resort to subgraphs...
        > https://sourceforge.net/p/rdkit/mailman/message/34549645/
        
        :return: dict mapping fragment atoms to target atoms
        """
        # todo: make sure hydrogens aren't being stripped, this will mess up mapping
        # Import fragment and target ligand
        fragment_mol = Chem.MolFromPDBFile(fragment_pdb)
        target_mol = Chem.MolFromPDBBlock(target_string)

        # Generate a RDKit mol object of the common substructure so that I can map the fragment and target onto it
        substructure_match = rdFMCS.FindMCS([fragment_mol, target_mol])
        sub_mol = Chem.MolFromSmarts(substructure_match.smartsString)

        # Return atom indicies of substructure matches for fragment and target
        frag_matches = fragment_mol.GetSubstructMatch(sub_mol)
        target_matches = target_mol.GetSubstructMatch(sub_mol)

        # todo: figure out how to represent mapping efficiently
        # Maps fragment atom index to target atom index
        fragment_target_mapping = {f_idx: t_idx for f_idx, t_idx in zip(frag_matches, target_matches)}

        return fragment_target_mapping


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
        :return: String of HETAM and CONECT records from the input pdb for the given ligand
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

        ligand_io = ''.join(line_list)
        return ligand_io


    def clean_pdbs(self):
        """
        Clean up PDBs for fragment alignments.
        * Extract chain with target ligand bound
        
        :return:
        """
        pass


    def determine_rotation_and_translation(self, fragment_target_mapping, fragment_pdb, target_string):
        """
        Implementing the Kabsch algorithm for aligning all fragment-containing small molecules to the target ligand
        on the mapped atoms as determined by fragment_target_mapping()
        
        :return: Prody transformation object with rotation and translation matrix/vector
        """
        # Fragment and target as prody atom objects
        # Fortunately the atom indicies used in the atom mapping with RDKit match prody indicies
        # It's literally the line number of the HETATM records indexed starting at 0
        target_prody = prody.parsePDBStream(io.StringIO(target_string))
        fragment_prody = prody.parsePDB(fragment_pdb)

        frag_inx_string = [str(a) for a in fragment_target_mapping.keys()]
        trgt_inx_string = [str(a) for a in fragment_target_mapping.values()]

        target_align_atoms = target_prody.select('index ' + ' '.join(trgt_inx_string))
        fragment_align_atoms = fragment_prody.select('index ' + ' '.join(frag_inx_string))

        print(len(target_align_atoms))
        print(len(fragment_align_atoms))

        return prody.calcTransformation(target_align_atoms, fragment_align_atoms)