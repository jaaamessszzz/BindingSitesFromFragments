#!/usr/bin/env python3

import numpy as np
import pypdb
import pubchempy

class Fragments():
    """
    This is a class for all things fragment related for generating hypothetical
    binding pockets.
    
    * Generate Fragments from Ligands
    * Query Pubchem or PDB for ligands containing defined fragments
    
    """
    def __init__(self):
        self.fragment_list = None

    def generate_fragements_from_ligand(self, ligand=None, ligand_input_format='name', method='smiles'):
        """
        Generate fragments from a given ligand.
        
        I don't really have a nice way of doing this yet... for the time being I'm going to define fragments by hand
        using the Fragment-Based Drug Discovery Rule of Threes:
        * Molecular weight <300Da
        * Up to three hydrogen bond donors/acceptors
        * LogP <= 3 (hydrophobicity)

        :param ligand: input ligand format defined be <ligand_input_format>
        :param ligand_input_format: format of <ligand>, small molecule name by default. [name|CID|smiles]
        :param method:
        :return:
        """

        compound = pubchempy.get_compounds([ligand], ligand_input_format)
        print(compound[0].canonical_smiles)

    def search_for_fragment_containing_ligands(self, predefined_fragments=False):
        """
        Search PDB or Pubchem for fragment-containing small molecules.
        
        For each fragment:
        * Identify all small molecules containing the fragment
        * Filter small molecules with representative fragments
            Check connectivities so that only break points in fragment generation are connected to a greater
            superstructure. Once this is defined in the fragment generation step, I should be able to check the 
            connectivities by difference: {bonds in molecule} - {bonds in fragment}, where the difference should
            be equivalent to the bonds defined as break points. 
        * Sort fragment-containing compounds 
        
        
        :param predefined_fragments: user-definied fragments? Assumes fragments generated using 
        Fragments.generate_fragements_from_ligand()
        :return: 
        """

