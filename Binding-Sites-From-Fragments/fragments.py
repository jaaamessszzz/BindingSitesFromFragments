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
        
        todo:
        * File input (probably yaml) to define fragments and breaks in target ligand
        
        :param ligand: input ligand format defined be <ligand_input_format>
        :param ligand_input_format: format of <ligand>, small molecule name by default. [name|CID|smiles]
        :param method:
        :return:
        """
        import pprint

        compound = pubchempy.get_compounds([ligand], ligand_input_format)
        pprint.pprint(compound[0].to_dict())

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
            
            EDIT: LOL this is not necessary if you explicitly specify hydrogens in the substructure search.
            Unfortunately, the substructure search function on the PDB sucks. Pubchem works quite well, but now I need
            a method of mapping fragment-containing small molecules to molecules present in the PDB.
            
            http://ligand-expo.rcsb.org/ld-download.html
            > InChIKey (tab delimited text)
            > Tabulation of PDB entries containing each chemical component. (tab delimited text)
            
            The InChIKey tsv maps all PDB ligand three-letter codes to the InChIKey. I can do a set intersection of 
            a Pubchem substructure search to this list and get all PDB codes for structures with my fragment-containing
            small molecules bound.
            
        * Sort fragment-containing compounds 
        
        todo:
        * Map fragment atoms to matched ligands
            I think the best way to do this would be to define three atoms in each fragment to calculate the 
            rotation/translation matrices. This would ensure all structures can be aligned with minimal amount of 
            user input and calculation time. Unfortunately, this will still require user input for mapping fragment 
            atoms to each fragment-containing small molecule...
        
        :param predefined_fragments: user-definied fragments? Assumes fragments generated using 
        Fragments.generate_fragements_from_ligand()
        :return: 
        """
        pass

    def bootstrap_method(self):
        """
        Doing things by hand since I still need to define rigorous protocols for the above methods.
        
        :return: 
        """
        import xmltodict
        import yaml
        import pprint
        import urllib

        # Substructure search for each fragment, returning ligands
        search_dict = yaml.load_all(open("./Working_Files/Search_query.yaml"))
        # YAML file is split up into individual queries for each fragment
        # Ideally I would just use one skeleton query yaml, but I can't think of a good way of substituting in the
        # smiles string into the dict. The dicts are structured improperly as they're being converted into XMLs, which
        # means I can't manipulate them normally before I execute searches.

        for search_query in search_dict:
            search_xml = xmltodict.unparse(search_query, pretty=False)
            search_request = urllib.request.Request('http://www.rcsb.org/pdb/rest/search', data=search_xml.encode())
            search_result = urllib.request.urlopen(search_request).read()

            pprint.pprint(str(search_result))



