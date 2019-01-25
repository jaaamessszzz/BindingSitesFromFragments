#!/usr/bin/env python3

import pandas as pd
import pypdb
import json
import pubchempy
import xmltodict
import pprint
import urllib
import os

from .utils import *

# --- Silence ProDy --- #
prody.confProDy(verbosity='none')

class Fragments(object):
    """
    This is a class for all things fragment related for generating hypothetical
    binding pockets.
    
    * Generate Fragments from Ligands
    * Query Pubchem or PDB for ligands containing defined fragments
    
    * Directory Tree
        
    user_defined_dir
    |--Inputs
    |  |--Fragment_inputs.csv
    |  |--[Fragment 1].csv
    |  |--[Fragment 2].csv
    |  |--[Fragment 1].mol # TBD
    |  |--[Fragment 2].mol # TBD
    |
    |--Fragment_PDB_Matches
       |--[Fragment smiles 1]
       |  |--[Fragment smiles 1]_pdb.csv
       |  |--[Ligand PDBID 1]
       |  |  |--[PDB 1]
       |  |  |--[PDB 2]
       |  |  |--[PDB ...]
       |  |
       |  |--[Ligand PDBID 2]
       |     |--[PDB 1]
       |     |--[PDB 2]
       |     |--[PDB ...]
       |
       |--[Fragment smiles 2]
       |  |--[Fragment smiles 2]_pdb.csv
       |  |--[Ligand PDBID]
       |  |  |--[PDB 1]
       |  |  |--[PDB 2]
       |  |  |--[PDB ...]
       |  |
       |  |--[Ligand PDBID 2]
       |     |--[PDB 1]
       |     |--[PDB 2]
       |     |--[PDB ...]
                  
    """
    def __init__(self, user_defined_dir):
        self.user_defined_dir = user_defined_dir
        self.resources_dir = os.path.join(os.path.dirname(__file__), '..', 'Additional_Files')

    # NEVER IMPLEMENTED!
    def generate_fragements_from_ligand(self, ligand=None, ligand_input_format='name', method='smiles'):
        """
        Generate fragments from a given ligand.

        > Look into RDKit.Chem.Lipinski.NumHAcceptors(x) + RDKit.Chem.Lipinski.NumHDonors(x) if I want to systematically
        generate fragments.

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
            Unfortunately, the substructure search function on the PDB sucks. Pubchem works quite well. I have found
            that it is important to manually set the number of allowed neighbors/substituents for atoms so you don't
            get branching from the specified atoms. This is done with the query atom tool in the PubChem Sketcher.
            It seems that the all of this additional information is encoded in the SMILES string.
            
            https://pubchem.ncbi.nlm.nih.gov/sketch/sketchhelp.html#wp915615
            
            Now I need a method of mapping fragment-containing small molecules to molecules present in the PDB.
            
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
            
            Actually I just thought of something. I know Kale may not like this idea, but if I limit the user to
            defining fragments without any rotational DOFs, I should be able to define three points and the distances
            between each. This should be constant for all fragment-containing small molecule matches and can be applied
            to automatically aligning all matches. 

        :param predefined_fragments: user-definied fragments? Assumes fragments generated using 
        Fragments.generate_fragements_from_ligand()
        :return: 
        """

        InChIKey_to_PDB = pd.read_csv(os.path.join(self.resources_dir, "Components-inchikey.tsv"),
                                      delimiter="\t",
                                      names=["cmpdinchikey", "PDBID", "Name"])

        # Perform fragment search through PubChem
        # I wanted to do this through the pubchempy API, but I can't get it to work for some reason
        # Falling back to manually doing the searches on Pubchem and downloading the detailed search results
        # Name each search using the SMILES string to simplify things

        search_dir = os.path.join(f"{self.user_defined_dir}", "Inputs", "Fragment_Inputs")
        fragment_inputs = pd.read_csv(os.path.join(search_dir, "Fragment_inputs.csv"))

        pdb_dict = {}

        for index, search_query in fragment_inputs.iterrows():
            # todo: Figure out how the PubChem PUG REST API works...
            current_fragment = "Fragment_{}".format(search_query["Fragment"])
            pdb_dict[current_fragment] = {}

            # Make directory for each fragment
            os.makedirs(os.path.join(self.user_defined_dir, "Fragment_PDB_Matches", current_fragment), exist_ok=True)
            pubchem_results = pd.read_csv(os.path.join(search_dir, "{}.csv".format(current_fragment)))
            ligands_in_pdb_df = pubchem_results.merge(InChIKey_to_PDB, how='inner', on='cmpdinchikey')

            ligands_in_pdb_df.to_csv(os.path.join(self.user_defined_dir, "Fragment_PDB_Matches", current_fragment, '{}_pdb.csv'.format(current_fragment)))

            # Generate set of PDBIDs for fragment-containing ligands
            pdbid_set = set(ligands_in_pdb_df['PDBID'])

            # REPORT
            print(f'\n{current_fragment}\n')
            print(f'Ligands found to contain {current_fragment}:\n{pdbid_set}')

            # Download PDBs
            pdb_dict[current_fragment]['Ligands'] = list(pdbid_set)
            pdb_dict[current_fragment]['PDBs'] = self.search_PDBs(pdbid_set)

        with open(os.path.join(self.user_defined_dir, 'PDB_search_results.json'), 'w') as jsonfile:
            json.dump(pdb_dict, jsonfile)

    def search_PDBs(self, ligand_pdbid_list):
        """
        Download all PDBs with the provided fragment-containing ligand
        :param ligand_pdbid: three character PDBID for the fragment-containing ligand
        :return: 
        """

        REST_search_xml = f"""
        <orgPdbCompositeQuery version="1.0">
        <queryRefinement>
        <queryRefinementLevel>0</queryRefinementLevel>
        <orgPdbQuery>
        <queryType>org.pdb.query.simple.ChemCompIdQuery</queryType>
        <chemCompId>{" ".join(ligand_pdbid_list)}</chemCompId>
        <polymericType>Free</polymericType>
        </orgPdbQuery>
        </queryRefinement>
        <queryRefinement>
        <queryRefinementLevel>1</queryRefinementLevel>
        <conjunctionType>and</conjunctionType>
        <orgPdbQuery>
        <queryType>org.pdb.query.simple.HomologueReductionQuery</queryType>
        <identityCutoff>90</identityCutoff>
        </orgPdbQuery>
        </queryRefinement>
        <queryRefinement>
        <queryRefinementLevel>2</queryRefinementLevel>
        <conjunctionType>and</conjunctionType>
        <orgPdbQuery>
        <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
        <containsProtein>Y</containsProtein>
        <containsDna>N</containsDna>
        <containsRna>N</containsRna>
        <containsHybrid>N</containsHybrid>
        </orgPdbQuery>
        </queryRefinement>
        </orgPdbCompositeQuery>
        """

        search_xml = xmltodict.unparse(xmltodict.parse(REST_search_xml), pretty=False)
        search_request = urllib.request.Request('http://www.rcsb.org/pdb/rest/search', data=search_xml.encode())
        search_result = urllib.request.urlopen(search_request, data=None, timeout=10).read()

        # REPORT
        print('\nPDB records containing relevant ligands:', [pdb.decode('utf-8') for pdb in search_result.split()])

        return [pdb.decode('utf-8').lower() for pdb in search_result.split()]