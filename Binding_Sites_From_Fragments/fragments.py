#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pypdb
import pubchempy
from Bio import PDB
import xmltodict
import pprint
import urllib
import os

class Fragments():
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

        fragment_inputs = pd.read_csv("{}/Inputs/Fragment_inputs.csv".format(self.user_defined_dir))
        search_dir = "{}/Inputs".format(self.user_defined_dir)

        for index, search_query in fragment_inputs.iterrows():
            # Using Pubchem API (if I could get it to work...)
            # pubchem_search_results = pubchempy.get_compounds("N1C=CN=C1",
            #                                                  searchtype='substructure',
            #                                                  listkey_count=0,
            #                                                  as_dataframe=True)

            current_fragment = "Fragment_{}".format(search_query["Fragment"])

            # Make directory for each fragment
            os.makedirs(os.path.join(self.user_defined_dir,
                                     "Fragment_PDB_Matches",
                                     current_fragment),
                        exist_ok=True)
            # print(os.path.join("Fragment_{}".format(str(search_query["Fragment"])), fragment_smiles))

            pubchem_results = pd.read_csv(os.path.join(search_dir, "{}.csv".format(current_fragment)))

            ligands_in_pdb_df = pubchem_results.merge(InChIKey_to_PDB, how='left', on='cmpdinchikey')
            ligands_in_pdb_df = ligands_in_pdb_df.dropna(how='any')

            # todo: Add option to export PBD ligand matches
            ligands_in_pdb_df.to_csv(os.path.join(self.user_defined_dir, "Fragment_PDB_Matches", current_fragment, '{}_pdb.csv'.format(current_fragment)))

            # Generate set of PDBIDs for fragment-containing ligands
            pdbid_set = set(ligands_in_pdb_df['PDBID'])
            print(pdbid_set)
            for ligand_pdbid in pdbid_set:
                # Make directory to dump PDB files for each fragment-containing ligand
                fragment_ligand_path = os.path.join(self.user_defined_dir,
                                                    "Fragment_PDB_Matches",
                                                    current_fragment,
                                                    ligand_pdbid)

                # Download PDBs
                self.download_PDBs(ligand_pdbid, fragment_ligand_path)

    def download_PDBs(self, ligand_pdbid, fragment_ligand_path):
        """
        Download all PDBs with the provided fragment-containing ligand
        :param ligand_pdbid: three character PDBID for the fragment-containing ligand
        :return: 
        """

        REST_search_xml = """
        <orgPdbCompositeQuery version="1.0">
        <queryRefinement>
        <queryRefinementLevel>0</queryRefinementLevel>
        <orgPdbQuery>
        <queryType>org.pdb.query.simple.ChemCompIdQuery</queryType>
        <chemCompId>{}</chemCompId>
        <polymericType>Free</polymericType>
        </orgPdbQuery>
        </queryRefinement>
        <queryRefinement>
        <queryRefinementLevel>1</queryRefinementLevel>
        <conjunctionType>and</conjunctionType>
        <orgPdbQuery>
        <queryType>org.pdb.query.simple.HomologueReductionQuery</queryType>
        <identityCutoff>70</identityCutoff>
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
        """.format(ligand_pdbid)

        search_xml = xmltodict.unparse(xmltodict.parse(REST_search_xml), pretty=False)
        search_request = urllib.request.Request('http://www.rcsb.org/pdb/rest/search', data=search_xml.encode())
        search_result = urllib.request.urlopen(search_request, data=None, timeout=300).read()

        result_pdbs = search_result.split()

        if len(result_pdbs) > 0:
            os.makedirs(fragment_ligand_path, exist_ok=True)
            print(result_pdbs)
            print("Created:", fragment_ligand_path)

            # Download PDB files
            for pdb in result_pdbs:
                pdb_name = pdb.decode("utf-8")

                if not os.path.exists(os.path.join(fragment_ligand_path, "{}.pdb".format(pdb_name))):
                    try:
                        print("Downloading {}...".format(pdb_name))
                        pdb_string = pypdb.get_pdb_file(pdb_name, filetype='pdb')
                        current_download = open(os.path.join(fragment_ligand_path, "{}.pdb".format(pdb_name)), 'w')
                        current_download.write(pdb_string)
                        current_download.close()
                    except Exception as e:
                        print(e)
                else:
                    print(os.path.join(fragment_ligand_path, "{}.pdb".format(pdb_name)), "exists! Moving on...")







