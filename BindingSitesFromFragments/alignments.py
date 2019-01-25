#!/usr/bin/env python3

import os
import io
import re
import sys
import json
from pprint import pprint

import prody
import numpy as np
import pandas as pd

import urllib
import requests
import xmltodict

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdFMCS
from rdkit.Chem.Fingerprints import FingerprintMols

from .utils import *

# --- Silence ProDy --- #
prody.confProDy(verbosity='none')

class Align_PDB_Factory(object):
    """
    This class does all the gruntwork for finding all the information required to process a PDB bound to a fragment-
    containing compound. All of the organized information is passed to Align_PDB.
    """

    def __init__(self, user_defined_dir):

        # --- Paths --- #
        self.user_defined_dir = user_defined_dir
        self.pdb_bank_dir = os.path.join(self.user_defined_dir, 'PDB_Bank')
        self.fragment_inputs_path = os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs')
        self.rejected_list_path = None  # Assigned in load_previously_rejected_pdbs()

        # --- Load data from files --- #
        self.ligands_to_exclude = self._init_ligands_to_exclude()
        self.pdb_ligand_json = self._init_pubchem_search_results()
        self.fragment_df = pd.read_csv(os.path.join(self.fragment_inputs_path, 'Fragment_inputs.csv'))

        # --- Other variables --- #
        self.sanitized_smiles_dict = {f"Fragment_{row['Fragment']}": re.sub(r'[^\[\]]+(?=\])', lambda x: f"{x.group().split(';')[0]}",
                                                                            row['SMILES_fragment']) for index, row in self.fragment_df.iterrows()}

    def alignment_monstrosity(self, rmsd_cutoff=0.5, use_cluster_pdb_database=False):
        """
        Consequences of not thinking ahead...
        For each fragment, align all fragment-containing ligands to fragment
        Generate PDBs with aligned coordinate systems
        :param args:
        :param rmsd_cutoff: fragment alignment RMSD cutoff, anything higher gets rejected
        :return:
        """
        # If use_cluster_pdb_database=False, use PDB FTP to download all structures
        # Otherwise, all relevant structures should be found in the local PDB database on Netapp
        if not use_cluster_pdb_database:
            prody.pathPDBFolder(folder=self.pdb_bank_dir)

            for current_fragment in self.pdb_ligand_json:

                # Only download PDBs that aren't already in PDB bank directory
                existing_PDBs = [pdb[:4].lower() for pdb in os.listdir(self.pdb_bank_dir)]
                PDBs_to_download = list(set(self.pdb_ligand_json[current_fragment]['PDBs']) - set(existing_PDBs))

                if len(PDBs_to_download) > 0:
                    print(f'Downloading PDBs for {current_fragment}...\n')
                    prody.fetchPDBviaFTP(*PDBs_to_download)
                else:
                    print(f'All relevant PDBs for {current_fragment} found in {self.pdb_bank_dir}!\n')

        # Create directory for processed PDBs
        processed_PDBs_path = os.path.join(self.user_defined_dir, 'Transformed_Aligned_PDBs')
        rejected_list = self.load_previously_rejected_pdbs(processed_PDBs_path)

        # Create directories...
        os.makedirs(self.pdb_bank_dir, exist_ok=True)
        os.makedirs(processed_PDBs_path, exist_ok=True)

        # Fragment_1, Fragment_2, ...
        for current_fragment in self.pdb_ligand_json:

            # Create directory for processed PDBs
            processed_dir = os.path.join(processed_PDBs_path, current_fragment)
            os.makedirs(processed_dir, exist_ok=True)

            # Save ideal_ligand_containers for each fragment so things are only downloaded once
            ideal_ligand_dict = dict()
            ideal_ligand_dict['Ligands'] = dict()
            ideal_ligand_dict['Failed'] = list()

            # Align_PDB class holds all information for the current fragment
            align = Align_PDB(self.user_defined_dir, current_fragment, self.sanitized_smiles_dict[current_fragment])

            # For each PDB containing a fragment-containing compound
            for pdbid in self.pdb_ligand_json[current_fragment]['PDBs']:

                # Return path of PDB file to use for processing
                found_pdb, pdb_path = self.return_PDB_to_use_for_alignments(pdbid, self.pdb_bank_dir, use_cluster_pdb_database=False)

                if not found_pdb:
                    continue

                # Proceed with processing if the current PDB passes all filters
                if not processed_check(processed_dir, pdbid, rejected_list):
                    print("\n\nProcessing {}...".format(pdbid))

                    # --- Check which ligands contain relevant fragments --- #

                    relevant_ligands = self.return_substructure_containing_ligands(pdb_path, self.pdb_ligand_json, current_fragment)

                    # Set things up! Get ligands from Ligand Expo if haven't already tried and failed
                    for ligand in relevant_ligands:

                        if not ideal_ligand_dict['Ligands'].get(ligand) and ligand not in ideal_ligand_dict['Failed']:
                            ideal_ligand_container = Ideal_Ligand_PDB_Container(ligand)

                            if ideal_ligand_container.success:
                                ideal_ligand_dict['Ligands'][ligand] = ideal_ligand_container
                            else:
                                ideal_ligand_dict['Failed'].append(ligand)

                    # Create a temp list for ligands that will be pulled from the current PDB
                    ligand_container_dict_for_current_pdb = {lig: ideal_ligand_dict['Ligands'][lig] for lig in ideal_ligand_dict['Ligands'] if lig in relevant_ligands}
                    relevant_ligands_prody_dict = align.extract_ligand_records(pdb_path, ligand_container_dict_for_current_pdb)

                    # Reject if no ligands with all atoms represented can be found for the given PDB
                    # todo: only reject ligands where fragment atoms are missing...
                    if len(relevant_ligands_prody_dict) < 1:
                        with open(self.rejected_list_path, 'a+') as reject_list:
                            reject_list.write('{}\n'.format(pdbid))
                        print('REJECTED - no target ligands were fully represented in the PDB')
                        continue

                    # --- Perform alignment of PDB fragment substructure (mobile) onto defined fragment (target) --- #

                    # ...if PDB has not been processed, rejected, or excluded by the user

                    else:

                        # Iterate over ligands found to contain fragments as substructures
                        for ligand_resname, ligand_chain, ligand_resnum in relevant_ligands_prody_dict:

                            # Mapping of fragment atoms to target ligand atoms
                            target_ligand_ideal_smiles = ligand_container_dict_for_current_pdb[ligand_resname].smiles
                            
                            target_ligand_pdb_string = io.StringIO()
                            target_ligand_prody = relevant_ligands_prody_dict[(ligand_resname, ligand_chain, ligand_resnum)].select('not hydrogen')
                            prody.writePDBStream(target_ligand_pdb_string, target_ligand_prody)

                            mapping_successful, fragment_target_map = align.fragment_target_mapping(target_ligand_ideal_smiles, target_ligand_pdb_string)

                            if not mapping_successful:
                                with open(self.rejected_list_path, 'a+') as reject_list:
                                    reject_list.write('{}\n'.format(pdbid))
                                print('REJECTED - failed atom mapping between target and reference fragment')
                                continue

                            print(f'\n{len(fragment_target_map)} possible mapping(s) of fragment onto {pdbid}:{ligand} found...\n')

                            # Iterate over possible mappings of fragment onto current ligand
                            for count, mapping in enumerate(fragment_target_map):

                                # Determine translation vector and rotation matrix
                                target_coords_and_serials, frag_atom_coords, transformation_matrix = align.determine_rotation_and_translation(mapping, target_ligand_prody)
                                trgt_atom_coords, target_fragment_atom_serials = target_coords_and_serials

                                # Apply transformation to protein_ligand complex if rmsd if below cutoff
                                # todo: Use information from PubChem fragment SMILES in determining correct mappings
                                rmsd = prody.calcRMSD(frag_atom_coords, prody.applyTransformation(transformation_matrix, trgt_atom_coords))
                                print('RMSD of target onto reference fragment:\t{}'.format(rmsd))

                                if rmsd < rmsd_cutoff:
                                    transformed_pdb = align.apply_transformation(pdb_path, ligand_resnum, target_fragment_atom_serials, transformation_matrix)
                                    transformed_pdb_name = f'{pdbid}_{ligand_resname}_{ligand_chain}_{ligand_resnum}-{count}.pdb'
                                    prody.writePDB(os.path.join(processed_dir, transformed_pdb_name), transformed_pdb)

                                else:
                                    with open(self.rejected_list_path, 'a+') as reject_list:
                                        reject_list.write('{}\n'.format(pdbid))
                                    print('REJECTED - high RMSD upon alignment to reference fragment')

                else:
                    print('{} has been processed!'.format(pdbid))

    # --- New methods after refactoring clusterfuck --- #

    def _init_ligands_to_exclude(self):
        """
        Parses the user-defined Exclude_Ligands.txt if it exists

        :return: List of three-letter compound identifier codes to exclude if Exclude_Ligands.txt exists, else None
        """
        exclude_txt = os.path.join(self.user_defined_dir, 'Inputs', 'User_Inputs', 'Exclude_Ligands.txt')
        if os.path.exists(exclude_txt):
            with open(exclude_txt, 'r') as exlude_ligands:
                return [lig.strip() for lig in exlude_ligands]
        else:
            return None

    def _init_pubchem_search_results(self):
        """
        Loads PubChem search results from PDB_search_results.json

        :return: Dict containing PubChem search results from PDB_search_results.json
        """
        with open(os.path.join(self.user_defined_dir, 'PDB_search_results.json'), 'r') as jsonfile:
            return json.load(jsonfile)

    def load_previously_rejected_pdbs(self, processed_PDBs_path):
        """
        Loads list of previously rejected pdb files

        :return:
        """
        self.rejected_list_path = os.path.join(processed_PDBs_path, 'Rejected_PDBs.txt')

        if os.path.exists(self.rejected_list_path):
            with open(self.rejected_list_path, 'r') as rejected_PDBs:
                return [pdb.strip() for pdb in rejected_PDBs]

        else:
            return list()

    def return_PDB_to_use_for_alignments(self, pdbid, pdb_bank_dir, use_cluster_pdb_database=False):
        """
        Returns the path of the PDB to use for alignment processing.

        :param pdbid: 4-character PDB identifier code
        :param pdb_bank_dir: directory to dump pdb.gz files
        :param use_cluster_pdb_database: if True, use the local PDB database on the cluster
        :return:
        """
        # --- Download PDB via RCSB FTP or find locally --- #
        if f'{pdbid}.pdb.gz' in os.listdir(pdb_bank_dir):

            pdb_path = os.path.join(pdb_bank_dir, f'{pdbid}.pdb.gz')
            return True, pdb_path

        elif use_cluster_pdb_database:

            pdb_path = f'/netapp/database/pdb/remediated/pdb/{pdbid[1:-1]}/pdb{pdbid}.ent.gz'
            pdb_path_exists = os.path.exists(pdb_path)

            print(f'Path on netapp: {pdb_path}\nExists: {pdb_path_exists}')

            if pdb_path_exists:
                return True, pdb_path
            else:
                return False, None

        else:

            try:
                pdb_path = prody.fetchPDBviaFTP(pdbid)
                return True, [pdb_path] if len(pdb_path) == 1 else pdb_path

            except Exception as e:
                print(f'\nThere was an error retrieving {pdbid.upper()} from the RCSB FTP server:\n{e}')
                return False, None

    def return_substructure_containing_ligands(self, pdb_path, pdb_ligand_json, current_fragment):
        """
        Parses the PDB header to find ligands that were found to contain the current fragment as a substructure
        through the PubChem substructure search

        :param pdb_path: path to PDB to use for alignment processing
        :param pdb_ligand_json: parsed JSON->dict containing PDB-Ligand search results from PubChem
        :param current_fragment: String for current fragment e.g. 'Fragment_1'

        :return: list of ligands found in pdb that containing fragment substructure(s)
        """
        # Use PDB header to identify ligands in pdb_ligand_json[current_fragment]['Ligands']
        pdb_header = prody.parsePDB(pdb_path, header=True, model=0)
        pdb_ligand_set = set([lig.resname for lig in pdb_header['chemicals']])

        return list(pdb_ligand_set & set(pdb_ligand_json[current_fragment]['Ligands']))


class Align_PDB(object):
    """
    Base class with functions for aligning individual PDBs
    NOTE: FRAGMENT NEEDS TO BE SUBSTRUCTURE OF TARGET FOR PROPER MAPPING!!!

    :param pdb_file: path to input PDB containing protein bound to fragment-containing target molecule (MOBILE)
    :param target_string: string of PDB for target ligand, fetched from LigandExpo (MOBILE)
    :param fragment_path: path to PDB of current fragment (TARGET)
    :param lig_request_suffix: file name of target ligand PDB fetched from LigandExpo
    :param ligand_chain: Chain ID of target ligand used for alignment
    :param ligand_ResSeq_ID: Resnum for target ligand used for alignment
    :param target_fragment_atom_names: names of atoms in target fragment
    """
    def __init__(self, user_defined_dir, current_fragment, fragment_smiles_sanitized):  # todo: fragment_string really shouldn't be a string...

        # --- Paths --- #
        self.current_fragment = current_fragment
        self.frag_rigid_pdb_name = f'{current_fragment}-rigid.pdb'

        self.user_defined_dir = user_defined_dir
        self.fragment_inputs_path = os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs')

        # Select rigid atoms from fragment map for alignments if rigid atoms are defined in the Fragment_Inputs directory
        self.frag_inputs_dir = os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs', 'Rigid_Fragment_Atoms')

        # --- Fragment-related variables --- #

        self.fragment_string = open(os.path.join(self.fragment_inputs_path, '{}.pdb'.format(current_fragment))).read()
        self.fragment_smiles_sanitized = fragment_smiles_sanitized
        self.fragment_prody = prody.parsePDBStream(io.StringIO(self.fragment_string)).select('not hydrogen')
        self.ideal_fragment_mapping = None  # Assigned in _init_map_fragment_smiles_onto_ideal()

        # RDKit Mol representations of fragment
        self.fragment_ideal_rdkit_mol = Chem.MolFromSmiles(self.fragment_smiles_sanitized)
        self.fragment_pdb_rdkit_mol = Chem.MolFromPDBBlock(self.fragment_string, removeHs=False)

        # Assert that mol objects were loaded properly
        assert not any([self.fragment_pdb_rdkit_mol is None, self.fragment_ideal_rdkit_mol is None])

        # --- Run initialization functions --- #

        self._init_map_fragment_smiles_onto_ideal()

    def _init_map_fragment_smiles_onto_ideal(self):
        """
        Maps the fragment SMILES string onto the user-defined fragment PDB
        :return:
        """
        # --- Map Fragment SMILES onto fragment PDB --- #

        fragment_mcs = rdFMCS.FindMCS([self.fragment_ideal_rdkit_mol, self.fragment_pdb_rdkit_mol], bondCompare=rdFMCS.BondCompare.CompareAny)
        fragment_substructure_mol = Chem.MolFromSmarts(fragment_mcs.smartsString)

        fragment_ideal_match = self.fragment_ideal_rdkit_mol.GetSubstructMatch(fragment_substructure_mol)
        fragment_pdb_match = self.fragment_pdb_rdkit_mol.GetSubstructMatch(fragment_substructure_mol)

        # Assert 1:1 mapping excluding hydrogens
        assert len(fragment_ideal_match) == len(fragment_pdb_match)

        self.ideal_fragment_mapping = {i_idx: p_idx for i_idx, p_idx in zip(fragment_ideal_match, fragment_pdb_match)}

    def extract_ligand_records(self, pdb_file, relevant_ligands):
        """
        Extract all instances of fragment-containing small molecules bound to a PDB

        :param pdb_file: Path to the PDB file containing the fragment-containing ligand
        :param relevant_ligands: dict of Ideal_Ligand_PDB_Containers for ligands to find in pdb_file
        :return: relevant_ligands_prody_dict - structured [resname, chain, resnum] : ligand prody
        """
        pdb_header = prody.parsePDB(pdb_file, model=0, header=True)

        # Find resnums for all instances of fragment-containing small molecules bound to this PDBrelevant_ligands
        relevant_ligand_resnums = [(res.chain, res.resnum, res) for res in pdb_header['chemicals'] if res.resname in relevant_ligands.keys()]

        # Pull all relevant ligands from PDB as a prody AtomGroup objects
        pdb_prody_hv = prody.parsePDB(pdb_file, altloc=True).getHierView()
        relevant_ligands_prody_dict = dict()

        # todo: only consider biological assemblies, we want to avoid duplicate observations in the fuzzball
        for ligand_chain, ligand_resnum, res in relevant_ligand_resnums:

            ligand_pdb_prody = pdb_prody_hv.getResidue(ligand_chain, ligand_resnum)

            # Issues occur when alternate coordinates are parsed for a ligand... we don't want those anyways
            if ligand_pdb_prody is None:
                continue

            """
            Still not sure what metrics work the best for all PDBs.
            
            I have encountered some structures where atom occupancy is 0 but coordinates are extrapolated from fully 
            resolved atoms. For now I'm just going to assume these instances are fine...
            
            20180114 - I'm starting to think that these metrics aren't necessary... see my group meeting for details. 
            Actually I'll just copy/paste here after I finish the slide (let's hope I remember...).
            """

            relevant_ligands_prody_dict[(res.resname, ligand_chain, ligand_resnum)] = ligand_pdb_prody

            # Check for multiple/missing occupancies
            # ligand_pdb_prody_without_hydrogens = ligand_pdb_prody.select('not hydrogen')
            # full_occupancy = [atom for atom in ligand_pdb_prody_without_hydrogens.getOccupancies() if atom == 1]
            # missing_occupancy = [atom for atom in ligand_pdb_prody_without_hydrogens.getOccupancies() if atom == 0]

            # fail_conditions_list = [ligand_pdb_prody_without_hydrogens.numCoordsets() > 1,              # Multiple occupancies
            #                         len(full_occupancy) + len(missing_occupancy) < len(ligand_pdb_prody_without_hydrogens),
            #                         len(ligand_pdb_prody_without_hydrogens) != len(relevant_ligands[res.resname].ideal_ligand_prody.select('not hydrogen'))  # Missing atoms (compared to ideal structure)
            #                         ]
            #
            # if not any(fail_conditions_list):
            #     relevant_ligands_prody_dict[(res.resname, ligand_chain, ligand_resnum)] = ligand_pdb_prody

            # DEBUGGING
            #
            # else:
            #     conditions = ['Multiple Occupancies', 'Missing Atoms', 'Missing Atoms']
            #     pprint(ligand_pdb_prody_without_hydrogens.numCoordsets())
            #     pprint([(a, b) for a, b in zip(conditions, fail_conditions_list)])
            #     pprint([len(ligand_pdb_prody_without_hydrogens), len(relevant_ligands[res.resname].ideal_ligand_prody.select('not hydrogen'))])
            #     pprint([a for a in ligand_pdb_prody_without_hydrogens])
            #     pprint([a for a in relevant_ligands[res.resname].ideal_ligand_prody.select('not hydrogen')])
            #
            #     if ligand_pdb_prody_without_hydrogens.numCoordsets() > 1:
            #         for coordset in ligand_pdb_prody_without_hydrogens.iterCoordsets():
            #             print(coordset)

        return relevant_ligands_prody_dict

    def fragment_target_mapping(self, target_ligand_ideal_smiles, target_ligand_pdb_string):
        """
        I talked to Chris today and he showed me this chemoinformatics package called RDKit. This may be a much
        more streamlined method of doing this without needing to resort to subgraphs...
        > https://sourceforge.net/p/rdkit/mailman/message/34549645/

        I need to pass in...
        * target_ligand_ideal_smiles: SMILES for ligand puled from LigandExpo                   -> self.target_ligand_dict['smiles']
        * target_ligand_pdb_string: Ligand pulled from the current PDB. HYDROGENS REMOVED.      -> self.target_string

        self.fragment_smiles_sanitized: Sanitized SMILES string for the user-defined fragment
        self.fragment_string: String representation of user-defined fragment in PDB format

        target_ligand_pdb_string: String representation of the target ligand pulled from current PDB in PDB format. HYDROGENS REMOVED.

        :return: dict mapping fragment atoms to target atoms
        """

        """
        20180516
        
        # todo: update fragment_target_mapping to use ideal ligands for mapping
        
        I've found that I lose atom valency information when using ligands pulled straight from the target PDBs. For
        instance, all information about double bonds are lost... obviously this is important information for finding
        correct substructures. FindMCS has a matchValences option that I want to use, but MolFromPDBBlock does not
        recover double bonds and hydrogens, leading to incorrect valences despite correct heavy atom connectivity.
        I should use the ideal ligand PDB from LigandExpo to do the mapping, and then map the correct atoms from the
        ideal PDB onto the ligand pulled from the target PDB.

        - Ideal ligand from SMILES string
        - Fragment from SMILES string
        - PDB ligand from coordinates
        """

        # --- Map ideal ligand onto ligand pulled from current PDB --- #

        ligand_ideal = Chem.MolFromSmiles(target_ligand_ideal_smiles)
        ligand_pdb = Chem.MolFromPDBBlock(target_ligand_pdb_string.getvalue(), removeHs=False)

        if ligand_pdb is None:
            print('\nUnable to load ligand PDB as an RDKit mol object...\n')
            return False, None

        # I need to generalize bond types since importing mol from a PDB fucks up valence information and things don't
        # map properly: https://github.com/rdkit/rdkit/issues/943

        # todo: convert this into a static method so I can use it in Generate_Fuzzball_using_PyRosetta.assemble_defined_fuzzball
        ligand_mcs = rdFMCS.FindMCS([ligand_ideal, ligand_pdb], bondCompare=rdFMCS.BondCompare.CompareAny)
        ligand_substructure_mol = Chem.MolFromSmarts(ligand_mcs.smartsString)

        ligand_ideal_match = ligand_ideal.GetSubstructMatch(ligand_substructure_mol)
        ligand_pdb_match = ligand_pdb.GetSubstructMatch(ligand_substructure_mol)

        # Assert 1:1 mapping excluding hydrogens
        if len(ligand_ideal_match) != len([a for a in ligand_ideal.GetAtoms()]):
            print('\nFull ligand was not mapped!\n')
            return False, None

        if len(ligand_ideal_match) != len(ligand_pdb_match):
            return False, None

        ideal_pdb_mapping = {i_idx: p_idx for i_idx, p_idx in zip(ligand_ideal_match, ligand_pdb_match)}

        # --- Map fragment SMILES onto ideal ligand SMILES --- #

        # Generate a RDKit mol object of the common substructure so that I can map the fragment and target onto it
        substructure_match = rdFMCS.FindMCS([self.fragment_ideal_rdkit_mol, ligand_ideal])
        sub_mol = Chem.MolFromSmarts(substructure_match.smartsString)

        # Fall back onto mapping without considering bond order
        # This should catch instances where the full fragment can't be mapped due to inconsistencies with the provided
        # SMILES strings... I've found that this is particularly an issue for matching aromatic bonds since the SMILES
        # strings are converted to kekule structures.

        if len(sub_mol.GetAtoms()) != len(self.ideal_fragment_mapping):
            print('\nUnable to match fragment with correct bond orders, let\'s try ignoring them...\n')
            substructure_match = rdFMCS.FindMCS([self.fragment_ideal_rdkit_mol, ligand_ideal], bondCompare=rdFMCS.BondCompare.CompareAny)
            sub_mol = Chem.MolFromSmarts(substructure_match.smartsString)

            if len(sub_mol.GetAtoms()) != len(self.ideal_fragment_mapping):
                print('\nFailed to correctly map fragment onto target molecule...\n')
                return False, None

        # Return atom indicies of substructure matches for fragment and target
        frag_matches = self.fragment_ideal_rdkit_mol.GetSubstructMatch(sub_mol)
        target_matches = ligand_ideal.GetSubstructMatches(sub_mol)

        # DEBUGGING
        # from rdkit.Chem import Draw
        # from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
        # DrawingOptions.includeAtomNumbers = True
        #
        # Draw.MolToFile(self.fragment_ideal_rdkit_mol, 'fragment_ideal.png')
        # Draw.MolToFile(self.fragment_pdb_rdkit_mol, 'fragment_pdb.png')
        # Draw.MolToFile(sub_mol, 'fragment_substructure_mol.png')
        # Draw.MolToFile(ligand_ideal, 'ligand_ideal.png')
        # Draw.MolToFile(ligand_pdb, 'ligand_pdb.png')

        # --- Map fragment atom index to target atom index --- #
        fragment_target_map = [[(self.ideal_fragment_mapping[f_idx], ideal_pdb_mapping[t_idx]) for f_idx, t_idx in zip(frag_matches, target_match)] for target_match in target_matches]

        return True, fragment_target_map

    def process_atom_mappings_into_coordinate_sets(self, fragment_target_map, target_prody):
        """

        :param fragment_target_map:
        :param target_prody:
        :return:
        """
        # Retrieve fragment and target atom indicies to align
        fragment_atom_indices = [a[0] for a in fragment_target_map]
        target_atom_indices = [a[1] for a in fragment_target_map]

        # Convert atom indicies into atom objects
        frag_atom_selections = [self.fragment_prody.select('index {}'.format(index)) for index in fragment_atom_indices]

        ligand_atoms = [atom for atom in target_prody]
        trgt_atom_selections = [ligand_atoms[index] for index in target_atom_indices]
        # trgt_atom_selections = [target_prody.select('index {}'.format(index)) for index in target_atom_indices]

        # Save atom names for mapped fragment atoms in target ligand
        # target_fragment_atom_names = ' '.join([atom.getNames()[0] for atom in trgt_atom_selections if atom != None])  # This isn't used anywhere... probably debugging at one point
        target_fragment_atom_serials = [str(atom.getSerial()) for atom in trgt_atom_selections if atom != None]

        # Get atom coordinates out of atom objects ({atom != None} because hydrogens were removed)
        frag_atom_coords = np.asarray([atom.getCoords()[0] for atom in frag_atom_selections if atom != None])
        trgt_atom_coords = np.asarray([atom.getCoords() for atom in trgt_atom_selections if atom != None])

        return frag_atom_coords, trgt_atom_coords, target_fragment_atom_serials

    def determine_rotation_and_translation(self, fragment_target_map, target_prody):
        """
        Implementing the Kabsch algorithm for aligning all fragment-containing small molecules to the target ligand
        on the mapped atoms as determined by fragment_target_mapping()

        :return: Prody transformation object with rotation and translation matrix/vector
        """
        frag_atom_coords, trgt_atom_coords, target_fragment_atom_serials = self.process_atom_mappings_into_coordinate_sets(fragment_target_map, target_prody)

        if os.path.exists(self.frag_inputs_dir) and self.frag_rigid_pdb_name in os.listdir(self.frag_inputs_dir):
            print('* Using defined rigid atoms for alignment *')
            frag_atom_rigid, trgt_atom_rigid = self.return_rigid_atoms(frag_atom_coords, trgt_atom_coords)
            return (trgt_atom_rigid, target_fragment_atom_serials), frag_atom_rigid, prody.calcTransformation(trgt_atom_rigid, frag_atom_rigid)

        else:
            return (trgt_atom_coords, target_fragment_atom_serials), frag_atom_coords, prody.calcTransformation(trgt_atom_coords, frag_atom_coords)

    def return_rigid_atoms(self, frag_atom_coords, trgt_atom_coords):
        """
        Returns rigid atoms as defined in Fragment_Inputs
        
        This relies on assumption that fragment and rigid atoms are in same coordinate frame
        USE SAME INPUT PDB to generate fragments and rigid atom selections!!!
        
        Saved Atoms     V       V   V   V           V       V       V
        Rigid Atoms   | * |   | * | * | * |   |   | * |   | * |   | * |
        Full Fragment | * | * | * | * | * | * | * | * | * | * | * | * |
        Mapped Target | * | * | * | * | * | * | * | * | * | * | * | * |
        
        :param current_fragment: 
        :param frag_atom_coords:
        :param trgt_atom_coords:
        :return: 
        """
        frag_inputs_dir = os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs', 'Rigid_Fragment_Atoms')
        frag_rigid_pdb_name = '{}-rigid.pdb'.format(self.current_fragment)

        frag_atom_rigid = []
        trgt_atom_rigid = []

        # Import rigid atoms from fragment inputs directory
        fragment_rigid_prody = prody.parsePDB(os.path.join(frag_inputs_dir, frag_rigid_pdb_name))

        # Get coordinates from rigid atoms
        fragment_rigid_coords = [atom.getCoords() for atom in fragment_rigid_prody]

        # For each point in fragment
        for frag_coord, trgt_coord in zip(frag_atom_coords, trgt_atom_coords):

            # If point present in rigid atoms, add frag atom and corresponding mapped atom to new coord array
            if any([np.array_equal(frag_coord, rigid_coord) for rigid_coord in fragment_rigid_coords]):
                frag_atom_rigid.append(frag_coord)
                trgt_atom_rigid.append(trgt_coord)

        return np.asarray(frag_atom_rigid), np.asarray(trgt_atom_rigid)

    def apply_transformation(self, pdb_file, ligand_resnum, target_fragment_atom_serials, transformation_matrix):
        """
        Apply transformation to the target ligand-protein complex.

        Also considering:
        * Only work with residues with CA within 8A of ligand
        * Write all transformed PDBs to a new working directory?

        :param transformation_matrix: 
        :param target_pdb: 
        :return: 
        """
        # Only work with residues within 12A of target ligand
        target_pdb = prody.parsePDB(pdb_file)
        target_shell = target_pdb.select('(protein and within 12 of (serial {0}) and not resnum {1}) or (serial {0})'.format(' '.join(target_fragment_atom_serials), ligand_resnum))

        transformed_pdb = prody.applyTransformation(transformation_matrix, target_shell)

        return transformed_pdb

class Ideal_Ligand_PDB_Container(object):
    """
    This class is responsible for aligning all fragment-containing small molecule structures to the corresponding
    fragments.

    :param resname: three-letter chemical component identifier for the target ligand
    :param processed_PDBs_path: directory where all processed aligned PDBs will be deposited
    """

    def __init__(self, resname):
        self.resname = resname.upper()
        self.success = True

        # Initialized later
        self.target_ligand_dict = None
        self.ideal_ligand_prody = None
        self.smiles = None

        self.fetch_ideal_ligand_info()

    def __repr__(self):
        return f'{self.resname}'

    def fetch_ideal_ligand_info(self):
        """
        Get Ligand information from the PDB (SMILES strings specifically)
        Note: you can also parse .cif representation of ligand to find SMILES {http://files.rcsb.org/ligands/view/{ligand_ID}.cif}
        """
        pdb_ligand_dict = prody.fetchPDBLigand(self.resname)

        # Set ideal ligand prody
        if pdb_ligand_dict.get('ideal', f'Ideal ligand representation for {self.resname} was not found!'):
            self.ideal_ligand_prody = pdb_ligand_dict.get('ideal')
        else:
            self.success = False

        # Set SMILES
        if pdb_ligand_dict.get('OpenEye_OEToolkits_SMILES_CANONICAL', 'OpenEye SMILES representation was not found!'):
            self.smiles = pdb_ligand_dict.get('OpenEye_OEToolkits_SMILES_CANONICAL')
        elif pdb_ligand_dict.get('CACTVS_SMILES_CANONICAL'):
            self.smiles = pdb_ligand_dict.get('CACTVS_SMILES_CANONICAL', 'CACTVS SMILES representation was not found!')
        else:
            self.success = False

        # Save for now because why not
        self.target_ligand_dict = pdb_ligand_dict
