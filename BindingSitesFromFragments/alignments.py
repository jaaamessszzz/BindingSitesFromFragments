#!/usr/bin/env python3

import os
import prody
import io
import requests
import sys
import time
import numpy as np
import pprint

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdFMCS
from rdkit.Chem.Fingerprints import FingerprintMols

class Align_PDB():
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
    def __init__(self, user_defined_dir, pdb_file=None, target_string=None, fragment_string=None):
        self.user_defined_dir = user_defined_dir
        self.pdb_file = pdb_file # Mobile PDB File to be aligned to a target
        self.target_string = target_string
        self.fragment_string = fragment_string
        self.target_prody = None
        self.fragment_prody = None
        self.lig_request_suffix = None
        self.ligand_chain = None
        self.ligand_ResSeq_ID = None
        self.target_fragment_atom_names = None
        self.target_fragment_atom_serials = None

    def extract_ligand_records(self, pdb_file):
        """
        Ligand records for all instances of a fragment-containing small molecule bound to a protein

        :param ligand: Three letter code for fragment-containing ligand
        :param pdb_file: Path to the PDB file containing the fragment-containing ligand
        :return: String of HETAM and CONECT records from the input pdb for the given ligand
        """
        # Identify all instances of target ligand in PDB
        pdb_header = prody.parsePDBHeader(pdb_file)
        pdbid = pdb_header['identifier']

        # Report
        print(self.ligand)
        print(pdb_file)

        # Start dict with ChID_resnum as keys
        # Get ligand ChID and resnum from prody header dict
        target_dict = {}
        for ligand in pdb_header['chemicals']:
            if ligand.resname == self.ligand:
                target_dict['{}_{}'.format(ligand.chain, ligand.resnum)] = {}
                target_dict['{}_{}'.format(ligand.chain, ligand.resnum)]['atom_records'] = []
                target_dict['{}_{}'.format(ligand.chain, ligand.resnum)]['atom_indices'] = []

        # Create new list for Conect records
        conect_record_list = []

        # Go through PDB file ONCE
        with open(pdb_file, 'r') as current_pdb_file:
            for line in current_pdb_file:

                # If ChID and resnum encountered for a ligand of interest
                if line[:6] == 'HETATM':
                    if '{}_{}'.format(line[21:22].strip(), line[22:26].strip()) in target_dict.keys():
                        atom_chain = line[21:22].strip()
                        atom_resnum = line[22:26].strip()

                        # Record atom record in list in my ChID_resnum dict
                        target_dict['{}_{}'.format(atom_chain, atom_resnum)]['atom_records'].append(line)
                        # Record atom indicies in list in my ChID_resnum dict
                        target_dict['{}_{}'.format(atom_chain, atom_resnum)]['atom_indices'].append(line[7:11].strip())

                # Add list of conect records for each line
                if line[:6] == "CONECT":
                    conect_record_list.append(line)

        # Evaluate each ChID_resnum
        # todo: update to extract all matching ligands from a PDB, currently only extracts the first ligand that passes
        for potential_ligand in target_dict:

            passing = self._verify_pdb_quality(target_dict[potential_ligand]['atom_records'], pdb_file)
            print("Passing: {}".format(passing))

            if passing:
                # If ligand passes evaluation...
                # Go through conect records and add if all atoms in atom index list
                for conect_record in conect_record_list:
                    if any([index in conect_record.split() for index in target_dict[potential_ligand]['atom_indices']]):
                        target_dict[potential_ligand]['atom_records'].append(conect_record)

                # If everything is all good, assign variables
                self.pdb_file = pdb_file
                self.target_string = ''.join(target_dict[potential_ligand]['atom_records'])
                self.ligand_chain = potential_ligand.split('_')[0]
                self.ligand_ResSeq_ID = potential_ligand.split('_')[1]

                # http://ligand-expo.rcsb.org/ld-download.html
                self.lig_request_suffix = '{}_{}_{}_{}'.format(pdbid,
                                                               self.ligand,
                                                               self.ligand_chain,
                                                               self.ligand_ResSeq_ID
                                                               )  # PDBID, Ligand, ChID, Resnum

                self.ligand_atom_indicies = target_dict[potential_ligand]['atom_indices']

                return False

        # If I can't find any ligands fully represented, reject this PDB
        print('\n{} REJECTED! SAD!\n'.format(pdbid))
        return True
        
    def _verify_pdb_quality(self, ligand_record_list, pdb_file):
        """
        Verifies PDB quailty
        * ligand does not have atoms with multiple occupancies
        * ligand and source PDB can be parsed with prody
         
        :param potential_ligand_block: string of HETATM PDB records
        :return: True if the ligand passes QC, else False
        """
        multiple_occupancies = self._check_multiple_occupancies(ligand_record_list)

        print('Multiple Occupancies: {}'.format(multiple_occupancies))
        if not multiple_occupancies:

            # Compare ideal to ligand pulled from Ligand Expo
            # prody is unable to parse empty pdb files, although who would even think to test that
            atom_count = 0

            for line in ligand_record_list:
                if line[:6] == 'HETATM' and line[77].strip() != 'H':
                    atom_count += 1

            # If the number of atoms are the same, good enough for me (at the moment)
            ideal_ligand_prody = prody.parsePDBStream(self.ideal_ligand_pdb)
            if atom_count == ideal_ligand_prody.select('not hydrogen').numAtoms():

                # todo: find a more elegant way of testing if prody can open these
                # Verify whether I can actually open the ligand pdb with Prody...
                # Test the source protein-ligand complex PDB at the same time...
                try:
                    prody.parsePDBStream(io.StringIO(''.join(ligand_record_list)))
                    prody.parsePDB(pdb_file)
                    return True

                except Exception as e:
                    print(e)
                    return False

        return False

    def _check_multiple_occupancies(self, ligand_record_list):
        """
        Checks to make sure there aren't multiple occupancies for any atoms in the ligand

        :param target_string: string of ligand pdb contents
        :return: 
        """
        multiple_occupancies = False
        for line in ligand_record_list:
            if line[:6] == 'HETATM' and line[16:17] != ' ':
                multiple_occupancies = True
                break

        return multiple_occupancies

    def fragment_target_mapping(self, fragment_smiles):
        """
        I talked to Chris today and he showed me this chemoinformatics package called RDKit. This may be a much
        more streamlined method of doing this without needing to resort to subgraphs...
        > https://sourceforge.net/p/rdkit/mailman/message/34549645/

        :return: dict mapping fragment atoms to target atoms
        """
        # 20180516
        # todo: update fragment_target_mapping to use ideal ligands for mapping
        # I've found that I lose atom valency information when using ligands pulled straight from the target PDBs. For
        # instance, all information about double bonds are lost... obviously this is important information for finding
        # correct substructures. FindMCS has a matchValences option that I want to use, but MolFromPDBBlock does not
        # recover double bonds and hydrogens, leading to incorrect valences despite correct heavy atom connectivity.
        # I should use the ideal ligand PDB from LigandExpo to do the mapping, and then map the correct atoms from the
        # ideal PDB onto the ligand pulled from the target PDB.

        # Ideal ligand from SMILES string
        # Fragment from SMILES string
        # PDB ligand from coordinates

        # Map ideal ligand to pdb ligand, assert 1:1 mapping excluding hydrogens
        ligand_ideal = Chem.MolFromSmiles(self.target_ligand_dict['smiles'])
        ligand_pdb = Chem.MolFromPDBBlock(self.target_string, removeHs=True)

        # I need to generalize bond types since importing mol from a PDB fucks up valence information and things don't
        # map properly: https://github.com/rdkit/rdkit/issues/943

        ligand_mcs = rdFMCS.FindMCS([ligand_ideal, ligand_pdb], bondCompare=rdFMCS.BondCompare.CompareAny)
        ligand_substructure_mol = Chem.MolFromSmarts(ligand_mcs.smartsString)

        ligand_ideal_match = ligand_ideal.GetSubstructMatch(ligand_substructure_mol)
        ligand_pdb_match = ligand_pdb.GetSubstructMatch(ligand_substructure_mol)

        if len(ligand_ideal_match) != len(ligand_pdb_match):
            return False

        ideal_pdb_mapping = {i_idx: p_idx for i_idx, p_idx in zip(ligand_ideal_match, ligand_pdb_match)}

        # Map Fragment SMILES onto fragment PDB
        fragment_ideal = Chem.MolFromSmiles(fragment_smiles)
        fragment_pdb = Chem.MolFromPDBBlock(self.fragment_string, removeHs=True)

        fragment_mcs = rdFMCS.FindMCS([fragment_ideal, fragment_pdb], bondCompare=rdFMCS.BondCompare.CompareAny)
        fragment_substructure_mol = Chem.MolFromSmarts(fragment_mcs.smartsString)

        fragment_ideal_match = fragment_ideal.GetSubstructMatch(fragment_substructure_mol)
        fragment_pdb_match = fragment_pdb.GetSubstructMatch(fragment_substructure_mol)

        if len(fragment_ideal_match) != len(fragment_pdb_match):
            return False

        ideal_fragment_mapping = {i_idx: p_idx for i_idx, p_idx in zip(fragment_ideal_match, fragment_pdb_match)}

        # Map fragment SMILES onto ideal ligand SMILES
        # Generate a RDKit mol object of the common substructure so that I can map the fragment and target onto it
        substructure_match = rdFMCS.FindMCS([fragment_ideal, ligand_ideal])
        sub_mol = Chem.MolFromSmarts(substructure_match.smartsString)

        # Return atom indicies of substructure matches for fragment and target
        frag_matches = fragment_ideal.GetSubstructMatch(sub_mol)
        target_matches = ligand_ideal.GetSubstructMatches(sub_mol)

        # Maps fragment atom index to target atom index
        # If there is more than one substructure match for the target, find the one with the lowest RMSD to the fragment
        # todo: use all mappings since each will have a unique local environment in the binding site, will also solve symmetric fragment issue

        # This will map fragment atoms to pdb atom indicies
        # fragment_target_map = [[(f_idx, t_idx) for f_idx, t_idx in zip(frag_matches, target_match)] for target_match in target_matches]
        fragment_target_map = [[(ideal_fragment_mapping[f_idx], ideal_pdb_mapping[t_idx]) for f_idx, t_idx in zip(frag_matches, target_match)] for target_match in target_matches]

        # if len(target_matches) > 1:
        #     fragment_target_map = self.identify_best_substructure(frag_matches, target_matches)
        #
        # else:
        #     fragment_target_map = [(f_idx, t_idx) for f_idx, t_idx in zip(frag_matches, target_matches[0])]

        # Assign successful mappings to self
        self.fragment_target_map = fragment_target_map

        # Fragment and target as prody atom objects
        self.target_prody = prody.parsePDBStream(io.StringIO(self.target_string)).select('not hydrogen')
        self.fragment_prody = prody.parsePDBStream(io.StringIO(self.fragment_string)).select('not hydrogen')

        return True

    def identify_best_substructure(self, frag_matches, target_matches):
        """
        Identifies the "correct" substructure from a RDKit GetSubstructMatches search.

        For each substructure match, retrieve the target atoms and calculate the RMSD to the fragment atoms. The target
        substructure with the lowest RMSD to the fragment is the one that will be used downstream for alignment.

        This method will not help discriminate decoys for molecules with several repeating substructures... 
        like sugars... like in ONPF...

        :param frag_matches: tuple with atom indicies for the fragment substructure match
        :param target_matches: iterable of tuples containing atom indicies for each substructure hit
        :param fragment_pdb: path to fragment pdb
        :param target_string: string containing contents of target pdb
        :return: fragment mapping for the "correct" substructure
        """
        rmsd_dict = {}
        for match in target_matches:
            fragment_target_map = [(f_idx, t_idx) for f_idx, t_idx in zip(frag_matches, match)]
            frag_atom_coords, trgt_atom_coords = self.process_atom_mappings_into_coordinate_sets(fragment_target_map)
            transformation_matrix = prody.calcTransformation(trgt_atom_coords, frag_atom_coords)
            aligned_trgt = prody.applyTransformation(transformation_matrix, trgt_atom_coords)

            rmsd_dict[prody.calcRMSD(frag_atom_coords, aligned_trgt)] = fragment_target_map

        return rmsd_dict[min(rmsd_dict.keys())]

    def process_atom_mappings_into_coordinate_sets(self, fragment_target_map):
        # Retrieve fragment and target atom indicies to align
        fragment_atom_indices = [a[0] for a in fragment_target_map[0]]
        target_atom_indices = [a[1] for a in fragment_target_map[0]]

        # Convert atom indicies into atom objects
        frag_atom_selections = [self.fragment_prody.select('index {}'.format(index)) for index in fragment_atom_indices]
        trgt_atom_selections = [self.target_prody.select('index {}'.format(index)) for index in target_atom_indices]

        # DEBUGGING - Export PDBs
        # os.makedirs(os.path.join('ONPF', 'Mapped_atoms'), exist_ok=True)
        # trgt_mapped_pdb = self.target_prody.select('index {}'.format(' '.join([str(a) for a in target_atom_indices])))
        # frag_mapped_pdb = self.fragment_prody.select('index {}'.format(' '.join([str(a) for a in fragment_atom_indices])))
        #
        # prody.writePDB(os.path.join('ONPF', 'Mapped_atoms', '{}-fragment.pdb'.format(self.pdb_file)), frag_mapped_pdb)
        # prody.writePDB(os.path.join('ONPF', 'Mapped_atoms', '{}-target.pdb'.format(self.pdb_file)), trgt_mapped_pdb)

        # Save atom names for mapped fragment atoms in target ligand
        self.target_fragment_atom_names = ' '.join([atom.getNames()[0] for atom in trgt_atom_selections if atom != None])
        self.target_fragment_atom_serials = ' '.join([str(atom.getSerials()[0]) for atom in trgt_atom_selections if atom != None])

        # Get atom coordinates out of atom objects
        # None is because hydrogens were removed
        frag_atom_coords = np.asarray([atom.getCoords()[0] for atom in frag_atom_selections if atom != None])
        trgt_atom_coords = np.asarray([atom.getCoords()[0] for atom in trgt_atom_selections if atom != None])

        return frag_atom_coords, trgt_atom_coords

    def determine_rotation_and_translation(self, current_fragment=None):
        """
        Implementing the Kabsch algorithm for aligning all fragment-containing small molecules to the target ligand
        on the mapped atoms as determined by fragment_target_mapping()

        :return: Prody transformation object with rotation and translation matrix/vector
        """
        # todo: implement an option for RMSD cutoff where mapped atoms do not necessarily have the same conformations, e.g. sugar pucker
        frag_atom_coords, trgt_atom_coords = self.process_atom_mappings_into_coordinate_sets(self.fragment_target_map)

        # Select rigid atoms from fragment map for alignments if rigid atoms are defined in the Fragment_Inputs directory
        frag_inputs_dir = os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs', 'Rigid_Fragment_Atoms')
        frag_rigid_pdb_name = '{}-rigid.pdb'.format(current_fragment)

        if os.path.exists(frag_inputs_dir) and frag_rigid_pdb_name in os.listdir(frag_inputs_dir):
            print('* Using defined rigid atoms for alignment *')
            frag_atom_rigid, trgt_atom_rigid = self.return_rigid_atoms(current_fragment, frag_atom_coords, trgt_atom_coords)
            return trgt_atom_rigid, frag_atom_rigid, prody.calcTransformation(trgt_atom_rigid, frag_atom_rigid)

        else:
            return trgt_atom_coords, frag_atom_coords, prody.calcTransformation(trgt_atom_coords, frag_atom_coords)

    def return_rigid_atoms(self, current_fragment, frag_atom_coords, trgt_atom_coords):
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
        frag_rigid_pdb_name = '{}-rigid.pdb'.format(current_fragment)

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


    def apply_transformation(self, transformation_matrix):
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
        target_pdb = prody.parsePDB(self.pdb_file)
        ligand_resnum = target_pdb.select('serial {0}'.format(self.target_fragment_atom_serials)).getResnums()[0]
        target_shell = target_pdb.select('(protein and within 12 of (serial {0}) and not resnum {1}) or (serial {0})'.format(self.target_fragment_atom_serials, ligand_resnum))

        transformed_pdb = prody.applyTransformation(transformation_matrix, target_shell)
        prody.writePDB(os.path.join(self.processed_PDBs_path, '{}-processed.pdb'.format(self.lig_request_suffix)), transformed_pdb)

        return transformed_pdb


class Fragment_Alignments(Align_PDB):
    """
    This subclass class is responsible for aligning all fragment-containing small molecule structures to the corresponding 
    fragments.
    """

    def __init__(self, user_defined_dir, ligand='DERP', processed_PDBs_path=None, pdb_file=None, target_string=None, fragment_string=None):
        Align_PDB.__init__(self, user_defined_dir, pdb_file=pdb_file, target_string=target_string, fragment_string=fragment_string)
        self.ligand = ligand.upper()
        self.processed_PDBs_path = processed_PDBs_path
        # Initialized later
        self.ideal_ligand_pdb = None
        self.target_ligand_dict = None

    def fetch_records(self):
        """
        Only download ligands from Ligand Expo if and when required...
        """
        self.ideal_ligand_pdb = self.fetch_ideal_pdb(self.ligand)

    def fetch_ideal_pdb(self, ligand):
        """
        Get the ideal pdb for the target ligand. This is used to check that all atoms are present in target ligands
        extracted from bound proteins

        :param ligand: Three letter code for the desired ligand (upper case!!)
        :return: io.StringIO for pdb of ligand retrieved from LigandExpo
        """
        ideal_pdb_text = requests.get(
            'http://ligand-expo.rcsb.org/reports/{}/{}/{}_ideal.pdb'.format(ligand[0], ligand, ligand),
            stream=False).text
        if ideal_pdb_text != "":
            return io.StringIO(ideal_pdb_text)

        # For whenever I figure out how to install ProDy 1.9 without getting stupid build errors...
        # else:
        #     ideal_cif_text = requests.get(
        #         'https://files.rcsb.org/ligands/view/{}.cif'.format(ligand), stream=False).text
        #     return io.StringIO(ideal_cif_text)


class Second_Shell_Alignments(Align_PDB):
    """
    Class for handling prerequisite information on motif residues for determining second shell interactions 
    """

    def __init__(self, motif_prody_path):
        Align_PDB.__init__(self, pdb_file=None, target_string=None, fragment_string=None)
        self.motif_prody = prody.parsePDB(motif_prody_path)
