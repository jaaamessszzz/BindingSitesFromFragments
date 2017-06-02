#!/usr/bin/env python3

import os
import prody
from rdkit import Chem
from rdkit.Chem import rdFMCS
import io
from bs4 import BeautifulSoup
import requests
import sys
import time
import numpy as np
import pprint


class Alignments():
    """
    This base class is responsible for aligning all fragment-containing small molecule structures to the corresponding 
    fragments.
    """

    def __init__(self, ligand='DERP', processed_PDBs_path=None):
        self.ligand = ligand.upper()
        self.processed_PDBs_path = processed_PDBs_path
        # Initialized later
        self.ligexpo_list = None
        self.ideal_ligand_pdb = None

    def fetch_records(self):
        """
        Only download ligands from Ligand Expo if and when required...
        """
        self.ligexpo_list = self.fetch_ligand_records(self.ligand)
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
        return prody.parsePDBStream(io.StringIO(ideal_pdb_text))

    def fetch_ligand_records(self, ligand):
        """
        Retrieves ligand records from LigandExpo

        20170515: So I was writing up my quals proposal and I realized that I could download specific ligands for any 
        structure in the PDB using LigandExpo as outlined here:
        http://ligand-expo.rcsb.org/ld-download.html >> http://ligand-expo.rcsb.org/files/<HASH>/<CC_ID>/<FORMAT>/

        This might be the way to go since it removes the need for all this parsing...

        20170516: So I've found that free amino acids such as HIS are not retrievable using LigandExpo in this manner...
        I'm going to need to find a solution for this. 

        :param ligand: Three letter code for fragment-containing ligand
        :param pdb_file: Path to the PDB file containing the fragment-containing ligand
        :return: List of all ligand HETATM records
        """
        lig_request = requests.get('http://ligand-expo.rcsb.org/files/{}/{}/ipdb/'.format(ligand[0], ligand),
                                   stream=False)
        soupy_soup = BeautifulSoup(lig_request.text, 'html.parser')

        return [link.get('href') for link in soupy_soup.find_all('a') if 'ipdb' in link.get('href')]


class Align_PDB(Alignments):
    """
    Subclass with functions for aligning individual PDBs
    :param pdb_file: path to input PDB containing protein bound to fragment-containing target molecule (MOBILE)
    :param target_string: string of PDB for target ligand, fetched from LigandExpo (MOBILE)
    :param fragment_path: path to PDB of current fragment (TARGET)
    :param lig_request_suffix: file name of target ligand PDB fetched from LigandExpo
    :param ligand_chain: Chain ID of target ligand used for alignment
    :param ligand_ResSeq_ID: Resnum for target ligand used for alignment
    :param target_fragment_atom_names: names of atoms in target fragment
    """
    def __init__(self, Alignments, pdb_file=None, target_string=None, fragment_string=None):
        self.__dict__ = Alignments.__dict__
        self.pdb_file = pdb_file # Mobile PDB File to be aligned to a target
        self.target_string = target_string
        self.fragment_string = fragment_string
        self.lig_request_suffix = None
        self.ligand_chain = None
        self.ligand_ResSeq_ID = None
        self.target_fragment_atom_names = None

    def fetch_specific_ligand_record(self, pdb_file):
        """
        Fetch a specific ligand record from ligand expo

        :param pdb_file: path to the specific pdb file containing the desired ligand
        :return: 
        """
        # I have determined this is definitely not okay...
        # Case in point: http://ligand-expo.rcsb.org/files/4/4MZ/ipdb/1keq_4MZ_1_A_350__F__D_.ipdb
        # It is empty. What. Why.

        reject = False

        # Go through all ligands bound to PDBs
        pdbid = os.path.basename(os.path.normpath(pdb_file)).split('.')[0]
        for lig_request_suffix in self.ligexpo_list:

            # For my specific ligand of interest bound to the correct PDB...
            if pdbid.lower() in lig_request_suffix:

                # Report
                print("\n Checking: {}".format(lig_request_suffix))

                # Get ligand record from LigExpo
                full_lig_request = requests.get(
                    'http://ligand-expo.rcsb.org/files/{}/{}/ipdb/{}'.format(self.ligand[0], self.ligand, lig_request_suffix),
                    stream=False).text

                # Determine whether there are multiple occupancies for the target ligand
                # I am going to avoid using ligands with multiple occupancies as this suggests binding conformations
                # with low affinity and specificity
                multiple_occupancies = self.clean_ligand_HETATM_records(full_lig_request)

                print('Multiple Occupancies: {}'.format(multiple_occupancies))
                if not multiple_occupancies:

                    # ligexpo_ligand = prody.parsePDBStream(io.StringIO(full_lig_request))
                    # Literally going line by line to count atoms...
                    # prody is unable to parse empty pdb files, although who would even think to test that
                    atom_count = 0
                    for line in full_lig_request.split('\n'):
                        if line[:6] == 'HETATM' and line[77].strip() != 'H':
                            atom_count += 1

                    # Compare ideal to ligand pulled from Ligand Expo
                    # If the number of atoms are the same, good enough for me (at the moment)
                    if atom_count == self.ideal_ligand_pdb.select('not hydrogen').numAtoms():
                        # todo: find a more elegant way of testing if prody can open these
                        # Verify whether I can actually open the ligand pdb with Prody...
                        try:
                            prody.parsePDBStream(io.StringIO(full_lig_request))
                            # Test the source protein-ligand complex PDB at the same time...
                            prody.parsePDB(pdb_file)
                            # If everything is all good, assign variables
                            self.pdb_file = pdb_file
                            self.target_string = full_lig_request
                            self.lig_request_suffix = lig_request_suffix.split('.')[0].upper()
                            self.ligand_chain = lig_request_suffix.split('_')[3]
                            self.ligand_ResSeq_ID = lig_request_suffix.split('_')[4]

                            return reject

                        except Exception as e:
                            print(e)
                            continue

        # If I can't find any ligands fully represented, reject this PDB
        print('\n{} REJECTED! SAD!\n'.format(pdbid))
        return True

    def clean_ligand_HETATM_records(self, full_lig_request):
        """
        Removes alternate location indicators...
        UPDATE: now this function just checks to make sure there aren't multiple occupancies for any atoms in the ligand

        :param target_string: string of ligand pdb contents
        :return: 
        """
        multiple_occupancies = False
        target_string_split = full_lig_request.split('\n')
        for line in target_string_split:
            if line[:6] == 'HETATM' and line[16:17] != ' ':
                multiple_occupancies = True
                break

        return multiple_occupancies

    def fragment_target_mapping(self):
        """
        I talked to Chris today and he showed me this chemoinformatics package called RDKit. This may be a much
        more streamlined method of doing this without needing to resort to subgraphs...
        > https://sourceforge.net/p/rdkit/mailman/message/34549645/

        :return: dict mapping fragment atoms to target atoms
        """
        # Import fragment and target ligand
        # todo: standardize this so that both target and fragment inputs are the same type...
        fragment_mol = Chem.MolFromPDBBlock(self.fragment_string, removeHs=False)
        target_mol_H = Chem.MolFromPDBBlock(self.target_string, removeHs=False)

        # Some ligands retrieved from LigandExpo have weird valence issues with RDKit... need to look into that
        # 3dyb_AD0_1_A_500__B___.ipdb
        # Oh hey this ligand was an issue anyways. Yay. I guess.
        # Anyways, here's another if statement
        if target_mol_H == None:
            return False

        target_mol = Chem.AddHs(target_mol_H)
        target_mol_PDB_Block = Chem.MolToPDBBlock(target_mol)

        # Debugging
        # For some reason 3dyb_AD0_1_A_500__B___.ipdb cannot be imported with RDKit...
        # Skipping for now...
        if target_mol == None or fragment_mol == None:
            return False

        # Generate a RDKit mol object of the common substructure so that I can map the fragment and target onto it
        substructure_match = rdFMCS.FindMCS([fragment_mol, target_mol])
        sub_mol = Chem.MolFromSmarts(substructure_match.smartsString)

        # Return atom indicies of substructure matches for fragment and target
        frag_matches = fragment_mol.GetSubstructMatch(sub_mol)
        target_matches = target_mol.GetSubstructMatches(sub_mol)

        # Maps fragment atom index to target atom index
        # If there is more than one substructure match for the target, find the one with the lowest RMSD to the fragment
        if len(target_matches) > 1:
            fragment_target_mapping = self.identify_best_substructure(frag_matches, target_matches)
            # fragment_target_mapping = [(f_idx, t_idx) for f_idx, t_idx in zip(frag_matches, target_matches[0])]

        else:
            fragment_target_mapping = [(f_idx, t_idx) for f_idx, t_idx in zip(frag_matches, target_matches[0])]

        # Assign successful mappings to self
        self.fragment_target_map = fragment_target_mapping
        self.target_mol_PDB_Block = target_mol_PDB_Block

        # Debugging
        print(self.fragment_target_map)

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
            fragment_target_mapping = [(f_idx, t_idx) for f_idx, t_idx in zip(frag_matches, match)]
            frag_atom_coords, trgt_atom_coords = self.process_atom_mappings_into_coordinate_sets(fragment_target_mapping)
            transformation_matrix = prody.calcTransformation(trgt_atom_coords, frag_atom_coords)
            aligned_trgt = prody.applyTransformation(transformation_matrix, trgt_atom_coords)

            rmsd_dict[prody.calcRMSD(frag_atom_coords, aligned_trgt)] = fragment_target_mapping

        return rmsd_dict[min(rmsd_dict.keys())]

    def process_atom_mappings_into_coordinate_sets(self, fragment_target_mapping):
        # Fragment and target as prody atom objects (including Hydrogens)
        target_prody_H = prody.parsePDBStream(io.StringIO(self.target_string))
        fragment_prody_H = prody.parsePDBStream(io.StringIO(self.fragment_string))

        # Debugging
        print('trgt')
        pprint.pprint([atom for atom in target_prody_H])
        print('frag')
        pprint.pprint([atom for atom in fragment_prody_H])

        # Fragment and target prody selections (excluding hydrogens)
        target_prody = target_prody_H.select('not hydrogen')
        fragment_prody = fragment_prody_H.select('not hydrogen')

        # Retrieve fragment and target atom indicies to align
        fragment_atom_indices = [a[0] for a in fragment_target_mapping]
        target_atom_indices = [a[1] for a in fragment_target_mapping]

        # Convert atom indicies into atom objects
        frag_atom_selections = [fragment_prody.select('index {}'.format(index)) for index in fragment_atom_indices]
        trgt_atom_selections = [target_prody.select('index {}'.format(index)) for index in target_atom_indices]

        # Save atom names for mapped fragment atoms in target ligand
        self.target_fragment_atom_names = ' '.join([atom.getNames()[0] for atom in trgt_atom_selections if atom != None])

        # Get atom coordinates out of atom objects
        # None is because hydrogens were removed
        frag_atom_coords = np.asarray([atom.getCoords()[0] for atom in frag_atom_selections if atom != None])
        trgt_atom_coords = np.asarray([atom.getCoords()[0] for atom in trgt_atom_selections if atom != None])

        return frag_atom_coords, trgt_atom_coords

    def determine_rotation_and_translation(self):
        """
        Implementing the Kabsch algorithm for aligning all fragment-containing small molecules to the target ligand
        on the mapped atoms as determined by fragment_target_mapping()

        :return: Prody transformation object with rotation and translation matrix/vector
        """
        # todo: implement an option for RMSD cutoff where mapped atoms do not necessarily have the same conformations, e.g. sugar pucker
        frag_atom_coords, trgt_atom_coords = self.process_atom_mappings_into_coordinate_sets(self.fragment_target_map)

        # Debugging
        print(frag_atom_coords)
        print(trgt_atom_coords)

        return prody.calcTransformation(trgt_atom_coords, frag_atom_coords)

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
        # Only work with residues within 8A of target ligand
        target_pdb = prody.parsePDB(self.pdb_file)

        # todo: somehow convert target fragment indicies to atom names that I can use for selection
        # So the reason for these super long selection strings is that the serial numbers from LigandExpo don't match
        # up with the serial numbers fromthe PDBs for the same chain/residue/atoms...
        print(self.target_fragment_atom_names)

        # This step is so slow...
        target_shell = target_pdb.select('protein and within 12 of\
        (resname {0} and chain {1} and resnum {2} and name {3}) or\
        (resname {0} and chain {1} and resnum {2} and name {3})'
                                         .format(self.ligand,
                                                 self.ligand_chain,
                                                 self.ligand_ResSeq_ID,
                                                 self.target_fragment_atom_names
                                                 )
                                         )

        transformed_pdb = prody.applyTransformation(transformation_matrix, target_shell)
        prody.writePDB(os.path.join(self.processed_PDBs_path, '{}processed.pdb'.format(self.lig_request_suffix)), transformed_pdb)

        return transformed_pdb
