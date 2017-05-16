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


class Alignments():
    """
    This class is responsible for aligning all fragment-containing small molecule structures to the corresponding 
    fragments.
    
    Chris suggests looking into a package called RDKit for manipulating SMILES strings and structures
    
    """
    def __init__(self, user_defined_dir, fragment, ligand):
        self.user_defined_dir = user_defined_dir
        self.fragment = fragment
        self.ligand = ligand.upper()
        self.ligexpo_list = self.fetch_ligand_records(self.ligand)
        self.ideal_ligand_pdb = self.fetch_ideal_pdb(self.ligand)

        # http://ligand-expo.rcsb.org/reports/4/4MZ/4MZ_ideal.pdb


    def fetch_ligand_records(self, ligand):
        """
        Retrieves ligand records from LigandExpo
        
        20170515: So I was writing up my quals proposal and I realized that I could download specific ligands for any 
        structure in the PDB using LigandExpo as outlined here:
        http://ligand-expo.rcsb.org/ld-download.html >> http://ligand-expo.rcsb.org/files/<HASH>/<CC_ID>/<FORMAT>/
        
        This might be the way to go since it removes the need for all this parsing...
        
        :param ligand: Three letter code for fragment-containing ligand
        :param pdb_file: Path to the PDB file containing the fragment-containing ligand
        :return: List of all ligand HETATM records
        """
        time.sleep(0.1)
        lig_request = requests.get('http://ligand-expo.rcsb.org/files/{}/{}/ipdb/'.format(ligand[0], ligand), stream=False)
        soupy_soup = BeautifulSoup(lig_request.text, 'html.parser')


        return [link.get('href') for link in soupy_soup.find_all('a') if 'ipdb' in link.get('href')]


    def fetch_ideal_pdb(self, ligand):
        """
        Get the ideal pdb for the target ligand. This is used to check that all atoms are present in target ligands
        extracted from bound proteins
        
        :param ligand: 
        :return: 
        """
        time.sleep(0.1)
        ideal_pdb_text = requests.get('http://ligand-expo.rcsb.org/reports/{}/{}/{}_ideal.pdb'.format(ligand[0], ligand, ligand), stream=False).text
        return prody.parsePDBStream(io.StringIO(ideal_pdb_text))


    def fetch_specific_ligand_record(self, pdb_file):
        """
        Fetch a specific ligand record from ligand expo
        
        :return: 
        """
        # todo: determine whether this is okay after assuming all PDBs are cleaned and determined usable...
        # Grab the first ligand from the pdb

        # I have determined this is definitely not okay...
        # Case in point: http://ligand-expo.rcsb.org/files/4/4MZ/ipdb/1keq_4MZ_1_A_350__F__D_.ipdb
        # It is empty. What. Why.

        reject = False

        # Go through all ligands bound to PDBs
        pdbid = os.path.basename(os.path.normpath(pdb_file)).split('.')[0]
        for item in self.ligexpo_list:

            # For my specific ligand of interest bound to the correct PDB...
            if pdbid.lower() in item:

                # Report
                lig_request_suffix = item
                print("\n Checking: {}".format(lig_request_suffix))

                # Get ligand record from LigExpo
                time.sleep(0.1)
                full_lig_request = requests.get('http://ligand-expo.rcsb.org/files/{}/{}/ipdb/{}'.format(self.ligand[0], self.ligand, lig_request_suffix),
                                                stream=False).text
                chain = lig_request_suffix.split('_')[3]

                print(full_lig_request)

                # Determine whether there are multiple occupancies for the target ligand
                # I am going to avoid using ligands with multiple occupancies as this suggests binding conformations
                # with low affinity and specificity
                cleaned_target, multiple_occupancies = self.clean_ligand_HETATM_records(full_lig_request)

                print('Multiple Occupancies')
                print(multiple_occupancies)
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

                        # Verify whether I can actually open it with Prody...
                        try:
                            prody.parsePDBStream(io.StringIO(full_lig_request))
                            return reject, full_lig_request, chain

                        except Exception as e:
                            print(e)
                            continue

        # If I can't find any ligands fully represented, reject this PDB
        print('\n\n{} REJECTED! SAD!\n\n'.format(pdbid))
        return True, None, None


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
        # Debugging
        print(target_string)
        print(fragment_pdb)

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
        fragment_target_mapping = [(f_idx, t_idx) for f_idx, t_idx in zip(frag_matches, target_matches)]

        return fragment_target_mapping


    def generate_graph_from_pdb(self, pdb_file):
        """
        Generates a graph representation of a molecule using NetworkX
        
        :param pdb_file: Iterable containing information on molecule connectivity in PDB format. Iterable should either
        be a list or an opened instance of a .pdb file
        :return: 
        """
        pass


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

        print(fragment_target_mapping)

        # For Reference:
        # fragment_target_mapping = [(f_idx, t_idx) for f_idx, t_idx in zip(frag_matches, target_matches)]
        frag_inx_string = [str(a[0]) for a in fragment_target_mapping]
        trgt_inx_string = [str(b[1]) for b in fragment_target_mapping]

        # Debugging
        print(frag_inx_string)
        print(trgt_inx_string)

        # So after inspecting the aligned structures, I've found that the wrong atoms are being aligned...
        # The atoms are currently ordered from highest to lowest according to index... need to change that so that the
        # orders instead agree with the atom mappings
        # I'll need to assemble the matrices by hand it seems.

        frag_atom_selections = [fragment_prody.select('index {}'.format(index)) for index in frag_inx_string]
        trgt_atom_selections = [target_prody.select('index {}'.format(index)) for index in trgt_inx_string]

        print(frag_atom_selections)
        print(trgt_atom_selections)

        frag_atom_coords = np.asarray([atom.getCoords()[0] for atom in frag_atom_selections])
        trgt_atom_coords = np.asarray([atom.getCoords()[0] for atom in trgt_atom_selections])

        print(frag_atom_coords)
        print(trgt_atom_coords)

        # fragment_align_atoms = fragment_prody.select('index ' + ' '.join(frag_inx_string))
        # for atom in fragment_align_atoms:
        #     print(atom)

        return prody.calcTransformation(trgt_atom_coords, frag_atom_coords)


    def clean_ligand_HETATM_records(self, target_string):
        """
        Removes alternate location indicators...
        :param target_string: 
        :return: 
        """
        cleaned_target = ''
        multiple_occupancies = False
        target_string_split = target_string.split('\n')
        for line in target_string_split:
            if line[:6] == 'HETATM':

                if line[16:17] != ' ':
                    multiple_occupancies = True
                    break

                HETATM = line[:6]
                Serial = line[7:11].strip()
                Atom_Name = line[13:16].strip()
                Residue_Name = line[17:20].strip()
                Chain_ID = line[21:22].strip()
                ResSeq_ID = line[23:26].strip()
                x = line[31:38].strip()
                y = line[39:46].strip()
                z = line[47:54].strip()
                Occupancy = line[55:60].strip()
                Temperature = line[61:66].strip()
                Element_Symbol = line[77].strip()

                cleaned_line = '{0:<6} {1:4} {2:>4} {3:3} {4}{5:>4}    {6:>8}{7:>8}{8:>8}{9:>6}{10:>6}         {11}\n'\
                    .format(HETATM, Serial, Atom_Name, Residue_Name, Chain_ID, ResSeq_ID, x, y, z, Occupancy,
                            Temperature, Element_Symbol)
                cleaned_target += cleaned_line

            else:
                cleaned_target += line
                cleaned_target += '\n'

        return cleaned_target, multiple_occupancies



    def apply_transformation(self, transformation_matrix, target_pdb_path, ligand, ligand_chain):
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
        target_pdb = prody.parsePDB(target_pdb_path)
        target_shell = target_pdb.select('within 8 of (resname {} and chain {})'.format(ligand, ligand_chain))

        # todo: implement this in usage so it only has to be done once...
        # Create directory for processed PDBs
        processed_PDBs_path = os.path.join(self.user_defined_dir, 'Transformed_Aligned_PDBs', self.fragment)
        os.makedirs(processed_PDBs_path, exist_ok=True)

        PDBID = os.path.basename(os.path.normpath(target_pdb_path)).split('.')[0]

        transformed_pdb = prody.applyTransformation(transformation_matrix, target_shell)
        prody.writePDB(os.path.join(processed_PDBs_path, '{}_processed.pdb'.format(PDBID)), transformed_pdb)

        print(target_shell)


    ####################################################################################################################
    #
    #  Everything below this block is depreciated or has been replaced!
    #
    ####################################################################################################################


    def extract_atoms_and_connectivities(self, ligand, pdb_file):
        """
        Extract information on atoms and connectivities for a fragment-containing small molecule.

        Biopython is a pain in the ass
        BioPandas doesn't handle CONECT records
        Prody doesn't handle CONECT records

        20170515: So I was writing up my quals proposal and I realized that I could download specific ligands for any 
        structure in the PDB using LigandExpo as outlined here:
        http://ligand-expo.rcsb.org/ld-download.html >> http://ligand-expo.rcsb.org/files/<HASH>/<CC_ID>/<FORMAT>/

        This might be the way to go since it removes the need for all this parsing...

        :param ligand: Three letter code for fragment-containing ligand
        :param pdb_file: Path to the PDB file containing the fragment-containing ligand
        :return: String of HETAM and CONECT records from the input pdb for the given ligand
        """

        # Debugging
        print(pdb_file)

        # todo: remove redundant parsePDB()... don't need to be parsing the same PDB multiple times
        # Determine which chains have my target ligand
        ag = prody.parsePDB(pdb_file)
        ag_ligands = ag.select('resname {}'.format(ligand))

        # Pick a chain and extract the target ligand... arbitrarily picks the first chain for now
        chain = ag_ligands.getChids()[0]

        # Debugging
        print(ag_ligands)
        print(chain)
        for atom in ag_ligands:
            print(atom.getChid())
        print(ligand)
        print(type(ligand))

        # I only need CONECT and HETATM records for my ligand
        # This function assumes that all CONECT records will be at the very end of the PDB file
        pdb = open(pdb_file)
        line_list = []

        # Store HETAM lines for the ligand
        HETAM_num_list = []

        for line in pdb:
            split_line = line.split()

            if split_line[0] == 'HETATM' and split_line[3] == ligand and split_line[4][0] == chain:
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

        # Debugging
        print(ligand_io)

        return ligand_io, chain