#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
import prody
from rdkit import Chem
import pprint
import yaml
import io
import re
import itertools
import subprocess
import sqlite3
import MySQLdb
import collections
from pathos.multiprocessing import ProcessingPool as Pool
from .alignments import Align_PDB
from .utils import *

class Generate_Motif_Residues():
    """
    This is a class for generating and manipulating motif residues for hypothetical ligand binding sites
    
    For each ligand conformer during motif residue generation, I will calculate and output:
    1.  A sparse matrix of motif residues that clash
    2.  Residue-ligand interaction scores as calculated with Rosetta
    3.  Constraint file for each residue-ligand interaction
    """
    # todo: accommodate C-term residues... could just add 1 to everything and say either or... sloppy though
    expected_atoms = {'ALA': 5, 'CYS': 6, 'ASP': 8, 'GLU': 9, 'PHE': 11, 'GLY': 4, 'HIS': 10, 'ILE': 8,
                      'LYS': 9, 'LEU': 8, 'MET': 8, 'ASN': 8, 'PRO': 7, 'GLN': 9, 'ARG': 11, 'SER': 6,
                      'THR': 7, 'VAL': 7, 'TRP': 14, 'TYR': 12, 'MSE': 8, 'SEC': 6}

    def __init__(self, user_defined_dir, motif_cluster_yaml):
        self.user_defined_dir = user_defined_dir
        self.fragment_cluster_list = motif_cluster_yaml
        self.residue_ligand_interactions_dir = os.path.join(self.user_defined_dir, 'Motifs', 'Residue_Ligand_Interactions')
        self.score_interactions_list_path = os.path.join(self.residue_ligand_interactions_dir, 'PDBs_to_score.txt')
        self.rosetta_inputs_path = os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs')

    def generate_motif_residues(self):
        """
        Define motif residues from clusters outlined in a user defined yaml file for a single ligand conformation.
        A yaml file will also be generated to keep track of which residues clash with each conformer. This way I can 
        output all representative motif residues once and just ignore the ones that clash for each conformer
        * Inputs/User_Inputs/Motif_Clusters.yaml
        """

        # Assemble dict to store list of prody residue instances for each cluster
        fragment_prody_dict = collections.OrderedDict()

        # For each cluster in the cluster list
        for fragment in self.fragment_cluster_list:
            fragment_prody_dict[fragment] = collections.OrderedDict()
            fragment_cluster_path = os.path.join(self.user_defined_dir, 'Cluster_Results', fragment)

            # Make a list for each cluster for a given fragment
            for cluster_number in self.fragment_cluster_list[fragment]:
                fragment_prody_dict[fragment][int(cluster_number)] = []

            # Check pdbs for given fragment and add to appropriate cluster list
            for fragment_pdb in pdb_check(fragment_cluster_path, base_only=True):
                cluster_number = int(re.split('-|_', fragment_pdb)[1])

                if cluster_number in self.fragment_cluster_list[fragment]:
                    fragment_prody_dict[fragment][cluster_number].append(prody.parsePDB(os.path.join(fragment_cluster_path, fragment_pdb)).select('not hydrogen'))

        # Okay, cool. Now that I have all of my clusters neatly organized in a dict, I can go through and generate motif
        # residues from each of the clusters based on residue type

        fragment_output_dir = os.path.join(self.user_defined_dir, 'Motifs', 'Representative_Residue_Motifs')
        os.makedirs(fragment_output_dir, exist_ok=True)

        # Start going through fragments and clusters
        total_motif_count = 1
        for fragment in fragment_prody_dict:
            motif_count = 1

            for cluster in fragment_prody_dict[fragment]:
                # For each type of residue in cluster
                residue_types = sorted(list(set([residue.getResnames()[0] for residue in fragment_prody_dict[fragment][cluster]])))

                for res in residue_types:
                    # Calculate average residue coordinates
                    cluster_residues = [residue for residue in fragment_prody_dict[fragment][cluster]
                                        if all([residue.getResnames()[0] == res, len(residue) == self.expected_atoms[residue.getResnames()[0]]])]

                    if len(cluster_residues) > 0:
                        representative_residue = self._select_respresentative_residue(cluster_residues)
                        residue_type = representative_residue.getResnames()[0]

                        # Output residue
                        prody.writePDB(os.path.join(fragment_output_dir, '{4}-{5}-{0}-Cluster_{1}-Motif_{2}-{3}'.format(fragment, cluster, motif_count, residue_type, total_motif_count, len(fragment_prody_dict[fragment][cluster]))),
                                       representative_residue)
                        total_motif_count += 1

                motif_count += 1

    def _select_respresentative_residue(self, cluster_residues):
        """
        Select representative residues from a given cluster 
        The method for determining representative residues depends on the contact type and residue ID
        * Work in progress *
        
        :param cluster_residues: list of residues in a cluster
        :return: representative residue prody
        """
        average_coordinates = np.mean([a.getCoords() for a in cluster_residues], axis=0)

        # Select residue closest to average
        rmsd_list = [prody.calcRMSD(average_coordinates, residue.getCoords()) for residue in cluster_residues]
        representative_residue = cluster_residues[rmsd_list.index(min(rmsd_list))]


        # Check for MSE and SEC
        # Change MSE into MET, SEC into CYS
        representative_residue_resname = representative_residue.getResnames()[0]
        if representative_residue_resname == 'MSE' or representative_residue_resname == 'SEC':
            representative_residue = self._fix_mse_sec(representative_residue, representative_residue_resname)

        return representative_residue

    def _fix_mse_sec(self, representative_residue, resname):
        """
        Takes MSE/SEC and turns it into MET/CYS
        :param representative_residue: representative_residue with Se
        :return: representative_residue with S
        """
        # Find index of SE
        res_elements = representative_residue.getElements()
        seleno_index = [e for e in res_elements].index('SE')

        # Set SE to S
        res_elements[seleno_index] = 'S'

        # Set elements to MET
        representative_residue.setElements(res_elements)

        # Set resnames MSE>MET, SEC>CYS
        # Set seleno atom name to met/cys atom sulfur atom name
        if resname == 'MSE':
            representative_residue.setResnames(['MET'] * len(representative_residue))
            res_atom_names = representative_residue.getNames()
            res_atom_names[seleno_index] = 'SD'
            representative_residue.setNames(res_atom_names)
        elif resname == 'SEC':
            representative_residue.setResnames(['CYS'] * len(representative_residue))
            res_atom_names = representative_residue.getNames()
            res_atom_names[seleno_index] = 'SG'
            representative_residue.setNames(res_atom_names)
        else:
            raise Exception('Either MSE or SEC had to be set!')

        return representative_residue

    def prepare_motifs_for_conformers(self):
        """
        Prepare things for conformer binding site generation
        :return: 
        """
        os.makedirs(self.residue_ligand_interactions_dir, exist_ok=True)
        # self.generate_residue_ligand_clash_list(os.path.join(self.user_defined_dir, 'Motifs', 'Representative_Residue_Motifs'))
        # self.generate_residue_residue_clash_matrix()
        # self.score_residue_ligand_interactions()
        self.generate_residue_ligand_constraints(torsion_constraint_sample_number=1, angle_constraint_sample_number=1, distance_constraint_sample_number=0)

    def generate_residue_ligand_clash_list(self, motif_residue_dir, cutoff_distance=2):
        """
        Filters representative motif residues that are clashing with the ligand for each conformer in Inputs/Rosetta_Inputs
        Outputs a yaml file for each conformer with a list of motif residue indicies that clash
        Clashing is determined by a distance cutoff for the closest atom-atom distance between each residue and the ligand
        
        :param motif_residue_dir: directory with representative motif residue PDBs
        :param clashing_cutoff: distance cutoff for clashing residues in angstroms. Default set to 2.
        :return: 
        """

        # Generate and output a yaml file for each conformer and list which residues clash
        # This way I can output all representative motif residues once and just ignore the ones that clash for each conformer
        # YAML will be a dict with conformer pdbs as keys and list of residue indices as values

        residue_ligand_clash_dict = {}

        target_molecule = os.path.normpath(self.user_defined_dir)
        for conformer in pdb_check(os.path.join(target_molecule, 'Inputs', 'Rosetta_Inputs'), conformer_check=True):
            target_molecule_prody = prody.parsePDB(conformer).select('not hydrogen')
            clashing_residue_indices = []
            for motif_residue in pdb_check(motif_residue_dir):
                residue_prody = prody.parsePDB(motif_residue).select('not hydrogen')
                if minimum_contact_distance(residue_prody, target_molecule_prody) <= cutoff_distance:
                    residue_index = os.path.basename(motif_residue).split('-')[0]
                    clashing_residue_indices.append(residue_index)
            residue_ligand_clash_dict[os.path.basename(conformer)] = clashing_residue_indices

        yaml.dump(residue_ligand_clash_dict, open(os.path.join(target_molecule, 'Inputs', 'User_Inputs', 'Residue_Ligand_Clash_List.yml'), 'w'))

    def transform_motif_residues_onto_fragments(self):
        """
        Exactly what is says!
        :return: 
        """
        # Precalculate transformation matrices for each fragment for each conformer
        conformer_transformation_dict = {}

        # Open fragment prodys and stick them in a dicitonary
        fragment_prody_dict = {}

        # First need to do translation and rotation of all motif residues relative to fragments in new conformations
        for conformer in pdb_check(os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs'), conformer_check=True):
            conformer_name = os.path.basename(os.path.normpath(conformer)).split('.')[0]
            conformer_transformation_dict[conformer_name] = {}

            for fragment in pdb_check(os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs')):
                # Map fragment onto conformer
                conformer_pdb_string = open(conformer).read()
                fragment_pdb_string = open(fragment).read()

                # Add fragment prody to fragment_prody_dict
                fragment_name = os.path.basename(os.path.normpath(fragment)).split('.')[0]
                fragment_prody_dict[fragment_name] = prody.parsePDB(fragment)

                con_name = os.path.basename(os.path.normpath(conformer))
                frag_name = os.path.basename(os.path.normpath(fragment))

                # Determine transformation for fragment onto conformer
                align = Align_PDB(self.user_defined_dir,
                                  pdb_file='{}_{}'.format(con_name, frag_name),
                                  target_string=conformer_pdb_string,
                                  fragment_string=fragment_pdb_string)
                align.fragment_target_mapping()

                frag_atom_coords, trgt_atom_coords = align.process_atom_mappings_into_coordinate_sets(align.fragment_target_map)

                # Align only on rigid atoms (if they are defined in Rigid_Fragment_Atoms dir)
                frag_inputs_dir = os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs', 'Rigid_Fragment_Atoms')

                if os.path.exists(frag_inputs_dir):
                    frag_rigid_pdb_name = '{}-rigid.pdb'.format(fragment_name)

                    if frag_rigid_pdb_name in os.listdir(frag_inputs_dir):
                        frag_atom_rigid, trgt_atom_rigid = align.return_rigid_atoms(fragment_name, frag_atom_coords, trgt_atom_coords)

                        transformation = prody.calcTransformation(frag_atom_rigid, trgt_atom_rigid)

                    else:
                        transformation = prody.calcTransformation(frag_atom_coords, trgt_atom_coords)

                else:
                    transformation = prody.calcTransformation(frag_atom_coords, trgt_atom_coords)

                conformer_transformation_dict[conformer_name][fragment_name] = transformation

                # Debugging: Outputs PDBs of fragment transformed onto corresponding atoms on conformer
                # import io
                # os.makedirs(os.path.join(self.user_defined_dir, 'Conformer_Tests'), exist_ok=True)
                # transformed_fragment = prody.applyTransformation(transformation, prody.parsePDBStream(io.StringIO(fragment_pdb_string)))
                # prody.writePDB(os.path.join(self.user_defined_dir, 'Conformer_Tests', '{}-{}.pdb'.format(con_name, frag_name)), transformed_fragment)

        return conformer_transformation_dict

    def generate_residue_residue_clash_matrix(self, clashing_cutoff=2):
        """
        Generates a sparse matrix for all representative motif residues that clash with each other. This is determined
        based on the distance of the closest atom-atom interaction between two residues.
        Matrix is output as a .csv that can be imported with numpy.
        
        :param motif_residue_dir: directory with representative motif residue PDBs
        :param clashing_cutoff: distance cutoff for clashing residues in angstroms. Default set to 2.
        :return: 
        """
        
        residue_residue_clash_dict = {}

        conformer_transformation_dict = self.transform_motif_residues_onto_fragments()
        
        # touch text file to keep track of paths for PDBs to score
        open(self.score_interactions_list_path, 'w').close()

        # Okay so for the actually residue-residue clashing stuff for each conformer
        # For each conformer I want to determine motif clashes with...
        for conformer in pdb_check(os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs'), conformer_check=True):
            conformer_name = os.path.basename(os.path.normpath(conformer)).split('.')[0]

            # For each residue...
            # Make a list of all transformed prody motif residues, then pass to minimum_contact_distance()
            motif_residue_list = []
            for motif_residue in pdb_check(os.path.join(self.user_defined_dir, 'Motifs', 'Representative_Residue_Motifs')):

                # Parse file name to get fragment, import prody
                motif_prody = prody.parsePDB(motif_residue)
                motif_filename = os.path.basename(os.path.normpath(motif_residue)).split('.')[0]
                motif_source_fragment = motif_filename.split('-')[2]

                # Map fragment to conformer
                target = open(conformer).read()
                mobile = os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs', '{}.pdb'.format(motif_source_fragment))

                # def __init__(self, Alignments, pdb_file=None, target_string=None, fragment_path=None):
                align = Align_PDB(self.user_defined_dir,
                                  pdb_file=motif_residue,
                                  target_string=open(mobile).read(),
                                  fragment_string=target
                                  )

                align.fragment_target_mapping()

                # Get translation and rotation for fragment onto conformer
                transformation_matrix = conformer_transformation_dict[conformer_name][motif_source_fragment]

                # transformation_matrix.setTranslation(t_vector_end)
                transformed_motif = prody.applyTransformation(transformation_matrix, motif_prody)

                motif_residue_list.append((motif_filename, transformed_motif))

                # # Debugging: Outputs PDBs of motif transformed relative to corresponding fragment atoms on conformer
                # motif_residue_name = motif_filename.split('.')[0]
                # conformer_residues_path = os.path.join(self.user_defined_dir, 'Hypothetical_Binding_Sites', conformer_name)
                # os.makedirs(conformer_residues_path, exist_ok=True)
                # prody.writePDB(os.path.join(conformer_residues_path, ('-'.join([motif_residue_name, conformer_name]) + '.pdb')), transformed_motif)

            # Generate residue-conformer PDBs for scoring here so I don't have to regenerate the transformed residues...
            self.generate_residue_ligand_pdbs(conformer, motif_residue_list)

            residue_residue_clash_set = set()
            for outer_index, outer_motif_tuple in enumerate(motif_residue_list):
                for inner_index, inner_motif_tuple in enumerate(motif_residue_list[outer_index + 1:]):
                    outer_motif_index = outer_motif_tuple[0].split('-')[0]
                    inner_motif_index = inner_motif_tuple[0].split('-')[0]
                    if minimum_contact_distance(outer_motif_tuple[1], inner_motif_tuple[1]) < clashing_cutoff:

                        # If clashing, append tuples in both orders... look up times? Whatever!
                        if outer_motif_index != inner_motif_index:
                            residue_residue_clash_set.add((outer_motif_index, inner_motif_index))
                            residue_residue_clash_set.add((inner_motif_index, outer_motif_index))

            residue_residue_clash_dict[conformer_name] = list(residue_residue_clash_set)

            yaml.dump(residue_residue_clash_dict, open(os.path.join(self.user_defined_dir, 'Inputs', 'User_Inputs', 'Residue_Residue_Clash_COO.yml'), 'w'))

    def generate_residue_ligand_pdbs(self, conformer, motif_residue_list):
        """
        Generate residue-ligand PDBs used for scoring
        :param conformer: path to pdb of conformer
        :param motif_residue_list: list of (motif_filename, motif_prody) tuples
        :return: 
        """
        
        # Write all PDBs to a text file so I can score them all with score_jd2 -in:file:l <text_file_with_list.txt>
        os.makedirs(self.residue_ligand_interactions_dir, exist_ok=True)

        ligand_ID = os.path.basename(os.path.normpath(conformer)).split('.')[0]
        ligand_prody = prody.parsePDB(conformer)

        conformer_residue_output_dir = os.path.join(self.residue_ligand_interactions_dir, ligand_ID)
        os.makedirs(conformer_residue_output_dir, exist_ok=True)

        # For each motif residue in Representative_Residue_Motifs directory, combine ligand and residue
        for motif_residue_name, residue_prody in motif_residue_list:
            residue_ligand_prody = residue_prody + ligand_prody
            motif_residue_index = motif_residue_name.split('-')[0]

            # Output ligand-residue pairs to a new Resdiue_Ligand_Interactions directory under Inputs
            pdb_output_path = os.path.join(conformer_residue_output_dir, '{}-{}.pdb'.format(ligand_ID, motif_residue_index))
            prody.writePDB(pdb_output_path, residue_ligand_prody)
            with open(self.score_interactions_list_path, 'a+') as pdb_list:
                pdb_list.write('{}\n'.format(pdb_output_path))

    def score_residue_ligand_interactions(self):
        """
        Score each unique residue-ligand interaction with Rosetta and export to .csv (fa_atr, hbond_sc, fa_elec)
        
        :return: 
        """
        # Calculate scores for all PDBs in this new directory (batch process, don't forget .params)
        current_ligand = os.path.basename(os.path.normpath(self.user_defined_dir))
        scores_dir = os.path.join(self.residue_ligand_interactions_dir, '{}_scores_raw.txt'.format(current_ligand))
        run_jd2_score = subprocess.Popen(['/Users/jameslucas/Rosetta/main/source/bin/score_jd2.macosclangrelease',
                                          '-l',
                                          self.score_interactions_list_path,
                                          '-extra_res_fa',
                                          os.path.join(self.rosetta_inputs_path, '{}.params'.format(current_ligand)),
                                          # '-scorefile_format',
                                          # 'json',
                                          '-out:file:scorefile',
                                          scores_dir,
                                          '-out:file:score_only',
                                          '-overwrite'
                                          ])
        run_jd2_score.wait()

        # Process score file into a .csv that can be imported by numpy as a matrix
        print('Raw score file outputted as {}!'.format(scores_dir))

        score_df = pd.read_csv(scores_dir,
                               delim_whitespace=True,
                               index_col=4,
                               usecols=['hbond_sc',
                                        'fa_elec',
                                        'fa_atr',
                                        'fa_rep', # Commented out for now so I can just do a quick aggregate for total scores
                                        'description'],
                               skiprows=0,
                               header=1)

        # Output only necessary scores to a .csv
        score_df.to_csv(os.path.join(self.residue_ligand_interactions_dir, '{}_scores_df.csv'.format(current_ligand)))

    def generate_residue_ligand_constraints(self, distance_tolerance_d=0.5, angle_A_tolerance_d=10, angle_B_tolerance_d=10, torsion_A_tolerance_d=10, torsion_AB_tolerance_d=10, torsion_B_tolerance_d=10, torsion_constraint_sample_number=2, angle_constraint_sample_number=2, distance_constraint_sample_number=0):
        """
        Generate matcher constraint for a given residue-ligand interaction using information from clusters
        res1 = residue
        res2 = ligand
        :return: 
        """
        # todo: update to generate constraints for all conformers
        ligand = os.path.basename(os.path.normpath(self.user_defined_dir))

        # Get dict of all representative residue filenames
        representative_residue_dir = os.path.join(self.user_defined_dir, 'Motifs', 'Representative_Residue_Motifs')
        representative_residue_info = {filename.split('-')[0]: filename for filename in pdb_check(representative_residue_dir, base_only=True)}

        # For each directory in Residue_Ligand_Interactions
        for source_conformer in directory_check(self.residue_ligand_interactions_dir):
            print('\n\n{}'.format(source_conformer))

            # For each pdb of a transformed residue-conformer interaction
            for conformer_residue_pdb in pdb_check(source_conformer, base_only=True):
                print(conformer_residue_pdb)

                residue_index = re.split('-|\.|_', conformer_residue_pdb)[2]

                # Get residue prody
                conformer_residue_prody = prody.parsePDB(os.path.join(source_conformer, conformer_residue_pdb))
                residue_prody = conformer_residue_prody.select('protein')

                # Get conformer prody
                # todo: fix conformer_io indicies
                # Could have done a selection from conformer_residue_prody, but there's something off about the indexing
                # So I need to reset the indexing on the selection, the selection keeps indices from original object
                # fragment_source_conformer = conformer_residue_prody.select('resname {}'.format(ligand))
                fragment_source_conformer_path = os.path.join(self.user_defined_dir,
                                                              'Inputs',
                                                              'Rosetta_Inputs',
                                                              '{}.pdb'.format(os.path.basename(os.path.normpath(source_conformer)))
                                                              )
                fragment_source_conformer = prody.parsePDB(fragment_source_conformer_path)

                # Write residue_prody to StringIO
                residue_io = io.StringIO()
                prody.writePDBStream(residue_io, residue_prody)

                # Write conformer to prody StringIO
                # conformer_io = io.StringIO()
                # prody.writePDBStream(conformer_io, fragment_source_conformer)
                # print(conformer_io.getvalue())

                # todo: update how I get this information
                # Get information about residue source from filename... uh...
                # This is an artifact from updating to support conformers...
                residue_split = re.split('-|\.', representative_residue_info[residue_index])
                fragment = residue_split[2]
                cluster = residue_split[3]
                resname = residue_split[5]

                RD_residue = Chem.MolFromPDBBlock(residue_io.getvalue(), removeHs=False)
                RD_ligand = Chem.MolFromPDBBlock(open(fragment_source_conformer_path, 'r').read(), removeHs=False) # Trouble?

                residue_index_atom_map = {atom.getIndex(): atom.getName() for atom in residue_prody.select('not hydrogen')}
                residue_atom_index_map = {v: k for k, v in residue_index_atom_map.items()}
                ligand_index_atom_map = {atom.getIndex(): atom.getName() for atom in fragment_source_conformer.select('not hydrogen')}

                # I need to determine the closest atom-atom contacts and two additional atoms for determining bond torsions and angles
                # NOTE: Contact distance and indicies are for residue and ligand with hydrogens stripped!

                # So it occurred to me that contraints would be more meaningful if the residue was constrained to atoms in
                # the ligand from its source fragment...

                # Map fragment onto conformer
                # todo: fix conformer_io indicies
                # Trying to use conformer_io here causes issues with mapping due to indicies being off...
                align = Align_PDB(self.user_defined_dir,
                                  fragment_string=open(os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs', '{}.pdb'.format(fragment)), 'r').read(),
                                  target_string=open(fragment_source_conformer_path, 'r').read())

                align.fragment_target_mapping()

                conformer_fragment_indices = [pair[1] for pair in align.fragment_target_map]
                conformer_fragment_atoms = fragment_source_conformer.select('index {}'.format(' '.join([str(a) for a in conformer_fragment_indices])))

                contact_distace, residue_index_low, ligand_index_low= minimum_contact_distance(residue_prody, conformer_fragment_atoms, return_indices=True)

                import pprint
                pprint.pprint([atom for atom in residue_prody])
                pprint.pprint([atom for atom in conformer_fragment_atoms])

                # So I derped, residue_index_low and ligand_index_low are indices from the distance matrix and do not
                # necessarily correspond to residue/ligand atom indicies... just happened to work for residues since most
                # of them have hydrogens stripped already

                residue_atom_list = [atom for atom in residue_prody.select('not hydrogen')]
                ligand_atom_list = [atom for atom in conformer_fragment_atoms.select('not hydrogen')]

                residue_contact_atom = residue_atom_list[residue_index_low]
                ligand_contact_atom = ligand_atom_list[ligand_index_low]

                residue_first_atom = residue_contact_atom.getIndex()
                ligand_first_atom = ligand_contact_atom.getIndex()

                # Okay, now the hard part: selecting two additional atoms downstream from the contact atoms for meaningful constraints... programmatically...
                # Using RDKit to determine neighboring atoms, removes Hydrogens by default

                # RESIDUE CONSTRAINT ATOM INDICES
                # Special cases: O (O>C>CA), N (N>CA>C), C (C>CA>N), CA (CA>C>O)
                if residue_contact_atom.getName() in ['C', 'CA', 'CB', 'N', 'O']:
                    if residue_contact_atom.getName()[0] == 'CB':
                        residue_second_atom = residue_atom_index_map['CA']
                        residue_third_atom = residue_atom_index_map['C']
                    elif residue_contact_atom.getName()[0] == 'O':
                        residue_second_atom = residue_atom_index_map['C']
                        residue_third_atom = residue_atom_index_map['CA']
                    elif residue_contact_atom.getName()[0] == 'N':
                        residue_second_atom = residue_atom_index_map['CA']
                        residue_third_atom = residue_atom_index_map['C']
                    elif residue_contact_atom.getName()[0] == 'C':
                        residue_second_atom = residue_atom_index_map['CA']
                        residue_third_atom = residue_atom_index_map['CB']
                    elif residue_contact_atom.getName()[0] == 'CA':
                        residue_second_atom = residue_atom_index_map['C']
                        residue_third_atom = residue_atom_index_map['O']
                    else:
                        raise Exception('wut.')

                else:
                    residue_second_atom = self._determine_next_residue_constraint_atom(residue_contact_atom.getIndex(), RD_residue, residue_prody)
                    residue_third_atom = self._determine_next_residue_constraint_atom(residue_second_atom, RD_residue, residue_prody)

                # So I derped, residue_index_low and ligand_index_low are indices from the distance matrix and do not
                # necessarily correspond to residue/ligand atom indicies... just happened to work for residues since most
                # of them have hydrogens stripped already

                print('RESIDUE - FIRST ATOM: {} {}'.format(residue_contact_atom.getIndex(), residue_index_atom_map[residue_contact_atom.getIndex()]))
                print('RESIDUE - SECOND ATOM: {} {}'.format(residue_second_atom, residue_index_atom_map[residue_second_atom]))
                print('RESIDUE - THIRD ATOM: {} {}'.format(residue_third_atom, residue_index_atom_map[residue_third_atom]))

                # LIGAND CONSTRAINT ATOM INDICES
                # So as far as choosing ligand constraints go, I think it's safe to chose whatever atoms as long as the
                # terminal atom has >1 neighbor. The rationale being this will propogate atom selection toward less
                # flexible parts of the ligand... hopefully...

                # todo: may need to generate multiple ligand constraints depending on where rotatable bonds lie...
                # This shouldn't be a huge issue since I generate all of the transformed ligand-residue contacts
                # Just need to consider the conformer when I concat complete constraint files together

                ligand_second_atom, ligand_third_atom = self._determine_ligand_constraint_atoms(ligand_contact_atom.getIndex(), RD_ligand, fragment_source_conformer)

                print('LIGAND - FIRST ATOM: {} {}'.format(ligand_contact_atom.getIndex(), ligand_index_atom_map[ligand_contact_atom.getIndex()]))
                print('LIGAND - SECOND ATOM: {} {}'.format(ligand_second_atom, ligand_index_atom_map[ligand_second_atom]))
                print('LIGAND - THIRD ATOM: {} {}'.format(ligand_third_atom, ligand_index_atom_map[ligand_third_atom]))

                # Import residues from clusters
                cluster_residue_prody_list = []
                fragment_cluster_path = os.path.join(self.user_defined_dir, 'Cluster_Results', fragment)

                for cluster_residue in pdb_check(fragment_cluster_path, base_only=True):
                    if cluster_residue.split('-')[0] == cluster:
                        cluster_residue_prody = prody.parsePDB(os.path.join(fragment_cluster_path, cluster_residue))
                        residue_type = cluster_residue_prody.getResnames()[0]

                        # Need to filter for residues with missing atoms...
                        if residue_type == resname and cluster_residue_prody.numAtoms() == self.expected_atoms[residue_type]:
                            cluster_residue_prody_list.append(cluster_residue_prody)

                #########################################################################################
                # Get ideal distance, angle, and torsion values from representative motif residue prody #
                #########################################################################################

                ideal_distance = float(contact_distace)
                # 'angle_A' is the angle Res1:Atom2 - Res1:Atom1 - Res2:Atom1
                ideal_angle_A = prody.calcAngle(fragment_source_conformer.select('index {}'.format(ligand_second_atom)),
                                                fragment_source_conformer.select('index {}'.format(ligand_first_atom)),
                                                residue_prody.select('index {}'.format(residue_first_atom))
                                                )
                # 'angle_B' is the angle Res1:Atom1 - Res2:Atom1 - Res2:Atom2
                ideal_angle_B = prody.calcAngle(fragment_source_conformer.select('index {}'.format(ligand_first_atom)),
                                                residue_prody.select('index {}'.format(residue_first_atom)),
                                                residue_prody.select('index {}'.format(residue_second_atom))
                                                )
                # 'torsion_A' is the dihedral Res1:Atom3 - Res1:Atom2 - Res1:Atom1 - Res2:Atom1
                ideal_torsion_A = prody.calcDihedral(fragment_source_conformer.select('index {}'.format(ligand_third_atom)),
                                                     fragment_source_conformer.select('index {}'.format(ligand_second_atom)),
                                                     fragment_source_conformer.select('index {}'.format(ligand_first_atom)),
                                                     residue_prody.select('index {}'.format(residue_first_atom))
                                                     )
                # 'torsion_AB' is the dihedral Res1:Atom2 - Res1:Atom1 - Res2:Atom1 - Res2:Atom2
                ideal_torsion_AB = prody.calcDihedral(fragment_source_conformer.select('index {}'.format(ligand_second_atom)),
                                                      fragment_source_conformer.select('index {}'.format(ligand_first_atom)),
                                                      residue_prody.select('index {}'.format(residue_first_atom)),
                                                      residue_prody.select('index {}'.format(residue_second_atom))
                                                      )
                ideal_torsion_B = prody.calcDihedral(fragment_source_conformer.select('index {}'.format(ligand_first_atom)),
                                                     residue_prody.select('index {}'.format(residue_first_atom)),
                                                     residue_prody.select('index {}'.format(residue_second_atom)),
                                                     residue_prody.select('index {}'.format(residue_third_atom))
                                                     )

                #################################################################################################
                # Get tolerance from clusters; I'm going to try making the tolerance +/- 1 SD of cluster values #
                #################################################################################################

                if len(cluster_residue_prody_list) > 3:
                    distance_list = [prody.calcDistance(fragment_source_conformer.select('index {}'.format(ligand_first_atom)),
                                                        motif.select('index {}'.format(residue_first_atom)))
                                     for motif in cluster_residue_prody_list]
                    distance_tolerance_SD = np.std(distance_list)

                    # 'angle_A' is the angle Res1:Atom2 - Res1:Atom1 - Res2:Atom1
                    angle_A_list = [prody.calcAngle(fragment_source_conformer.select('index {}'.format(ligand_second_atom)),
                                                    fragment_source_conformer.select('index {}'.format(ligand_first_atom)),
                                                    motif.select('index {}'.format(residue_first_atom)))
                                    for motif in cluster_residue_prody_list]
                    angle_A_tolerance_SD = np.std(angle_A_list)

                    # 'angle_B' is the angle Res1:Atom1 - Res2:Atom1 - Res2:Atom2
                    angle_B_list = [prody.calcAngle(fragment_source_conformer.select('index {}'.format(ligand_first_atom)),
                                                    motif.select('index {}'.format(residue_first_atom)),
                                                    motif.select('index {}'.format(residue_second_atom)))
                                    for motif in cluster_residue_prody_list]
                    angle_B_tolerance_SD = np.std(angle_B_list)

                    # 'torsion_A' is the dihedral Res1:Atom3 - Res1:Atom2 - Res1:Atom1 - Res2:Atom1
                    torsion_A_list = [prody.calcDihedral(fragment_source_conformer.select('index {}'.format(ligand_third_atom)),
                                                         fragment_source_conformer.select('index {}'.format(ligand_second_atom)),
                                                         fragment_source_conformer.select('index {}'.format(ligand_first_atom)),
                                                         motif.select('index {}'.format(residue_first_atom)),
                                                         ) for motif in cluster_residue_prody_list]
                    torsion_A_tolerance_SD = np.std(torsion_A_list)

                    # 'torsion_AB' is the dihedral Res1:Atom2 - Res1:Atom1 - Res2:Atom1 - Res2:Atom2
                    torsion_AB_list = [prody.calcDihedral(fragment_source_conformer.select('index {}'.format(ligand_second_atom)),
                                                          fragment_source_conformer.select('index {}'.format(ligand_first_atom)),
                                                          motif.select('index {}'.format(residue_first_atom)),
                                                          motif.select('index {}'.format(residue_second_atom)),
                                                         ) for motif in cluster_residue_prody_list]
                    torsion_AB_tolerance_SD = np.std(torsion_AB_list)

                    # 'torsion_B' is the dihedral Res1:Atom1 - Res2:Atom1 - Res2:Atom2 - Res2:Atom3
                    torsion_B_list = [prody.calcDihedral(fragment_source_conformer.select('index {}'.format(ligand_first_atom)),
                                                         motif.select('index {}'.format(residue_first_atom)),
                                                         motif.select('index {}'.format(residue_second_atom)),
                                                         motif.select('index {}'.format(residue_third_atom)),
                                                          ) for motif in cluster_residue_prody_list]
                    torsion_B_tolerance_SD = np.std(torsion_B_list)

                    # Set lower and upper bounds for tolerances
                    distance_tolerance = 0.5 if distance_tolerance_SD > 0.5 else distance_tolerance_SD
                    angle_A_tolerance = (360 / (2 * angle_constraint_sample_number + 1) * angle_constraint_sample_number) if angle_A_tolerance_SD > 120 else angle_A_tolerance_SD
                    angle_B_tolerance = (360 / (2 * angle_constraint_sample_number + 1) * angle_constraint_sample_number) if angle_B_tolerance_SD > 120 else angle_B_tolerance_SD
                    torsion_A_tolerance = (360 / (2 * torsion_constraint_sample_number + 1) * torsion_constraint_sample_number) if torsion_A_tolerance_SD > 120 else torsion_A_tolerance_SD
                    torsion_B_tolerance = (360 / (2 * torsion_constraint_sample_number + 1) * torsion_constraint_sample_number) if torsion_B_tolerance_SD > 120 else torsion_B_tolerance_SD
                    torsion_AB_tolerance = (360 / (2 * torsion_constraint_sample_number + 1) * torsion_constraint_sample_number) if torsion_AB_tolerance_SD > 120 else torsion_AB_tolerance_SD

                else:
                    distance_tolerance = distance_tolerance_d
                    angle_A_tolerance = angle_A_tolerance_d
                    angle_B_tolerance = angle_B_tolerance_d
                    torsion_A_tolerance = torsion_A_tolerance_d
                    torsion_AB_tolerance = torsion_AB_tolerance_d
                    torsion_B_tolerance = torsion_B_tolerance_d

                residue_resname = residue_prody.getResnames()[0]
                ligand_resname = fragment_source_conformer.getResnames()[0]

                constraint_block = ['  TEMPLATE::   ATOM_MAP: 1 atom_name: {} {} {}'.format(fragment_source_conformer.select('index {}'.format(ligand_contact_atom.getIndex())).getNames()[0],
                                                                                            fragment_source_conformer.select('index {}'.format(ligand_second_atom)).getNames()[0],
                                                                                            fragment_source_conformer.select('index {}'.format(ligand_third_atom)).getNames()[0]
                                                                                            ),
                                    '  TEMPLATE::   ATOM_MAP: 1 residue3: {}\n'.format(ligand_resname),
                                    '  TEMPLATE::   ATOM_MAP: 2 atom_name: {} {} {}'.format(residue_prody.select('index {}'.format(residue_contact_atom.getIndex())).getNames()[0],
                                                                                            residue_prody.select('index {}'.format(residue_second_atom)).getNames()[0],
                                                                                            residue_prody.select('index {}'.format(residue_third_atom)).getNames()[0]
                                                                                            ),
                                    '  TEMPLATE::   ATOM_MAP: 2 residue3: {}\n'.format(residue_resname),
                                    '  CONSTRAINT:: distanceAB: {0:7.2f} {1:6.2f} {2:6.2f}       0  {3:3}'.format(ideal_distance, distance_tolerance, 100, distance_constraint_sample_number),
                                    '  CONSTRAINT::    angle_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(float(ideal_angle_A), angle_A_tolerance, 100, angle_constraint_sample_number),
                                    '  CONSTRAINT::    angle_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(float(ideal_angle_B), angle_B_tolerance, 100, angle_constraint_sample_number),
                                    '  CONSTRAINT::  torsion_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(float(ideal_torsion_A), torsion_A_tolerance, 100, torsion_constraint_sample_number),
                                    '  CONSTRAINT::  torsion_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(float(ideal_torsion_B), torsion_B_tolerance, 100, torsion_constraint_sample_number),
                                    '  CONSTRAINT:: torsion_AB: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(float(ideal_torsion_AB), torsion_AB_tolerance, 100, torsion_constraint_sample_number)
                                    ]

                single_constraint_path = os.path.join(self.user_defined_dir, 'Motifs', 'Single_Constraints')
                os.makedirs(single_constraint_path, exist_ok=True)

                with open(os.path.join(single_constraint_path, '{}-{}.cst'.format(os.path.basename(os.path.normpath(source_conformer)), residue_index)), 'w') as constraint_file:
                    constraint_file.write('\n'.join(constraint_block))

    def _determine_next_residue_constraint_atom(self, current_atom_index, RD_residue, residue_prody):
        """
        Grabs the next residue atom down the line for residue constraint records
        :return: 
        """
        # Okay so I can't figure out for now how to retrieve atom names from RDKit atom objects,,,
        # I'm going to create a dict mapping atom indices to atom names from the prody object
        # Prody and RDKit both seem to derive atom indices in order of the PDB file
        residue_index_atom_map = {atom.getIndex(): atom.getName() for atom in residue_prody.select('not hydrogen')}
        residue_atom_index_map = {v: k for k, v in residue_index_atom_map.items()}
        
        # Setting up hierarchy for residue atom selections - want selections to propagate toward main chain
        # So.... H>Z>E>D>G>B>A
        # Okay so starting at any given atom I should be able to determine connectivity and how many atoms out I am from CA
        atom_hierarchy = ['H', 'Z', 'E', 'D', 'G', 'B', 'A']

        residue_contact_atom_neighbors = [atom.GetIdx() for atom in RD_residue.GetAtomWithIdx(int(current_atom_index)).GetNeighbors() if atom.GetSymbol != 'H']
        if len(residue_contact_atom_neighbors) == 1:
            next_atom_index = residue_contact_atom_neighbors[0]
            return next_atom_index
        
        else:
            first_atom_ID = residue_index_atom_map[current_atom_index]
            second_atom_IDs = [residue_index_atom_map[atom] for atom in residue_contact_atom_neighbors]

            for atom in second_atom_IDs:
                if atom[1] == atom_hierarchy[atom_hierarchy.index(first_atom_ID[1]) + 1]:
                    next_atom_index = residue_atom_index_map[atom]
                    return next_atom_index

    
    def _determine_ligand_constraint_atoms(self, current_atom_index, RD_ligand, ligand_prody):
        """
        Grabs the next ligand atoms down the line for ligand constraint records
        :return: 
        """
        ligand_index_atom_map = {atom.getIndex(): atom.getName() for atom in ligand_prody.select('not hydrogen')}

        # DEBUGGING
        pprint.pprint(ligand_index_atom_map)

        ligand_contact_atom_neighbors = [atom for atom in RD_ligand.GetAtomWithIdx(int(current_atom_index)).GetNeighbors()]

        # Check all atoms that are adjacent to ligand contact atom
        for second_atom in ligand_contact_atom_neighbors:

            # Only investigate potential second constraint atoms with more than one neighbor
            if len([atom.GetIdx() for atom in second_atom.GetNeighbors()]) > 1:
                ligand_second_atom_neighbors = [atom for atom in second_atom.GetNeighbors() if atom.GetSymbol() != 'H']
                second_atom_neighbors_set = set([ligand_index_atom_map[atom.GetIdx()] for atom in ligand_second_atom_neighbors])

                # Check all atoms that are adjacent to the potential second constraint atom
                for third_atom in ligand_second_atom_neighbors:
                    third_atom_neighbors_set = set([ligand_index_atom_map[atom.GetIdx()] for atom in third_atom.GetNeighbors() if atom.GetSymbol() != 'H'])
                    common_two_three = third_atom_neighbors_set & second_atom_neighbors_set

                    # Only investigate potential third constraint atoms with more than one neighbor, not including second atom
                    if len(third_atom_neighbors_set - common_two_three) > 1 and third_atom.GetIdx() != current_atom_index:
                        ligand_second_atom = second_atom.GetIdx()
                        ligand_third_atom = third_atom.GetIdx()

                        return ligand_second_atom, ligand_third_atom

    def define_second_shell_contraints(self):
        """
        Identify second shell contacts with motif residues
        :return: 
        """
        # This is most likely just going to take the PDB for a given representative motif residue and pull out a 10 A
        # shell or something like that... this is really going to be useful if I can use it on super tight clusters
        # so I can pull out second shell interaction patterns
        pass


class Generate_Binding_Sites():
    """
    This is a class for combining specified groups of binding motifs into hypothetical binding sites
    User input required to deal with residues that obviously do not get along with other motif residues or the ligand
    e.g. steric clashing, non-favorable interactions
    """

    ref2015_weights = {'fa_atr': 1,
                       'fa_elec': 1,
                       'hbond_sc': 1,
                       'fa_rep': 0.55
                       }

    def __init__(self, user_defined_dir, residue_groups=None, hypothetical_binding_sites=None):
        self.user_defined_dir = user_defined_dir
        self.residue_groups = residue_groups
        self.hypothetical_binding_sites = hypothetical_binding_sites

        self.constraints_path = os.path.join(self.user_defined_dir, 'Complete_Matcher_Constraints')
        self.binding_site_pdbs = os.path.join(self.constraints_path, 'Binding_Site_PDBs')
        self.complete_constraint_files = os.path.join(self.constraints_path, 'Constraint_Files')

    def calculate_energies_and_rank(self, use_mysql=True):
        """
        Calculate total binding site interaction energies for all possible conformations of representative binding 
        motifs and rank.
        
        Probs a good idea to use a SQLite3 DB for keeping track of all scores and binding site combintations generated
        
        :return: 
        """
        # Retrieve list of Residue-Residue and Residue-Ligand Clashes
        residue_ligand_clash_dict = yaml.load(open(os.path.join(self.user_defined_dir, 'Inputs', 'User_Inputs', 'Residue_Ligand_Clash_List.yml'), 'r'))
        residue_residue_clash_dict = yaml.load(open(os.path.join(self.user_defined_dir, 'Inputs', 'User_Inputs', 'Residue_Residue_Clash_COO.yml'), 'r'))

        # Import Rosetta Score df
        current_ligand = os.path.basename(os.path.normpath(self.user_defined_dir))
        score_df = pd.read_csv(os.path.join(self.user_defined_dir, 'Motifs', 'Residue_Ligand_Interactions', '{}_scores_df.csv'.format(current_ligand)), index_col='description')

        # Set up agg score dict (in the lamest fashion possible)
        score_agg_dict = {}
        for conformer in pdb_check(os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs'), base_only=True):
            if conformer != 'conf.lib.pdb':
                score_agg_dict[conformer.split('.')[0]] = {}

        for index, row in score_df.iterrows():
            current_conformer = index.split('-')[0]
            motif_res_index = re.split('-|_', index)[2]

            score_agg_dict[current_conformer][motif_res_index] = sum([row['fa_atr'] * self.ref2015_weights['fa_atr'],
                                                                      row['fa_elec'] * self.ref2015_weights['fa_elec'],
                                                                      row['hbond_sc'] * self.ref2015_weights['hbond_sc'],
                                                                      row['fa_rep'] * self.ref2015_weights['fa_rep']
                                                                      ])

        if use_mysql:
              self._use_mysql(residue_residue_clash_dict, residue_ligand_clash_dict, score_agg_dict)
        else:
              self._use_sqlite(residue_residue_clash_dict, residue_ligand_clash_dict, score_agg_dict)

    def _generate_sqlite_db(self):
        """
        Generate SQLite database + table for storing binding site scores
        :return: 
        """
        db_path = os.path.join(self.user_defined_dir, 'Inputs', 'Scored_Binding_Motifs-sqlite.db')
        db_existed = os.path.exists(db_path)
        connection = sqlite3.connect(db_path)
        cur = connection.cursor()

        if not db_existed:
            cur.execute("create table binding_motif_scores (conformer TEXT NOT NULL, first INTEGER NOT NULL, second INTEGER NOT NULL, third INTEGER NOT NULL, fourth INTEGER NOT NULL, score REAL NOT NULL, PRIMARY KEY (conformer, first, second, third, fourth))")

        return connection, cur

    def _use_sqlite(self, residue_residue_clash_dict, residue_ligand_clash_dict, score_agg_dict):
        # Generate SQLite DB
        sqlite_connection, sqlite_cusor = self._generate_sqlite_db()

        # List of residue indicies
        rep_motif_path = os.path.join(self.user_defined_dir, 'Representative_Residue_Motifs')
        representative_motif_residue_indices = [motif.split('-')[0] for motif in
                                                pdb_check(rep_motif_path, base_only=True)]

        # For every conformer I've generated
        for conformer, residue_list in residue_residue_clash_dict.items():

            # Get rid of any residues that clash with the ligand straight away
            useable_residues_distance = list(
                set(representative_motif_residue_indices) - set(residue_ligand_clash_dict[conformer + '.pdb']))
            useable_residues = list(set(useable_residues_distance) - set(
                [key for key, value in score_agg_dict[conformer].items() if value > 0]))

            # Generate all XXX choose 4 binding site configurations
            list_of_residue_combinations = itertools.combinations(useable_residues, 4)

            # For each combination of four residues...
            for combo in list_of_residue_combinations:
                # For each pair or residues, check if tuple is in residue_residue_clash_dict
                fail_residue_combo = any([(pair in residue_list) for pair in itertools.combinations(combo, 2)])

                # If there are no clashes, calculate Rosetta score (sum(fa_atr, fa_elec, hbond_sc))
                if not fail_residue_combo:
                    total_score = sum([score_agg_dict[conformer][res] for res in combo])

                    # Push to SQLite DB table if score is < 0
                    if total_score < 0:
                        sorted_combo = sorted(combo)
                        # print('Inserting {:10} {:3} {:3} {:3} {:3} {:7.3f} into SQLite3 DB'.format(conformer, sorted_combo[0], sorted_combo[1], sorted_combo[2], sorted_combo[3], total_score))
                        sqlite_cusor.execute(
                            "INSERT OR IGNORE INTO binding_motif_scores (conformer, first, second, third, fourth, score) VALUES (?,?,?,?,?,?)",
                            (str(conformer), int(sorted_combo[0]), int(sorted_combo[1]), int(sorted_combo[2]),
                             int(sorted_combo[3]), float(total_score)))

            sqlite_connection.commit()

        table = pd.read_sql_query("SELECT * from binding_motif_scores", sqlite_connection)
        table.to_csv("ASDF" + '.csv', index_label='index')
        sqlite_connection.close()

    def _generate_mysql_db(self):
        """
        Generate a MySQL database + table for storing binding site scores
        Switching to MySQL to support multiprocessing
        :return: 
        """
        # Requires initial setup of MySQL database
        # Create ~/.my.cnf with [client] header (username, password)
        # Ensure MySQL permissions are set so that user can read, write, and execute
        connection = MySQLdb.connect(host='localhost', read_default_file="~/.my.cnf")
        connection.query('CREATE DATABASE IF NOT EXISTS scored_binding_motifs_{}'.format(os.path.basename(os.path.normpath(self.user_defined_dir))))
        connection.query('USE scored_binding_motifs_{}'.format(os.path.basename(os.path.normpath(self.user_defined_dir))))
        connection.query('CREATE TABLE IF NOT EXISTS binding_motif_scores (conformer VARCHAR(10) NOT NULL, first INTEGER NOT NULL, second INTEGER NOT NULL, third INTEGER NOT NULL, fourth INTEGER NOT NULL, score REAL NOT NULL, PRIMARY KEY (conformer, first, second, third, fourth))')
        cursor = connection.cursor()
        return connection, cursor

    def _use_mysql(self, residue_residue_clash_dict, residue_ligand_clash_dict, score_agg_dict):
        # Generate SQLite DB
        mysql_connection, mysql_cursor = self._generate_mysql_db()

        # List of residue indicies
        rep_motif_path = os.path.join(self.user_defined_dir, 'Motifs', 'Representative_Residue_Motifs')
        representative_motif_residue_indices = [motif.split('-')[0] for motif in
                                                pdb_check(rep_motif_path, base_only=True)]

        def push_scores_to_db(residue_residue_clash_dict_tuple):
            connection_embed = MySQLdb.connect(host='localhost',
                                               db='scored_binding_motifs_{}'.format(os.path.basename(os.path.normpath(self.user_defined_dir))),
                                               read_default_file="~/.my.cnf")
            cursor_embed = connection_embed.cursor()

            score_list = []

            # For every conformer I've generated...
            conformer = residue_residue_clash_dict_tuple[0]
            residue_list = residue_residue_clash_dict_tuple[1]

            # Get rid of any residues that clash with the ligand straight away
            useable_residues_distance = list(set(representative_motif_residue_indices) - set(residue_ligand_clash_dict[conformer + '.pdb']))
            useable_residues = list(set(useable_residues_distance) - set([key for key, value in score_agg_dict[conformer].items() if value > 0]))

            # Generate all XXX choose 4 binding site configurations
            list_of_residue_combinations = itertools.combinations(useable_residues, 4)

            # For each combination of four residues...
            for combo in list_of_residue_combinations:
                # For each pair or residues, check if tuple is in residue_residue_clash_dict
                fail_residue_combo = any([(pair in residue_list) for pair in itertools.combinations(combo, 2)])

                # If there are no clashes, calculate Rosetta score (sum(fa_atr, fa_elec, hbond_sc))
                if not fail_residue_combo:
                    total_score = sum([score_agg_dict[conformer][res] for res in combo])

                    # Push to SQLite DB table if score is < 0
                    if total_score < 0:
                        sorted_combo = sorted(combo)

                        # Sanitary? No. Good enough for now? Yes.
                        # mysql_connection.query("INSERT OR IGNORE INTO binding_motif_scores (conformer, first, second, third, fourth, score) VALUES ({},{},{},{},{},{})".format(str(conformer), int(sorted_combo[0]), int(sorted_combo[1]), int(sorted_combo[2]), int(sorted_combo[3]), float(total_score)))
                        score_list.append((str(conformer), int(sorted_combo[0]), int(sorted_combo[1]), int(sorted_combo[2]), int(sorted_combo[3]), float(total_score)))

            insert_query = """INSERT IGNORE INTO binding_motif_scores (conformer, first, second, third, fourth, score) VALUES (%s,%s,%s,%s,%s,%s)"""
            cursor_embed.executemany(insert_query, score_list)
            connection_embed.commit()

        process = Pool()
        process.map(push_scores_to_db, residue_residue_clash_dict.items())
        process.close()
        process.join()

        mysql_connection.close()

    def generate_binding_site_constraints(self, score_cutoff=-15, secondary_matching=False, use_mysql=True):
        """
        Generate constraint files (and optionally binding site PDBs) for binding sites that pass score filters
        :return: 
        """
        if use_mysql:
            mysql_connection = MySQLdb.connect(host='localhost',
                                               db='scored_binding_motifs_{}'.format(os.path.basename(os.path.normpath(self.user_defined_dir))),
                                               read_default_file="~/.my.cnf")
            mysql_cursor = mysql_connection.cursor()
            mysql_cursor.execute("SELECT * FROM binding_motif_scores WHERE score < {}".format(score_cutoff))
            score_table_rows = mysql_cursor.fetchall()
        else:
            sqlite_connection, sqlite_cusor = self._generate_sqlite_db()
            sqlite_cusor.execute("SELECT * FROM binding_motif_scores WHERE score < {}".format(score_cutoff))
            score_table_rows = sqlite_cusor.fetchall()

        # Create directories and files
        os.makedirs(self.constraints_path, exist_ok=True)
        os.makedirs(self.binding_site_pdbs, exist_ok=True)
        os.makedirs(self.complete_constraint_files, exist_ok=True)

        # Touch
        open(os.path.join(self.binding_site_pdbs, 'Binding_sites_to_score.txt'), 'w').close()

        if len(score_table_rows) > 0:
            for row in score_table_rows:
                score_things = True
                row_conformer = row[0]
                row_motif_indicies = row[1:-1]
                row_score = row[5]

                # Generate binding site PDB
                self._generate_constraint_file_binding_site(row_conformer, row_motif_indicies)

                # Get individual motif residue constraint blocks
                motif_constraint_block_list = [open(os.path.join(self.user_defined_dir, 'Motifs', 'Single_Constraints', '{}.cst'.format(index))).read() for index in row_motif_indicies]

                # Write constraint blocks to file
                with open(os.path.join(self.complete_constraint_files, '-'.join([str(a) for a in row[:-1]]) + '.cst'), 'w') as complete_constraint_file:

                    # Options for secondary matching
                    for index, block in enumerate(motif_constraint_block_list):
                        complete_constraint_file.write('CST::BEGIN\n')
                        complete_constraint_file.write(motif_constraint_block_list[index])
                        complete_constraint_file.write('\n')
                        if secondary_matching and index != 0:
                            complete_constraint_file.write('  ALGORITHM_INFO:: match\n')
                            complete_constraint_file.write('    SECONDARY_MATCH: DOWNSTREAM\n')
                            complete_constraint_file.write('  ALGORITHM_INFO::END\n')
                        complete_constraint_file.write('CST::END\n')

            # Score complete binding sites
            self._score_constraint_file_binding_site()
            print('Done!')

        else:
            print('There are no motifs with a score less than {}!'.format(score_cutoff))

    def _generate_constraint_file_binding_site(self, conformer, motif_indicies):
        """
        Generates a binding site PDB based on a constraint file
        :param row_conformer: 
        :param row_motif_indicies: 
        :return: 
        """

        #todo: update so that conformers constraints are supported
        conformer_path = os.path.join(self.user_defined_dir, 'Motifs', 'Residue_Ligand_Interactions', conformer)

        # Use residue-ligand pairs in Inputs/Residue_Ligand_Interactions
        relevant_residue_list = [open(os.path.join(conformer_path, '{}-{}.pdb'.format(conformer, str(index)))) for index in motif_indicies]
        conformer_pdb = open(os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs', '{}.pdb'.format(conformer))).read()

        #  Ligand HETATM + Residue ATOM records
        binding_site_pdb_path = os.path.join(self.binding_site_pdbs, '{}-{}.pdb'.format(conformer, '-'.join([str(a) for a in motif_indicies])))
        with open(binding_site_pdb_path, 'w') as binding_site_pdb:
            binding_site_pdb.write('# {} {}\n'.format(conformer, ' '.join([str(a) for a in motif_indicies])))
            for res in relevant_residue_list:
                for line in res:
                    if line.split()[0] == 'ATOM':
                        binding_site_pdb.write(line)
            binding_site_pdb.write(conformer_pdb)

        # Add path to .txt for scoring later
        with open(os.path.join(self.binding_site_pdbs, 'Binding_sites_to_score.txt'), 'a') as score_list:
            score_list.write('{}\n'.format(binding_site_pdb_path))

    def _score_constraint_file_binding_site(self):
        """
        Scores a binding site as defined by a constraint file generated using calculate_energies_and_rank()
        This is to screen for residue-residue clashes that were not caught with the distance filter
        :return: 
        """

        # Calculate scores for all PDBs in this new directory (batch process, don't forget .params)
        current_ligand = os.path.basename(os.path.normpath(self.user_defined_dir))
        scores_dir = os.path.join(self.constraints_path, '{}_scores_raw.txt'.format(current_ligand))
        score_sites_list_path = os.path.join(self.binding_site_pdbs, 'Binding_sites_to_score.txt')

        run_jd2_score = subprocess.Popen(['/Users/jameslucas/Rosetta/main/source/bin/score_jd2.macosclangrelease',
                                          '-l',
                                          score_sites_list_path,
                                          '-extra_res_fa',
                                          os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs', '{}.params'.format(current_ligand)),
                                          '-out:file:scorefile',
                                          scores_dir,
                                          '-out:file:score_only',
                                          '-overwrite'
                                          ])
        run_jd2_score.wait()

        # Process final score file with weighted scores
        print('Raw score file outputted as {}!'.format(scores_dir))

        score_df = pd.read_csv(scores_dir,
                               delim_whitespace=True,
                               index_col=4,
                               usecols=['hbond_sc',
                                        'fa_elec',
                                        'fa_atr',
                                        'fa_rep',
                                        'description'],
                               skiprows=0,
                               header=1)

        for index, row in score_df.iterrows():
            score_df.loc[index:, 'total_score'] = sum([row['fa_atr'] * self.ref2015_weights['fa_atr'],
                                                       row['fa_elec'] * self.ref2015_weights['fa_elec'],
                                                       row['hbond_sc'] * self.ref2015_weights['hbond_sc'],
                                                       row['fa_rep'] * self.ref2015_weights['fa_rep']
                                                       ])

        # Output only necessary scores to a .csv
        score_df.to_csv(os.path.join(self.constraints_path, '{}_scores_df.csv'.format(current_ligand)))

    def generate_binding_sites_by_hand(self):
        """
        This method takes the user defined residue groups and combines them as specified in the hypothetical_binding_sites
        yaml file. This method will output all hypothetical binding sites into a new directory and rank them based on the
        sum total of cluster members from which the representative residues were derived.
        
        {Rank}-{res_#}_{res_#}_{res_#}_{res_#}-{Description}.pdb
        
        :return: 
        """
        # For each hypothetical binding site defined in the yaml file
        for hbs in self.hypothetical_binding_sites:
            # Generate output path
            output_path = os.path.join(self.user_defined_dir, 'Hypothetical_Binding_Sites', hbs)
            os.makedirs(output_path, exist_ok=True)

            # Generate all possible residue combinations for all positions defined
            list_of_residue_combinations = itertools.product(*[self.residue_groups[group] for group in self.hypothetical_binding_sites[hbs]])

            # For each unique binding site
            for residue_combination in list_of_residue_combinations:

                binding_residue_list = []
                # Open all residue pdbs from representative binding residues directory as defined in list
                # Calculate sum of all cluster representitives
                cluster_sum = 0

                for residue in residue_combination:
                    for pdb in pdb_check(os.path.join(self.user_defined_dir, 'Motifs', 'Representative_Residue_Motifs'), base_only=True):
                        if pdb.split('-')[0] == str(residue):
                            binding_residue_list.append(prody.parsePDB(os.path.join(self.user_defined_dir, 'Motifs', 'Representative_Residue_Motifs', pdb)))
                            cluster_sum += int(pdb.split('-')[1])

                # Simple distance check for all residues to prevent clashing
                clashing = False
                cutoff_distance = 2

                for outer_index, a in enumerate(binding_residue_list):
                    for inner_index, b in enumerate(binding_residue_list[(outer_index+1):]):
                        if minimum_contact_distance(a, b) < cutoff_distance:
                            clashing = True
                            print('CLASHING!!!!~!!!!!!!!')
                            break
                    if clashing:
                        break

                # Proceed to generating and outputting binding site if non-clashing
                if not clashing:
                    # Combine all residues into the same file, including the ligand
                    complete_binding_site = prody.parsePDB(os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs', '{}_0001.pdb'.format(os.path.normpath(self.user_defined_dir))))
                    for binding_residue in binding_residue_list:
                        complete_binding_site = complete_binding_site + binding_residue

                    # Generate Constraints
                    binding_site_description = '{0}-{1}.pdb'.format(cluster_sum, '_'.join([str(a) for a in residue_combination]))
                    self.generate_binding_site_constraints(binding_site_description, complete_binding_site)

                    # Output PDB
                    prody.writePDB(os.path.join(output_path, binding_site_description), complete_binding_site)