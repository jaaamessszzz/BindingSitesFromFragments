#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
import prody
import pprint
import yaml
import json
import re
import itertools
import subprocess
from .alignments import Alignments, Align_PDB
from .utils import *

class Generate_Motif_Residues():
    """
    This is a class for generating and manipulating motif residues for hypothetical ligand binding sites
    
    For each ligand conformer during motif residue generation, I will calculate and output:
    1.  A sparse matrix of motif residues that clash
    2.  Residue-ligand interaction scores as calculated with Rosetta
    3.  Constraint file for each residue-ligand interaction
    """
    def __init__(self, user_defined_dir, motif_cluster_yaml):
        self.user_defined_dir = user_defined_dir
        self.fragment_cluster_list = motif_cluster_yaml
        self.score_interactions_list_path = os.path.join(self.user_defined_dir, 'Inputs', 'Residue_Ligand_Interactions', 'PDBs_to_score.txt')
        self.rosetta_inputs_path = os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs')
        self.residue_ligand_interactions_dir = os.path.join(self.user_defined_dir, 'Inputs', 'Residue_Ligand_Interactions')
        # todo: accommodate C-term residues... could just add 1 to everything and say either or... sloppy though
        self.expected_atoms = {'ALA': 5, 'CYS': 6, 'ASP': 8, 'GLU': 9, 'PHE': 11, 'GLY': 4, 'HIS': 10, 'ILE': 8,
                               'LYS': 9, 'LEU': 8, 'MET': 8, 'ASN': 8, 'PRO': 7, 'GLN': 9, 'ARG': 11, 'SER': 6,
                               'THR': 7, 'VAL': 7, 'TRP': 14, 'TYR': 12}


    def generate_motif_residues(self):
        """
        Define motif residues from clusters outlined in a user defined yaml file for a single ligand conformation.
        A yaml file will also be generated to keep track of which residues clash with each conformer. This way I can 
        output all representative motif residues once and just ignore the ones that clash for each conformer
        * Inputs/User_Inputs/Motif_Clusters.yaml
        """

        # Assemble dict to store list of prody residue instances for each cluster
        fragment_prody_dict = {}

        # For each cluster in the cluster list
        for fragment in self.fragment_cluster_list:
            fragment_prody_dict[fragment] = {}
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

        fragment_output_dir = os.path.join(self.user_defined_dir, 'Representative_Residue_Motifs')
        os.makedirs(fragment_output_dir, exist_ok=True)

        # Start going through fragments and clusters
        total_motif_count = 1
        for fragment in fragment_prody_dict:
            motif_count = 1

            for cluster in fragment_prody_dict[fragment]:
                # For each type of residue in cluster
                residue_types = list(set([residue.getResnames()[0] for residue in fragment_prody_dict[fragment][cluster]]))

                for res in residue_types:
                    # Calculate average residue coordinates
                    cluster_residues = [residue for residue in fragment_prody_dict[fragment][cluster]
                                        if all([residue.getResnames()[0] == res, len(residue) == self.expected_atoms[residue.getResnames()[0]]])]

                    if len(cluster_residues) > 0:
                        average_coordinates = np.mean([a.getCoords() for a in cluster_residues], axis=0)

                        # Select residue closest to average
                        rmsd_list = [prody.calcRMSD(average_coordinates, residue.getCoords()) for residue in cluster_residues]
                        representative_residue = cluster_residues[rmsd_list.index(min(rmsd_list))]

                        # Output residue
                        prody.writePDB(os.path.join(fragment_output_dir, '{4}-{5}-{0}-Cluster_{1}-Motif-{2}-{3}'.format(fragment, cluster, motif_count, res, total_motif_count, len(fragment_prody_dict[fragment][cluster]))),
                                       representative_residue)
                        total_motif_count += 1

                motif_count += 1


    def prepare_motifs_for_conformers(self):
        """
        Prepare things for conformer binding site generation
        :return: 
        """
        os.makedirs(self.residue_ligand_interactions_dir, exist_ok=True)
        # self.generate_residue_ligand_clash_list(os.path.join(self.user_defined_dir, 'Representative_Residue_Motifs'))
        self.generate_residue_residue_clash_matrix()
        self.score_residue_ligand_interactions()

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
        for conformer in pdb_check(os.path.join(target_molecule, 'Inputs','Rosetta_Inputs')):
            target_molecule_prody = prody.parsePDB(conformer).select('not hydrogen')
            clashing_residue_indices = []
            for motif_residue in pdb_check(motif_residue_dir):
                residue_prody = prody.parsePDB(motif_residue).select('not hydrogen')
                if minimum_contact_distance(residue_prody.getCoords(), target_molecule_prody.getCoords()) <= cutoff_distance:
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
        
        # I derped designing the alignment classes so now this just exists. Hi.
        herp_derp = Alignments()

        for conformer in pdb_check(os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs')):
            conformer_name = os.path.basename(os.path.normpath(conformer)).split('.')[0]
            conformer_transformation_dict[conformer_name] = {}

            for fragment in pdb_check(os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs')):
                # todo: somehow define which atoms are absolutely rigid, flexible atoms mess up alignment
                # Map fragment onto conformer
                conformer_pdb_string = open(conformer).read()
                fragment_pdb_string = open(fragment).read()

                # Add fragment prody to fragment_prody_dict
                fragment_name = os.path.basename(os.path.normpath(fragment)).split('.')[0]
                fragment_prody_dict[fragment_name] = prody.parsePDB(fragment)

                # Debugging
                con_name = os.path.basename(os.path.normpath(conformer))
                frag_name = os.path.basename(os.path.normpath(fragment))

                # Determine transformation for fragment onto conformer
                align = Align_PDB(herp_derp,
                                  pdb_file='{}_{}'.format(con_name, frag_name),
                                  target_string=conformer_pdb_string,
                                  fragment_string=fragment_pdb_string)
                align.fragment_target_mapping()

                frag_atom_coords, trgt_atom_coords = align.process_atom_mappings_into_coordinate_sets(
                    align.fragment_target_map)
                transformation = prody.calcTransformation(frag_atom_coords, trgt_atom_coords)
                conformer_transformation_dict[conformer_name][fragment_name] = transformation

                # Debugging: Outputs PDBs of fragment transformed onto corresponding atoms on conformer
                # os.makedirs(os.path.join(self.user_defined_dir, 'Conformer_Tests'), exist_ok=True)
                # transformed_fragment = prody.applyTransformation(rotation, prody.parsePDBStream(io.StringIO(fragment_pdb_string)))
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
        
        # I derped designing the alignment classes so now this just exists. Hi.
        herp_derp = Alignments()
        
        conformer_transformation_dict = self.transform_motif_residues_onto_fragments()
        
        # touch text file to keep track of paths for PDBs to score
        self.score_interactions_list_path = os.path.join(self.user_defined_dir, 'Inputs', 'Residue_Ligand_Interactions', 'PDBs_to_score.txt')
        open(self.score_interactions_list_path, 'w').close()

        # Okay so for the actually residue-residue clashing stuff for each conformer
        # For each conformer I want to determine motif clashes with...
        for conformer in pdb_check(os.path.join(self.user_defined_dir, 'Inputs', 'Rosetta_Inputs')):
            if os.path.basename(os.path.normpath(conformer)) != 'conf.lib.pdb':
                conformer_name = os.path.basename(os.path.normpath(conformer)).split('.')[0]

                # For each residue...
                # Make a list of all transformed prody motif residues, then pass to minimum_contact_distance()
                motif_residue_list = []
                for motif_residue in pdb_check(os.path.join(self.user_defined_dir, 'Representative_Residue_Motifs')):

                    # Parse file name to get fragment, import prody
                    motif_prody = prody.parsePDB(motif_residue)
                    motif_filename = os.path.basename(os.path.normpath(motif_residue)).split('.')[0]
                    motif_source_fragment = motif_filename.split('-')[2]

                    # Map fragment to conformer
                    # def __init__(self, Alignments, pdb_file=None, target_string=None, fragment_path=None):
                    target = open(conformer).read()
                    mobile = os.path.join(self.user_defined_dir, 'Inputs', 'Fragment_Inputs', '{}.pdb'.format(motif_source_fragment))
                    align = Align_PDB(herp_derp, motif_residue, mobile, target)
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
                        if minimum_contact_distance(outer_motif_tuple[1].getCoords(), inner_motif_tuple[1].getCoords()) < clashing_cutoff:

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
        scores_dir = os.path.join(self.residue_ligand_interactions_dir, '{}_scores.txt'.format(current_ligand))
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
        print('Score file outputted as {}!'.format(scores_dir))

        score_df = pd.read_csv(scores_dir,
                               delim_whitespace=True,
                               index_col=4,
                               usecols=['hbond_sc', 'fa_rep', 'fa_elec', 'fa_atr', 'description'],
                               skiprows=0,
                               header=1)
        pprint.pprint(score_df)

    def generate_constraint(self, binding_site_description, prody_binding_site):
        """
        Generate matcher constraint for a given residue-ligand interaction using information from clusters

        :return: 
        """
        # For residue in binding site
        for residue in prody_binding_site:
            pass

            # I need to determine the closest atom-atom contacts and two additional atoms for determining bond torsions and angles
            # Get ideal distance, angle, and torsion values from prody residue
            # Get tolerance from clusters; I'm going to try making the tolerance +/- 1 SD of cluster values
            # Set penalty to whatever since it isn't used by the matcher
            # Set distance to 0, otherwise periodicity (per) to 360
            # Set sample number so that samples are made every 5 degrees or every 0.1A


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
    def __init__(self, user_defined_dir, residue_groups, hypothetical_binding_sites):
        self.user_defined_dir = user_defined_dir
        self.residue_groups = residue_groups
        self.hypothetical_binding_sites = hypothetical_binding_sites


    def calculate_energies_and_rank(self):
        """
        Calculate total binding site interaction energies for all possible conformations of representative binding 
        motifs and rank.
        
        :return: 
        """
        pass


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
                    for pdb in pdb_check(os.path.join(self.user_defined_dir, 'Representative_Residue_Motifs'), base_only=True):
                        if pdb.split('-')[0] == str(residue):
                            binding_residue_list.append(prody.parsePDB(os.path.join(self.user_defined_dir, 'Representative_Residue_Motifs', pdb)))
                            cluster_sum += int(pdb.split('-')[1])

                # Simple distance check for all residues to prevent clashing
                clashing = False
                cutoff_distance = 2

                for outer_index, a in enumerate(binding_residue_list):
                    for inner_index, b in enumerate(binding_residue_list[(outer_index+1):]):
                        if minimum_contact_distance(a.getCoords(), b.getCoords()) < cutoff_distance:
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
                    self.generate_constraints(binding_site_description, complete_binding_site)

                    # Output PDB
                    prody.writePDB(os.path.join(output_path, binding_site_description), complete_binding_site)

    def generate_constraints(self, binding_site_description, complete_binding_site):
        """
        Generate a constraint file for a complete hypothetical binding site by concatenating previously generated
        residue-ligand constraints
        :return: 
        """
        pass