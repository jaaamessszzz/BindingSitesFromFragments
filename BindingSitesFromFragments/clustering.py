#!/usr/bin/env python3

import os
import sys
import pickle
import operator
import itertools
from pprint import pprint

import prody
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster

from .utils import *

# --- Silence ProDy --- #
prody.confProDy(verbosity='none')

class Cluster(object):
    """
    This class is responsible for taking aligned fragment PDBs and identifying representative contacts. 
    """

    def __init__(self, processed_PDBs_dir):
        self.processed_PDBs_dir = processed_PDBs_dir
        self.pdb_object_list = self._import_pdbs()
        self.distance_cutoff = 1 - np.cos(20 * np.pi / 180)  # Maximum of 20 degree deviation between clusters
        self.clusters = dict()
        self.df = None

    def _import_pdbs(self):
        """
        For each fragment ensemble, converts each residue in all processed PDBs into objects with representative 
        matrices and other relevant information
        
        :return: 
        """
        processsed_residue_list = []
        for pdb in pdb_check(self.processed_PDBs_dir):

            # Make sure I can load things...
            try:

                prody_protein = prody.parsePDB(pdb)

                # Check that residues exist within cutoff distance provided in alignments, otherwise pass
                prody_protein_selection = prody_protein.select('protein and not hetatm')
                if prody_protein_selection == None:
                    continue
                else:
                    prody_protein_hv = prody_protein_selection.getHierView()

            except Exception as e:
                print(e)
                continue

            pdb_info = os.path.basename(os.path.normpath(pdb))
            prody_ligand = prody.parsePDB(pdb).select('hetatm and resname {}'.format(pdb_info.split('_')[1]))

            # todo: CATCH THIS!!!
            if prody_ligand is None: continue

            # Iterate over residues in contacts and generate representative vector with weights applied
            processsed_residue_list += [fragment_PDB(residue, pdb_info, prody_ligand,) for residue in prody_protein_hv.iterResidues()]

        processsed_residue_list_cleaned = [residue for residue in processsed_residue_list if residue.viable is not None]
        print(f'Unique processed and viable residues: {len(processsed_residue_list_cleaned)}')

        return processsed_residue_list_cleaned

    def cluster_scipy(self, display_dendrogram=False):
        """
        Using Hierarchical clustering through scipy to identify representative contacts.

        Clustering is separated into two steps since hierarchical clustering is an O(n^2) thing and blows up my RAM when
        clustering relatively common fragments... First, residues are separated based on categorical variables into a
        dict... think of it as clustering by hamming distance with a distance cutoff of one. Each of these groups are
        then clustered by cosine distance to group residues in space around the fragment.

        :return: 
        self.clusters: list of indicies corresponding to cluster numbers for each element in self.pdb_object_list
        """

        category_dict = dict()
        unique_category_list = list()

        for residue in self.pdb_object_list:
            cat_array_string = residue.categorical_array.tostring()

            if cat_array_string in category_dict.keys():
                category_dict[cat_array_string].append(residue)
            else:
                category_dict[cat_array_string] = [residue]
                unique_category_list.append(np.copy(residue.categorical_array))

        # CLUSTER BY CATEGORIES FIRST!!! ALLOW HAMMING DIST OF  1
        hamming_distance = pdist(unique_category_list, metric='hamming')  # Hamming distance normalized to vector length!
        cat_Z = linkage(hamming_distance, method='complete')
        cat_clusters = fcluster(cat_Z, 2/14, criterion='distance')
        mapped_categories = [(cluster, category) for cluster, category in zip(cat_clusters, unique_category_list)]

        pprint(mapped_categories)

        # Combine lists from category_dict into category-level clusters
        clustered_category_dict = dict()
        for cluster, category_vector in mapped_categories:

            print(cluster, category_vector, len(category_dict[category_vector.tostring()]))
            if cluster not in clustered_category_dict.keys():
                clustered_category_dict[cluster] = list()
            clustered_category_dict[cluster] = clustered_category_dict[cluster] + category_dict[category_vector.tostring()]

        pprint([(key, len(asdf)) for key, asdf in clustered_category_dict.items()])

        # Cluster spatial positions for categorically clustered residues
        for category, category_contacts in clustered_category_dict.items():
            vector_array = [residue.contact_vector for residue in category_contacts]

            print(f'Clustering {category} with {len(vector_array)} residue')

            if len(vector_array) > 1:
                cosine_distances = pdist(np.asarray(vector_array), metric='cosine')
                Z = linkage(cosine_distances, method='average')
                clusters = fcluster(Z, self.distance_cutoff, criterion='distance')
                self.clusters[category] = [(cluster, residue) for cluster, residue in zip(clusters, category_contacts)]

            # Handle singleton clusters by hand
            else:
                self.clusters[category] = [(0, category_contacts[0])]

    def generate_output_directories(self, base_dir, fragment_path, consolidate_clusters=True, output_source_pdb_info=True):
        """


        :param base_dir:
        :param fragment_path:
        :param consolidate_clusters:
        :param output_source_pdb_info:
        :return:
        """
        fragment_basepath, fragment = os.path.split(os.path.normpath(fragment_path))
        fragment_number = fragment.split('_')[1]

        # Generate list of fragment_PDB objects containing all relevant information on residue-ligand contacts in a cluster
        fragment_cluster_path = os.path.join(base_dir, 'Cluster_Results', fragment)
        os.makedirs(fragment_cluster_path, exist_ok=True)

        dict_list = []
        if output_source_pdb_info:
            source_pdb_dict = dict()

        for cluster_count, (category, cluster_annotated_residues) in enumerate(self.clusters.items(), start=1):

            source_pdb_dict[cluster_count] = dict()

            # todo: somehow separate annoted residue tuples into cluster lists
            # https://stackoverflow.com/questions/8092877/split-a-list-of-tuples-into-sub-lists-of-the-same-tuple-field
            cluster_annotated_residues_sorted = sorted(cluster_annotated_residues, key=operator.itemgetter(0))

            for cluster, cluster_residues in itertools.groupby(cluster_annotated_residues_sorted, operator.itemgetter(0)):

                # todo: simplify this jankiness
                cluster_residues = list(cluster_residues)

                # Generate report on cluster qualities
                dict_list.append(self.generate_report_row(cluster_residues, f'{cluster_count}-{cluster}'))

                # Add PDB sources to source_pdb_dict
                cluster_source_list = [res[1].pdb_info.split(".")[0].split('_')[0] for res in cluster_residues]
                cluster_source_set = list(set(cluster_source_list))
                source_pdb_dict[cluster_count][cluster] = {'set': cluster_source_set, 'list': cluster_source_list}

                # Output residue PDBs for each cluster
                # todo: add option to output clusters as PDB alongside ag.npz files
                if consolidate_clusters:

                    def tag_residue(cluster, res, index=2):
                        """"
                        Tag residue CA with source PDB string
                        :param cluster: cluster index
                        :param res: Fragment_PDB() container instance
                        """
                        residue_atomgroup = res.prody_residue.copy()
                        contact_source_filename = res.pdb_info.split(".")[0]
                        contact_source_pdbid = contact_source_filename.split('_')[0]
                        residue_source_string = f'F_{int(fragment_number)}-C_{int(cluster_count)}_{int(cluster)}-{res.prody_residue.getResname()}_{res.prody_residue.getResnum()}-{contact_source_filename}'

                        residue_atomgroup.setData('contact_source', [residue_source_string for atom in residue_atomgroup])
                        residue_atomgroup.setData('source_PDB', [contact_source_pdbid for atom in residue_atomgroup])
                        residue_atomgroup.setData('cluster', [f'{cluster_count}-{cluster}' for atom in residue_atomgroup])
                        residue_atomgroup.setData('fragment_id', [int(fragment_number) for atom in residue_atomgroup])

                        residue_atomgroup.setResnums([index] * len(residue_atomgroup))
                        residue_atomgroup.setChids(['A'] * len(residue_atomgroup))

                        return residue_atomgroup

                    # Start cluster_ensemble using first element... it seems you can't add residues to an empty AtomGroup
                    # todo: simplify this jankiness
                    cluster_ensemble = tag_residue(*list(cluster_residues)[0])

                    for index, (cluster_number, residue_container) in enumerate(cluster_residues[1:], start=3):
                        cluster_ensemble = cluster_ensemble + tag_residue(cluster_number, residue_container, index=index)

                    prody.saveAtoms(cluster_ensemble, filename=os.path.join(fragment_cluster_path, f'Cluster_{cluster_count}_{cluster}'))

                    # Write pdb
                    # todo: each residue should be in its own coordset (for visualization!!!)
                    prody.writePDB(os.path.join(fragment_cluster_path, f'Cluster_{cluster_count}_{cluster}'), cluster_ensemble)

                else:
                    for index, residue in enumerate(cluster_residues):
                        if residue[1].prody_residue.getResnames()[0] in ['ALA', 'CYS', 'SEC', 'ASP', 'GLU', 'PHE',
                                                                         'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET',
                                                                         'MSE', 'ASN', 'PRO', 'GLN', 'ARG', 'SER',
                                                                         'THR', 'VAL', 'TRP', 'TYR']:
                            prody.writePDB(os.path.join(fragment_cluster_path,
                                                        f'Cluster_{cluster_count}_{cluster}-{residue[1].prody_residue.getResnames()[0]}_{index}-{residue[1].pdb_info}'
                                                        ),
                                           residue[1].prody_residue)

        if output_source_pdb_info:
            with open(os.path.join(fragment_cluster_path, f'{fragment}-PDB_Sources.pickle'), 'wb') as pdb_sources:
                pickle.dump(source_pdb_dict, pdb_sources)

        # Export .csv
        self.df = pd.DataFrame(dict_list)
        self.df.to_csv(os.path.join(fragment_cluster_path,'{}_report.csv'.format(fragment)))

    def generate_report_row(self, residue_list, cluster_number):
        """
        Generates reports on each cluster
        I guess the first step would be to basically output the information contained in each of the representative vectors?
        :param residue_list: list of residue objects for a given cluster index
        :param cluster_number: cluster index, a bit redundant since I can get this from residue_list[0][0]
        :return: 
        """
        # Determine spread in min ligand-residue atom contact unit vectors
        # atom_atom_centroid_vector = np.mean([residue[1].contact_unit_vector for residue in residue_list], axis=0)
        # atom_atom_angle_deviations = [np.arccos(np.dot(atom_atom_centroid_vector, residue[1].contact_unit_vector)) for residue in residue_list]
        # atom_atom_mean = np.degrees(np.mean(atom_atom_angle_deviations))
        # atom_atom_SD = np.degrees(np.std(atom_atom_angle_deviations))

        # 20170901 - I've determined that the atom-centroid mean and SD aren't useful metrics...
        # Determine spread in unit vectors connecting min ligand contact atom and residue centroids
        # atom_centroid_vectors = [(residue[1].residue_center - residue[1].ligand_contact_atom.getCoords()[0]) /
        #                          np.linalg.norm(residue[1].residue_center - residue[1].ligand_contact_atom.getCoords()[0])
        #                          for residue in residue_list]
        # atom_centroid_centroid_vector = np.mean(atom_centroid_vectors, axis=0)
        #
        # atom_centroid_angle_deviations = [np.arccos(np.dot(atom_centroid_centroid_vector, atom_centroid_vector)) for atom_centroid_vector in atom_centroid_vectors]
        # atom_centroid_mean = np.degrees(np.mean(atom_centroid_angle_deviations))
        # atom_centroid_SD = np.degrees(np.std(atom_centroid_angle_deviations))

        return {'cluster_index': cluster_number,
                'cluster_members': len(residue_list),
                'cutoff': self.distance_cutoff,
                # 'atom-atom_mean': atom_atom_mean,
                # 'atom-atom_SD': atom_atom_SD
                }

    def automated_cluster_selection(self):
        """
        Automatically selects clusters based on number of members in each clusters. Clusters >1 SD above mean will be 
        selected for motif residue pool.
        :return: 
        """
        # Bruh. Paths.
        motif_yaml_path = os.path.join(self.processed_PDBs_dir.split('/')[0], 'Inputs', 'User_Inputs', 'Motif_Clusters.yml')
        current_fragment = os.path.basename(os.path.normpath(self.processed_PDBs_dir))

        # Select worthy clusters
        sorted_df = self.df.sort_values('cluster_members',ascending=False)

        cluster_size_mean = self.df.agg({'cluster_members': 'mean'})[0]
        cluster_size_SD = self.df.agg({'cluster_members': 'std'})[0]

        print('Cluster size mean:', cluster_size_mean)
        print('Cluster size SD:', cluster_size_SD)

        # Extended Mean - 0.5 SD
        selected_cluster_rows = sorted_df[sorted_df.cluster_members > (cluster_size_mean - 0.5 * cluster_size_SD)]

        # Mean + 1 SD
        # selected_cluster_rows = sorted_df[sorted_df.cluster_members > (cluster_size_mean + cluster_size_SD)]

        # Mean
        # selected_cluster_rows = sorted_df[sorted_df.cluster_members > cluster_size_mean]

        print('Selected clusters:')
        pprint(selected_cluster_rows)

        # with open(motif_yaml_path, 'a') as motif_yaml:
        #     motif_yaml.write('{}:\n'.format(current_fragment))
        #     for index, row in selected_cluster_rows.iterrows():
        #         motif_yaml.write('- {}\n'.format(int(row['cluster_index'])))

class fragment_PDB(object):
    """
    Class for holding all information related to a processed fragment PDB
    
    :param pdb_info:
    :param prody_pdb:
    :param prody_residue:
    :param pdbid:
    :param vector:
    """
    def __init__(self, prody_residue, pdb_info, prody_ligand):
        self.pdb_info = pdb_info
        self.prody_residue = prody_residue
        self.categorical_array = None
        self.contact_vector = None
        self.viable = self.process_residue_into_vector(prody_ligand, prody_residue)

    def process_residue_into_vector(self, prody_ligand, prody_residue, distance_cutoff_polar=3.5, distance_cutoff_greasy=4):
        """
        Converting each residue into a representative vector

        Metrics I care about for each residue:
        * Contact distance
        * Residue characteristics 
            * Amino acid identity or degenerate amino acid groups?
            * Amino acid chemical characteristics e.g. [Sigma Amino Acid Reference Chart]
              (http://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html)
        * Position of residue relative to fragment
            * Vector from fragment centroid to {residue centroid | closest residue atom }
        * Backbone or side chain
        
         {Angstroms} - X component, vector from fragment centroid to closest residue atom
         {Angstroms} - Y component, vector from fragment centroid to closest residue atom
         {Angstroms} - Z component, vector from fragment centroid to closest residue atom

         { 0 | 1 }   - Backbone contact OR Sidechain contact
         { 0 | 1 }   - Ligand Polar Contact OR Ligand Non-polar Contact
         
         { 0 | 1 }   - Side chain has hydrogen bond donor/acceptor (DEHKNQRSTY)
         
         { 0 | 1 }   - Hydrophobic, aliphatic (AILV)
         { 0 | 1 }   - Hydrophobic, aromatic (FWY)
         { 0 | 1 }   - Polar (NCQMST)
         { 0 | 1 }   - Charged, Acidic (DE)
         { 0 | 1 }   - Charged, Basic (HKR)
         { 0 | 1 }   - Glycine
         { 0 | 1 }   - Proline
         { 0 | 1 }   - Backbone carbonyl
         { 0 | 1 }   - Backbone amino
         { 0 | 1 }   - Backbone C/CA
        
        :return: 
        """

        min_contact_distance, row_index_low, column_index_low = minimum_contact_distance(prody_residue, prody_ligand, return_indices=True)
        polar_residues = ['ASP', 'GLU', 'HIS', 'LYS', 'ASN', 'GLN', 'ARG', 'SER', 'THR', 'TYR']

        if all([prody_residue.getResnames()[0] in polar_residues, min_contact_distance > distance_cutoff_polar]):
            return None

        elif min_contact_distance > distance_cutoff_greasy:
            return None

        else:
            residue_contact_atom = prody_residue.copy().select('index {}'.format(row_index_low))
            ligand_contact_atom = prody_ligand.copy().select('index {}'.format(column_index_low))

            # Save min contact residue and ligand atom indicies for evaluating cluster quality later
            # self.min_contact_distance = min_contact_distance

            # Residue Contact Type
            residue_contact_type = 0 if residue_contact_atom.getNames()[0] in ['C', 'CA', 'N', 'O'] else 1

            # Ligand Contact Type
            ligand_contact_type = 1 if ligand_contact_atom.getNames()[0][0] in ['C'] else 0

            # Vector from fragment centroid to closest residue atom
            contact_vector = (residue_contact_atom.getCoords() - prody.calcCenter(prody_ligand))[0]
            # contact_unit_vector = contact_vector / np.linalg.norm(contact_vector)

            # Side chain has hydrogen bond donor/acceptor (DEHKNQRSTY)
            h_bond_donor_acceptor = 1 if residue_contact_atom.getResnames()[0] in polar_residues else 0

            # Polar atom on residue is contacting ligand
            residue_polar_contact = 1 if residue_contact_atom.getNames()[0][0] in ['O', 'N'] else 0

            # Residue characteristics
            # todo: UPDATE so that only one of the below can hold value of 1 at any given time
            # {0 | 1} - Hydrophobic, aliphatic(AILVC)
            greasy_ali = 1 if residue_contact_atom.getResnames()[0] in ['ALA', 'ILE', 'LEU', 'VAL', 'CYS'] else 0
            # {0 | 1} - Hydrophobic, aromatic(FWY)
            greasy_aro = 1 if residue_contact_atom.getResnames()[0] in ['PHE', 'TYR', 'TRP'] else 0
            # {0 | 1} - Polar(NCQMST)
            polar = 1 if residue_contact_atom.getResnames()[0] in ['ASN', 'CYS', 'GLN', 'MET', 'SER', 'THR'] else 0
            # {0 | 1} - Charged, Acidic(DE)
            charged_acid = 1 if residue_contact_atom.getResnames()[0] in ['ASP', 'GLU'] else 0
            # {0 | 1} - Charged, Basic(HKR)
            charged_basic = 1 if residue_contact_atom.getResnames()[0] in ['HIS', 'LYS', 'ARG'] else 0
            # {0 | 1} - Glycine
            glycine = 1 if residue_contact_atom.getResnames()[0] in ['GLY'] else 0
            # {0 | 1} - Proline
            proline = 1 if residue_contact_atom.getResnames()[0] in ['PRO'] else 0
            # {0 | 1} - Backbone carbonyl
            bb_carbonyl = 1 if residue_contact_atom.getNames()[0] in ['O'] else 0
            # {0 | 1} - Backbone amino
            bb_amino = 1 if residue_contact_atom.getNames()[0] in ['N'] else 0
            # {0 | 1} - Backbone C / CA
            bb_c_ca = 1 if residue_contact_atom.getNames()[0] in ['C', 'CA'] else 0

            categorical_vector = [
                residue_contact_type,
                ligand_contact_type,
                h_bond_donor_acceptor,
                residue_polar_contact,
                greasy_ali,
                greasy_aro,
                polar,
                charged_acid,
                charged_basic,
                glycine,
                proline,
                bb_carbonyl,
                bb_amino,
                bb_c_ca
            ]

            self.categorical_array = np.asanyarray(categorical_vector)
            self.contact_vector = contact_vector

            return True
