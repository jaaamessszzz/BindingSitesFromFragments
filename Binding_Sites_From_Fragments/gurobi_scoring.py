#!/usr/bin/env python3

import sqlite3
import subprocess
import os
import sys
import pandas as pd
from gurobipy import *
import pprint
import itertools

class score_with_gurobi():
    """
    Class for determining best combinations of motif residues using Rosetta's feature reporter system
    """
    def __init__(self, user_defined_dir):
        self.user_defined_dir = user_defined_dir
        self.resources_dir = os.path.join(os.path.dirname(__file__), '..', 'Additional_Files')

    def complete_everything(self):
        pass

    def generate_feature_reporter_db(self):
        """
        Generates a SQLITE3 database with all necessary two-body residue scores
        :return: 
        """
        # todo: update with config file paths and options...
        run_feature_reporter = subprocess.Popen(['/Users/jameslucas/Rosetta/main/source/bin/rosetta_scripts.macosclangrelease', # UPDATE
                                                 '-parser:protocol',
                                                 os.path.join(self.resources_dir, 'RosettaScripts', 'Two_body_residue_feature_reporter.xml'),
                                                 '-out:nooutput',
                                                 '-parser:script_vars',
                                                 'target={}'.format(self.user_defined_dir),
                                                 '-l',
                                                 './{}/Motifs/Residue_Ligand_Interactions/Single_Poses/single_pose_list.txt'.format(self.user_defined_dir),
                                                 '-extra_res_fa',
                                                 './{}/Inputs/Rosetta_Inputs/{}.params'.format(self.user_defined_dir, self.user_defined_dir)
                                                 ])
        run_feature_reporter.wait()

    def consolidate_scores(self):
        """
        Consolidates two-body terms into single score for two-body interactions
        Dumps these values into a new table in the source SQLITE3 database generated by the feature reporter
        :return: 
        """
        connection = sqlite3.connect('./{}/two_body_terms.db'.format(self.user_defined_dir))
        cursor = connection.cursor()

        # Make new table for consolidated_2b_scores
        cursor.execute(
            """
            create table if not exists consolidated_2b_scores (
            batch_id INTEGER,
            struct_id INTEGER,
            resNum1 INTEGER,
            resNum2 INTEGER,
            score_total INTEGER,
            FOREIGN KEY (struct_id, resNum1) REFERENCES residues(struct_id, resNum) DEFERRABLE INITIALLY DEFERRED,
	        FOREIGN KEY (struct_id, resNum2) REFERENCES residues(struct_id, resNum) DEFERRABLE INITIALLY DEFERRED,
	        PRIMARY KEY (batch_id, struct_id, resNum1, resNum2)
	        )
            """)

        # Get all relevant information from residue_scores_2b
        table = pd.read_sql_query(
            """
            SELECT residue_scores_2b.batch_id, residue_scores_2b.struct_id, residue_scores_2b.resNum1,
            residue_scores_2b.resNum2, score_types.score_type_name, residue_scores_2b.score_value from residue_scores_2b
            left join score_types on residue_scores_2b.score_type_id == score_types.score_type_id where 
            (score_types.score_type_name = 'fa_atr' or 
            score_types.score_type_name = 'fa_rep' or 
            score_types.score_type_name = 'fa_elec' or
            score_types.score_type_name = 'hbond_sc' or
            score_types.score_type_name = 'hbond_bb_sc'
            )
            """, connection)

        # Commit every 10000
        execute_count = 0
        for value_tuple, residue_pair in table.groupby(['batch_id', 'struct_id', 'resNum1', 'resNum2']):

            ref2015_weights = {'fa_atr': 1,
                               'fa_elec': 1,
                               'hbond_sc': 1,
                               'hbond_bb_sc': 1,
                               'fa_rep': 0.55
                               }

            total_score = sum([ref2015_weights[row['score_type_name']] * row['score_value'] for index, row in residue_pair.iterrows()])

            cursor.execute(
                "INSERT OR IGNORE INTO consolidated_2b_scores VALUES (?,?,?,?,?)",
                (int(value_tuple[0]), int(value_tuple[1]), int(value_tuple[2]), int(value_tuple[3]), float(total_score),)
            )
            execute_count += 1

            if execute_count == 10000:
                connection.commit()
                print('Committed 10000...')
                execute_count = 0

        connection.commit()
        print('Completed!')

    def do_gurobi_things(self):
        connection = sqlite3.connect('./{}/two_body_terms.db'.format(self.user_defined_dir))
        cursor = connection.cursor()

        residue_table = pd.read_sql_query("SELECT * from residues", connection)
        score_table = pd.read_sql_query("SELECT * from consolidated_2b_scores", connection)

        for struct_id, table in residue_table.groupby(['struct_id']):
            print(struct_id)
            print(table)

            # Set up model
            residue_interactions = Model("residue_interactions")

            # Add MIP binary variables
            MIP_var_list = []
            for index, row in table.iterrows():
                MIP_var_list.append(residue_interactions.addVar(vtype=GRB.BINARY, name=str(row['resNum'])))

            # Set up dict with pairwise scores
            score_dict = {}
            relevant_scores = score_table.groupby(['struct_id']).get_group(struct_id)
            for index, row in relevant_scores.iterrows():
                score_dict[(row['resNum1'], row['resNum2'])] = row['score_total']

            # Set objective function
            # Get all possible reisude pairs
            residue_pairs = itertools.combinations(MIP_var_list, 2)
            pprint.pprint([res for res in residue_pairs])
            residue_interactions.setObjective(quicksum((MIP_var_list[int(key[0] - 1)] * MIP_var_list[int(key[1] - 1)] * value) for key, value in score_dict.items()), GRB.MINIMIZE)

            # Add constraints
            residue_interactions.addConstr(quicksum(MIP_var_list) == 5)
            residue_interactions.addConstr(MIP_var_list[0] == 1)
            for index, row in relevant_scores.iterrows():
                if row['score_total'] > 0:
                    residue_interactions.addConstr(MIP_var_list[int(row['resNum1'] - 1)] + MIP_var_list[int(row['resNum2'] - 1)] <= 1)

            # Set Parameters
            residue_interactions.Params.PoolSolutions = 1000

            # Optimzie
            residue_interactions.optimize()

            break

