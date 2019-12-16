#!/usr/bin/env python3

"""
Benchmark sequence recovery of an existing binding site with/without a complementary rotamerset at varying special_rot
weights.

"""

from multiprocessing import Pool
import os
import pickle
import random
import subprocess
import tarfile
import tempfile
from pprint import pprint

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import prody
from prody.sequence.sequence import Sequence
from scipy.spatial.distance import jensenshannon
import seaborn as sns

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align.Applications import ClustalOmegaCommandline

import pyrosetta
from pyrosetta import rosetta

from ..design import create_task_factory
from ..utils import *


def binding_interface_recovery_bm(args):
    """
    Generate HMMER for designs at varying special_rot weights and compare native to each profile

    For each design condition:
        Generate an MSA
        Build a profile with HMMbuild

    Usage:
      bsff benchmark binding_interface_recovery align <benchmark_dir> [options]
      bsff benchmark binding_interface_recovery hmmer <benchmark_dir> [options]
      bsff benchmark binding_interface_recovery jsd <benchmark_dir> [options]

    Arguments:
      <benchmark_dir>           Root Directory
      hmmer                     Generate HMMER Profiles

    Options:
      --rawfasta_dir=<fasta_dir>        Generate alignments using existing FASTAs in <fasta_dir>
      --alignedfasta_dir=<fasta_dir>    Use existing alignments in <fasta_dir>
    """
    # Directory with True Sequence FASTAs
    true_fastas = 'TrueSequences'
    benchmark_dir = args['<benchmark_dir>']

    if args['align']:
        mega_pickle, true_fastas = consolidate_design_sequences(benchmark_dir)
        fasta_path_dict = fasta_dump(mega_pickle)
        alignments_dir = generate_alignments(fasta_path_dict)

    if args['hmmer']:
        if args['--rawfasta_dir']:
            fasta_dir = args['--rawfasta_dir']
            fasta_path_dict = dict()
            for complex in os.listdir(fasta_dir):
                fasta_path_dict[complex] = list()
                for fasta in os.listdir(os.path.join(fasta_dir, complex)):
                    if fasta.endswith('.fasta'):
                        fasta_path_dict[complex].append(os.path.join(fasta_dir, complex, fasta))

            alignments_dir = generate_alignments(fasta_path_dict)
            profile_db_dir = generate_hmmer_profiles(alignments_dir)
            run_hmmscan(true_fastas, profile_db_dir)

        elif args['--alignedfasta_dir']:
            alignments_dir = args['--alignedfasta_dir']
            profile_db_dir = generate_hmmer_profiles(alignments_dir)
            run_hmmscan(true_fastas, profile_db_dir)

        else:
            mega_pickle, true_fastas = consolidate_design_sequences(benchmark_dir)
            fasta_path_dict = fasta_dump(mega_pickle)
            alignments_dir = generate_alignments(fasta_path_dict)
            profile_db_dir = generate_hmmer_profiles(alignments_dir)
            run_hmmscan(true_fastas, profile_db_dir)

    if args['jsd']:
        if args['--alignedfasta_dir']:
            alignments_dir = args['--alignedfasta_dir']
        else:
            mega_pickle, true_fastas = consolidate_design_sequences(benchmark_dir)
            fasta_path_dict = fasta_dump(mega_pickle)
            alignments_dir = generate_alignments(fasta_path_dict)

        jensen_shannon_divergence(benchmark_dir, alignments_dir)


def consolidate_design_sequences(root_dir, dump_dir=os.getcwd(), dump_true_sequences=True):

    # Walk through root directory and find all tar.gz or pickles
    # Multiprocess.Pool() on all tar.gz

    # Load TrueSequences.pickle if it already exists, otherwise load during walk
    # Convert add_sequences_to_dict() to take tarball path and add all pdbs within

    def add_sequences_to_dict(prody_hv, current_pdb):
        if prody_hv.numChains() > 1:
            # Consolidate symmetry
            consolidate = True
            design_chains = [a for a in prody_hv.iterChains()]
            chain1 = design_chains[0].getSequence()
            chain2 = design_chains[1].getSequence()
            true_sequence = true_sequences[current_pdb]

            # Workaround since I forgot to make sure dimers were completely symmertric...
            # There are zero occupancy and missing residues in dimers either at N/C term or loops distal to binding site
            if file.startswith('5J5G'):
                chain2 = chain2[:-5]
            elif file.startswith('5J5I'):
                chain1 = chain1[:-1]
                true_sequence = true_sequence[:-1]
            elif file.startswith('6m9b'):
                chain1 = chain1[1:]
                chain2 = chain2[:-1]
                true_sequence = true_sequence[1:]
            elif file.startswith('2IZL'):
                chain1 = chain1[:-1]
                true_sequence = true_sequence[:-1]
            elif file.startswith('4QAC'):
                chain1 = chain1[6:]
                chain2 = chain2[:157] + chain2[159:]
                true_sequence = true_sequence[6:]
            elif file.startswith('3CZ1'):
                chain1 = chain1[2:]
                true_sequence = true_sequence[2:]
            elif file.startswith('3oki'):
                consolidate = False
                design_sequence = [a for a in prody_hv.iterChains()][0].getSequence()
            elif file.startswith('2QRY'):
                chain2 = chain2[:-1]

            if not file.startswith('3oki'):
                assert len(chain1) == len(chain2)

            if consolidate:
                consolidated_sequence = list()
                for index, true_res in enumerate(true_sequence):
                    # print(index, true_res, chain1[index], chain2[index])
                    if chain1[index] == chain2[index]:
                        consolidated_sequence.append(chain1[index])
                    else:
                        mutation = {chain1[index], chain2[index]} - set(true_res)
                        consolidated_sequence.append(list(mutation)[0])
                design_sequence = ''.join(consolidated_sequence)

        else:
            # design_sequence = ''.join([THREE_TO_ONE[res.getResname()] for res in design_prody_hv.iterResidues()])
            design_sequence = [a for a in prody_hv.iterChains()][0].getSequence()

        # Add sequence to design_sequences_per_weight
        design_key = 'Control' if comprotset_off else special_rot_weight
        if design_sequences_per_weight[current_pdb].get(design_key):
            design_sequences_per_weight[current_pdb][design_key].append(design_sequence)
        else:
            design_sequences_per_weight[current_pdb][design_key] = [design_sequence]

    # True sequences should be found first before designs based on how os.walk traverses directories
    # Single chain
    true_sequences = dict()

    # Design sequences / weight / complex
    design_sequences_per_weight = dict()

    # Consolidate sequences per weight
    # This mess assumes you have all the RAM in the world to burn
    tarball_list = list()

    for root, dirs, files in os.walk(root_dir):

        path_head, dir_basename = os.path.split(root)

        for file in files:
            if file.endswith('-clean.pdb'):
                pdb_name = file.split('-')[0]
                true_prody = prody.parsePDB(os.path.join(root, file))
                true_prody_hv = true_prody.select('protein').getHierView()
                true_sequence = [a for a in true_prody_hv.iterChains()][0].getSequence()
                print(pdb_name, true_sequence)
                true_sequences[pdb_name] = true_sequence
                if design_sequences_per_weight.get(pdb_name) is None:
                    design_sequences_per_weight[pdb_name] = dict()

        if dir_basename.endswith('-clean'):
            dir_basename_pdb = dir_basename.split('-')[0]
            if design_sequences_per_weight.get(dir_basename_pdb) is None:
                design_sequences_per_weight[dir_basename_pdb] = dict()

        if dir_basename.startswith('Designs-taskid'):

            # Get weight
            comprotset_off = False

            if 'NoCompRotSet' in dir_basename:
                special_rot_weight = 0
                comprotset_off = True
            else:
                special_rot_weight = float(dir_basename.split('_')[-1])
                print(f'Current weight:\t{special_rot_weight}')

            if f'{dir_basename}.tar.gz' in files:
                with tarfile.open(os.path.join(root, f'{dir_basename}.tar.gz'), 'r|gz') as tar:
                    with tempfile.TemporaryDirectory() as tar_temp_dir:
                        tar.extractall(path=tar_temp_dir)

                        for root, dirs, files in os.walk(tar_temp_dir):
                            for file in files:
                                if file.endswith('.pdb'):
                                    complex_name = file.split('-')[0]
                                    design_prody = prody.parsePDB(os.path.join(root, file))
                                    design_prody_hv = design_prody.select('protein').getHierView()
                                    add_sequences_to_dict(design_prody_hv, complex_name)

            else:
                for file in files:

                    # DEBUGGING
                    # if file[:4] in ['5J5G', '5edb', '5HZ6', '5HZ8', '5J5G', '5J5I', '5ura', '1TOU', '6m9b', '1sri', '1lke',
                    #                 '1n0s', '2XN3', '2IZL', '4QAC', '4AFH', '3CZ1', '5t52', '3oki', '4AFG', '4B5D', '2QRY']:
                    #     continue

                    if file.endswith('.pdb'):
                        complex_name = file.split('-')[0]
                        design_prody = prody.parsePDB(os.path.join(root, file))
                        design_prody_hv = design_prody.select('protein').getHierView()
                        add_sequences_to_dict(design_prody_hv, complex_name)


    # Dump design_sequences_per_weight as pickle
    with open(os.path.join(dump_dir, 'ConsolidatedSpecialRotDesignSequences.pickle'), 'wb') as mega_pickle:
        pickle.dump(design_sequences_per_weight, mega_pickle)

    # Dump true_sequences as pickle
    with open(os.path.join(dump_dir, 'TrueSequences.pickle'), 'wb') as true_sequence_pickle:
        pickle.dump(true_sequences, true_sequence_pickle)

    # Dump true_sequences as individual FASTA files
    true_fastas = os.path.join(dump_dir, 'TrueSequences')
    if dump_true_sequences:
        os.makedirs(true_fastas, exist_ok=True)
        for complex in true_sequences:
            seqrecord = SeqRecord(Seq(true_sequences[complex], IUPAC.protein), id=f'{complex}', name=f'{complex}')
            fasta_path = os.path.join(true_fastas, f'{complex}.fasta')
            SeqIO.write(seqrecord, fasta_path, 'fasta')

    return design_sequences_per_weight, true_fastas


def fasta_dump(design_sequences_per_weight, dump_dir='DesignFASTA-Raw'):
    """
    Dump consolidated sequences into FASTA files / complex / weight
    :param design_sequences_per_weight: dictionary containing design sequences per complex per weight
    :param dump_dir:
    :return:
    """
    os.makedirs(dump_dir, exist_ok=True)
    fasta_path_dict = dict()

    # Generate FASTA files / weight / complex
    def create_seqrecord(sequence, complex, weight, index):
        return SeqRecord(Seq(sequence, IUPAC.protein), id=f'{complex}_{weight}_{index}', name=f'{complex}')

    for complex in design_sequences_per_weight:
        complex_subdir = os.path.join(dump_dir, complex)
        os.makedirs(complex_subdir, exist_ok=True)
        fasta_path_dict[complex] = list()

        for weight in design_sequences_per_weight[complex]:
            # Convert sequences into Seq Records
            seqrecords = [create_seqrecord(sequence, complex, weight, index) for index, sequence in enumerate(design_sequences_per_weight[complex][weight])]

            # Output FASTA
            fasta_path = os.path.join(complex_subdir, f'{complex}_{weight}.fasta')
            fasta_path_dict[complex].append(fasta_path)
            SeqIO.write(seqrecords, fasta_path, 'fasta')

    return fasta_path_dict


def generate_alignments(fasta_path_dict, dump_dir='DesignFASTA-Aligned'):
    """
    Technically not necessary since we're comparing designs on the same backbone, but included for completeness
    :param raw_fasta_paths: list of paths to unaligned FASTA files
    :return:
    """
    os.makedirs(dump_dir, exist_ok=True)

    # Generate alignment with ClustalOmega
    for complex in fasta_path_dict:
        complex_subdir = os.path.join(dump_dir, complex)
        os.makedirs(complex_subdir, exist_ok=True)

        for fasta_path in fasta_path_dict[complex]:
            path_base, fasta_file = os.path.split(fasta_path)
            out_file = os.path.join(complex_subdir, fasta_file)
            clustalomega_cline = ClustalOmegaCommandline(infile=fasta_path, outfile=out_file, verbose=True, auto=True, force=True)
            stdout, stderr = clustalomega_cline()

    return dump_dir

# --- HMMER Profiles --- #

def generate_hmmer_profiles(alignments_dir, dump_profiles_dir='DesignHMMERProfiles', dump_db_dir='DesignHMMERDB'):
    """
    hmmbuild [options] hmmfile msafile
      -n    <new_profile_path>
      -o    <summary_output_path>
      -O    Resave alignment in annotated Stockholm format

    hmmpress [options] hmmfile
      -f    Force output, overwrites existing profiles

    :param alignments_dir: Directory containing alignments
    :return:
    """

    os.makedirs(dump_profiles_dir, exist_ok=True)

    # Generate Profiles with HMMER
    for complex in os.listdir(alignments_dir):
        complex_dir = os.path.join(dump_profiles_dir, complex)
        os.makedirs(complex_dir, exist_ok=True)

        for alignment in os.listdir(os.path.join(alignments_dir, complex)):
            alignment_path = os.path.join(alignments_dir, complex, alignment)
            alignment_filename, _ = os.path.splitext(alignment)

            hmm_profile_filename = f'{alignment_filename}.hmm'
            hmm_profile_stdout = os.path.join(complex_dir, f'{alignment_filename}.out')
            hmm_profile_path = os.path.join(complex_dir, hmm_profile_filename)

            cmd = ['hmmbuild', '-n', f'{hmm_profile_filename}', '-o', f'{hmm_profile_stdout}', f'{hmm_profile_path}', f'{alignment_path}']
            subprocess.run(cmd)
            # f'hmmbuild -n {hmm_profile_filename} -o {hmm_profile_stdout} {hmm_profile_path} {alignment_path}'

    # Concat all HMMER Profiles
    os.makedirs(dump_db_dir, exist_ok=True)

    for complex in os.listdir(dump_profiles_dir):
        complex_profile_dir = os.path.join(dump_profiles_dir, complex)

        with open(os.path.join(dump_db_dir, f'{complex}-ProfileDB.hmm'), 'w') as hmmer_db:
            for profile in os.listdir(complex_profile_dir):
                if profile.endswith('.hmm'):
                    with open(os.path.join(complex_profile_dir, profile), 'r') as single_profile:
                        for line in single_profile:
                            hmmer_db.write(line)

    # Prepare profile databases for hmmscan
    for profile_db in os.listdir(dump_db_dir):
        if profile_db.endswith('-ProfileDB.hmm'):
            subprocess.run([f'hmmpress', '-f', f'{os.path.join(dump_db_dir, profile_db)}'])

    return dump_db_dir


def run_hmmscan(true_sequence_dir, profile_db_dir):
    """
    Compares true sequence to HMMER profiles for special_rot weights

    :param true_sequence_dir: Directory containing FASTA files of true sequences
    :return:
    """
    # Compare sequence to HMMER Profiles
    for true_sequence in os.listdir(true_sequence_dir):
        complex, _ = os.path.splitext(true_sequence)
        fasta_path = os.path.join(true_sequence_dir, true_sequence)
        complex_profile_db = os.path.join(profile_db_dir, f'{complex}-ProfileDB.hmm')
        subprocess.run(['hmmscan', '--tblout', f'{complex}-HMMscanResults.tsv', f'{complex_profile_db}', f'{fasta_path}'])


# --- Noah's CoupledMoves Native Sequence Recovery --- #

def jensen_shannon_divergence(benchmark_dir, aligned_fasta_dir, recovery_cutoff=0.5):
    """
    Calculate profile similarity between true sequence and designs using the phi-JS divergence as described in:

        Within the Twilight Zone: A Sensitive Profile-Profile Comparison Tool Based on Information Theory
        Golan Yona* and Michael Levitt

    where phi=0.5.

    :param aligned_fasta_dir: path to DesignFASTA-Aligned directory
    :param true_sequence_dir: path to TrueSequences directory
    :param phi: value of phi to use in phi-JS divergence
    :return:
    """

    # I messed up and forgot to map complex:ligand somewhere easy... so have this instead
    complex_ligand_map = dict()
    for root, dirs, files in os.walk(benchmark_dir):
        for file in files:
            if file.endswith('-clean.pdb'):
                complex_ligand_map[file.split('-')[0]] = root.split('/')[1]

    # Make directory tree to dump statistics
    root_dir = 'JensenShannonDivergence'
    os.makedirs(root_dir, exist_ok=True)

    # for loops all the way down...
    for complex in os.listdir(aligned_fasta_dir):

        complex_dir = os.path.join(aligned_fasta_dir, complex)
        if os.path.isdir(complex_dir):

            # Get designable positions from input complex pose
            benchmark_base = os.path.join(benchmark_dir, complex_ligand_map[complex], 'Design', f'{complex}-clean')
            params_path = os.path.join(benchmark_base, 'Conformers', f'{complex_ligand_map[complex]}.params')
            complex_pdb_path = os.path.join(benchmark_base, f'{complex}-clean.pdb')

            # Load params and pose as in 07.03-Ligand-Docking-PyRosetta.ipynb
            pyrosetta.init()
            complex_pose = pyrosetta.rosetta.core.pose.Pose()
            ligand_params = pyrosetta.Vector1([params_path])
            residue_type_set = complex_pose.conformation().modifiable_residue_type_set_for_conf()
            residue_type_set.read_files_for_base_residue_types(ligand_params)
            complex_pose.conformation().reset_residue_type_set_for_conf(residue_type_set)
            rosetta.core.import_pose.pose_from_file(complex_pose, complex_pdb_path)

            # Normal way
            # pyrosetta.init('-extra_res_fa', params_path)
            # complex_pose = rosetta.core.import_pose.pose_from_file(complex_pdb_path)

            ligand_name = complex_pose.residue(complex_pose.size()).name3()

            # Get task_factory and get designable positions for complex_pose
            task_factory = create_task_factory(complex_pose, complex_pdb_path)
            design_packer_task = task_factory.create_task_and_apply_taskoperations(complex_pose)
            designable_positions = [index for index, res in enumerate(design_packer_task.designing_residues(), start=1) if res is True]

            # todo: write out designable positions to a pickle somewhere...

            # Fix designable positions for ternary+ complexes, PART 1
            if complex_pose.num_chains() > 2 and complex_pose != '3oki':
                ligand = complex_pose.split_by_chain(3)

                chain1_pose = complex_pose.split_by_chain(1)
                chain1_pose.append_pose_by_jump(ligand, chain1_pose.size())
                chain1_task_factory = create_task_factory(chain1_pose, complex_pdb_path)
                chain1_packer_task = chain1_task_factory.create_task_and_apply_taskoperations(chain1_pose)

                chain2_pose = complex_pose.split_by_chain(2)
                chain2_pose.append_pose_by_jump(ligand, chain2_pose.size())
                chain2_task_factory = create_task_factory(chain2_pose, complex_pdb_path)
                chain2_packer_task = chain2_task_factory.create_task_and_apply_taskoperations(chain2_pose)

            # Get all interface positions
            complex_prody = prody.parsePDB(complex_pdb_path)
            contact_residue_selection = complex_prody.select(f'protein within 4 of resname {ligand_name}')

            # Dump data into df
            list_of_dicts = list()

            # Iterate through alignments and calculate JSD
            for fasta in os.listdir(complex_dir):

                fasta_name, fasta_ext = os.path.splitext(fasta)
                print(fasta_name, fasta_ext)
                if fasta_ext == '.fasta':

                    fasta_path = os.path.join(complex_dir, fasta)
                    fasta_split = fasta_name.split('_')
                    fasta_weight = 'Control' if fasta_split[1] == 'Control' else float(fasta_split[1])

                    # Align pose sequence to MSA to get "correct" numbering
                    true_sequence_msa = prody.parseMSA(fasta_path)
                    alignment, seq_indices, msa_indices = prody.alignSequenceToMSA(complex_prody, true_sequence_msa, f'{complex}_{fasta_weight}_0')
                    seq_to_msa = {seq: msa for seq, msa in zip(seq_indices, msa_indices)}  # seq is 0 indexed!!! position/pose is 1 indexed!!!

                    # todo: need to map msa to chain A in oligomers for correct designable position indices (messes up ligand contact annotation)
                    msa_to_seq = {v: k for k, v in seq_to_msa.items()}

                    # if complex_pose.num_chains() > 2 and complex_pose != '3oki':
                    #     msa_to_seq = {v: k for k, v in seq_to_msa.items() if k <= chain1_pose.size()}
                    # else:
                    #     msa_to_seq = {v: k for k, v in seq_to_msa.items()}

                    # Fix designable positions for ternary+ complexes, PART 2
                    if complex_pose.num_chains() > 2 and complex_pose != '3oki':
                        chain1_seq = complex_prody.select(f'chain {chain1_pose.pdb_info().chain(1)}').copy()
                        alignment, seq_indices, msa_indices = prody.alignSequenceToMSA(chain1_seq, true_sequence_msa, f'{complex}_{fasta_weight}_0')
                        chain1_seq_to_msa = {seq: msa for seq, msa in zip(seq_indices, msa_indices)}

                        chain2_seq = complex_prody.select(f'chain {chain2_pose.pdb_info().chain(1)}').copy()
                        alignment, seq_indices, msa_indices = prody.alignSequenceToMSA(chain2_seq, true_sequence_msa, f'{complex}_{fasta_weight}_0')
                        chain2_seq_to_msa = {seq: msa for seq, msa in zip(seq_indices, msa_indices)}

                        chain1_designable = [chain1_seq_to_msa[index] for index, res in enumerate(chain1_packer_task.designing_residues()) if res is True]
                        chain2_designable = [chain2_seq_to_msa[index] for index, res in enumerate(chain2_packer_task.designing_residues()) if res is True]

                        designable_positions = sorted([msa_to_seq[pos] + 1 for pos in (set(chain1_designable) | set(chain2_designable))])

                        print(chain1_seq_to_msa)
                        print(chain2_seq_to_msa)

                    print('DESIGNABLE:', designable_positions)

                    # Calculate JSD and "native sequence recovery" for design positions
                    design_position_info_raw = {position: list() for position in designable_positions}
                    for sequence in range(true_sequence_msa.numSequences()):
                        for position in designable_positions:
                            design_position_info_raw[position].append(str(true_sequence_msa[sequence, seq_to_msa[position - 1]]))

                    for position in design_position_info_raw:
                        true_resid = complex_pose.residue(position).name1()
                        unique_residues = list(set(design_position_info_raw[position] + [true_resid]))
                        seqeunce_count = len(design_position_info_raw[position])
                        resn_counts = {resn: design_position_info_raw[position].count(resn) for resn in unique_residues}
                        recovered = resn_counts[true_resid] / seqeunce_count >= recovery_cutoff

                        # phi-JSD where phi=0.5
                        p = [1 if true_resid == res else 0 for res in unique_residues]
                        q = [resn_counts[res] / seqeunce_count if res in resn_counts.keys() else 0 for res in unique_residues]
                        jsd = jensenshannon(p, q, base=2) ** 2

                        # Contact residues
                        contact_residues = set(designable_positions) & {atom.getResnum() for atom in contact_residue_selection}

                        df_row = {'complex': fasta_split[0],
                                  'weight': fasta_weight,
                                  'position': position,
                                  'truth': true_resid,
                                  'recovered': recovered,
                                  'jsd': jsd,
                                  'similarity': 1 - jsd,
                                  'contact': position in contact_residues,
                                  }

                        list_of_dicts.append(df_row)

            df = pd.DataFrame(list_of_dicts)
            df.to_csv(os.path.join(root_dir, f'{complex}_stats.csv'), index=False)

# --- Count things in dataframes --- #

def fetch_top_sequence_identity(alignment_dir, truesequence_dir):
    """
    Gets seqeunce with best sequence identity to truth for each alignment FASTA

    :param alignment_dir:
    :param truesequence_dir:
    :return:
    """

    # Load true sequence
    # Load alignment
    # Add true sequence to alignment
    # Calculate sequence identity matrix
    # Get max value from sequence identity matrix

    for complex in os.listdir(alignment_dir):

        complex_dir = os.path.join(alignment_dir, complex)
        true_sequence = prody.parseMSA(os.path.join(truesequence_dir, f'{complex}.fasta'))

        for alignment_fasta in os.listdir(complex_dir):

            alignment_fasta_path = os.path.join(complex_dir, alignment_fasta)
            weight_msa = prody.parseMSA(alignment_fasta_path)



# --- Figures --- #

def binding_interface_figures(args):
    """
    Generate HMMER profile similarity figures

    Usage:
      bsff benchmark binding_interface_figures <statistics_dir>

    Arguments:
      <statistics_dir>          CSV file generated by binding_interface_recovery_bm()

    Options:
      --derp                    asdf
    """
    statistics_dir = args['<statistics_dir>']

    table_headers = ['target name',
                     'accession',
                     'query name',
                     'accession',
                     'E - value',
                     'score',
                     'bias',
                     'E - value',
                     'score',
                     'bias',
                     'exp',
                     'reg',
                     'clu',
                     'ov',
                     'env',
                     'dom',
                     'rep',
                     'inc',
                     'description of target']

    consolidated_df = pd.DataFrame()

    for hmm_results in os.listdir(statistics_dir):
        print(os.path.join(statistics_dir, hmm_results))
        if hmm_results.endswith('.tsv'):
            df = pd.read_csv(os.path.join(statistics_dir, hmm_results), delim_whitespace=True, comment='#', header=None, names=table_headers)
            print(df)


# Look up Noah's CoupledMoves Native Sequence Recovery
# Look up Colin's PDZ doamin paper for metrics
def jsd_figures(args):
    """
    Generate figures for profile similarity and native sequence recovery

    Usage:
      bsff benchmark jsd_figures <statistics_dir>

    Arguments:
      <statistics_dir>          Directory generated by jensen_shannon_divergence()
    """
    stats_dir = args['<statistics_dir>']

    # Consolidate csvs into one dataframe
    all_csvs = [os.path.join(stats_dir, csv) for csv in os.listdir(stats_dir) if csv.endswith('.csv')]
    df_list = [pd.read_csv(csv) for csv in all_csvs]
    df_all = pd.concat(df_list)

    # --- Profile Similarity (Box plots) --- #
    sns.set(style="whitegrid", palette="hls")

    def make_boxplot(df, output_name=None):
        complex_name = df['complex'][0]

        plt.figure()
        boxplot = sns.boxplot(x=df['weight'], y=df['similarity'],
                              order=['Control', '0.0', '-0.5', '-1.0', '-1.5', '-2.0', '-2.5', '-3.0', '-3.5', '-4.0',
                                     '-4.5', '-5.0', '-6.0', '-8.0', '-10.0'])
        sns.despine(left=True)

        fig = boxplot.get_figure()
        plt.xlabel('Weights')
        plt.ylabel('Profile Similarity')
        fig.savefig(f'ProfileSimilarity-{output_name if output_name else complex_name}.png')
        plt.close()

    for df in df_list:
        make_boxplot(df)

    # Aggregate information from all complexes
    make_boxplot(df_all, output_name='Aggregate')
    df_weight_pivot = pd.pivot_table(df_all, values='similarity', index=['weight'],
                                     aggfunc={'similarity': [np.median, np.mean]})
    df_weight_pivot.to_csv('ProfileSimilarity-AggregateStats.csv')

    # --- Sequence Recovery (Bar Chart) --- #

    def make_barchart(df, output_name=None):
        complex_name = df['complex'][0]

        list_of_dicts = list()

        for weight, df_sub in df.groupby('weight'):
            list_of_dicts.append({'weight': weight, 'count': df_sub['recovered'].value_counts()[True] / len(df_sub['recovered'])})

        count_df = pd.DataFrame(list_of_dicts)

        plt.figure()
        barplot = sns.barplot(x=count_df['weight'], y=count_df['count'],
                              order=['Control', '0.0', '-0.5', '-1.0', '-1.5', '-2.0', '-2.5', '-3.0', '-3.5', '-4.0',
                                     '-4.5', '-5.0', '-6.0', '-8.0', '-10.0'],
                              linewidth=1.5, edgecolor=".15",)
        sns.despine(left=True)

        fig = barplot.get_figure()
        plt.ylim(0, 1)
        plt.xlabel('Weights')
        plt.ylabel('Fraction Recovered')
        fig.savefig(f'SeqeuenceRecovery-{output_name if output_name else complex_name}.png')
        plt.close()

    for df in df_list:
        make_barchart(df)

    # Aggregate information from all complexes
    make_barchart(df_all, output_name='Aggregate')
