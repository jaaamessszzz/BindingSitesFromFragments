#!/usr/bin/env python3

from ..clustering import Cluster
from ..utils import *

def clustering(args):
    """
    Perform clustering of aligned protein environments for all fragments

    Usage:
      bsff cluster <user_defined_dir>

    Arguments:
      <user_defined_dir>      Path to project root directory

    :return:
    """

    for fragment in directory_check(os.path.join(args['<user_defined_dir>'], 'Transformed_Aligned_PDBs')):

        fragment_name = os.path.basename(os.path.normpath(fragment))
        if os.path.exists(os.path.join(args['<user_defined_dir>'], 'Cluster_Results', fragment_name)):
            continue

        print(f'Processing {fragment_name}...')
        cluster = Cluster(fragment)

        if len(cluster.pdb_object_list) > 2:
            cluster.cluster_scipy()

        if cluster.clusters is not None:
            cluster.generate_output_directories(args['<user_defined_dir>'], fragment)