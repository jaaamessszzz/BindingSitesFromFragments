#!/usr/bin/env python3
import os

def new_project(args):
    """
    Generate a new project.

    This command will generate a new working directory for your target ligand named <chemical_id>. You can name your project
    anything, but the first three characters of your project directory need to be the chemical component identifier of your
    target ligand. A chemical component identifier is typically a three-character string that the PDB uses to uniquely
    identify ligands (e.g. ATP). For instance, acceptable project names can be simply 'ATP' or 'ATP-binders' if a descriptor
    is desired. If your ligand does not exist in the PDB, any three-character string will suffice. BSFF will use the
    chemical component identifier defined in the project name to keep track of your ligand throughout the protocol. This is
    to avoid the necessity of a config file...

    Within your new project directory, an `Inputs` directory is generated with the following subdirectories:

    `Rosetta_Inputs`:   All information related to the target molecule live here. To start, find a Sybyl Mol2 representation
                        of your target ligand (or a conformer library, if desired). We will use Rosetta’s
                        molfile_to_params.py to standardize atom names and prepare .params files so that Rosetta will know
                        how to handle your molecule. The `molfile_split_into_singles.py` found under `Additional_Files` will
                        take an input .mol2 file and populate `Rosetta_Inputs` with all the required .params and .pdb files.

    `Fragment_Inputs`:  All information related to the fragments that make up your target ligand live here. All fragments
                        need to be generated from a single .pdb file produced by Rosetta's `molfile_to_params.py` (this was
                        handled by `molfile_split_into_singles.py`.
                        'Fragment_inputs.csv` has been automatically generated for you to save the SMILES queries for each
                        fragment search performed through PubChem. The PubChem queries *will eventually be* used to properly
                        map your defined fragments onto the target ligand.

    `User_Inputs`:      `Exclude_Ligands.txt` can be used to specify ligands that you want to ignore.

    Usage: bsff new <chemical_id>

    """
    print('Starting a new project for {}'.format(args['<chemical_id>']))
    project_root_dir = os.path.join(os.getcwd(), args['<chemical_id>'])
    os.makedirs(project_root_dir)
    os.makedirs(os.path.join(project_root_dir, 'Inputs'))
    os.makedirs(os.path.join(project_root_dir, 'Inputs', 'Rosetta_Inputs'))
    os.makedirs(os.path.join(project_root_dir, 'Inputs', 'Fragment_Inputs'))
    os.makedirs(os.path.join(project_root_dir, 'Inputs', 'User_Inputs'))
    with open(os.path.join(project_root_dir, 'Inputs', 'Fragment_Inputs', 'Fragment_inputs.csv'), 'w') as fragment_inputs:
        fragment_inputs.write('Fragment,SMILES_fragment')
