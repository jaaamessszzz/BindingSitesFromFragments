****************
Generate Designs
****************

Design using BindingSitesFromFragments takes place in two two steps. First, we will need to generate a new contact pool
for contacts that may be recapitulated in the context of a protein-ligand complex. Second, rotamers need to be added to
Rosetta's Packer and flagged so that we bias incorporation of these residues with the special_rot score term.

.. _existingcomplexfuzzball:

Generating a Contact Pool for an Existing Complex
-------------------------------------------------

From the ``Compounds/`` directory, enter: ::

    bsff assemble <Compound_ID> existing <existing_complex_path> <ligand_params> <ligand_ref> [options]

Where `<existing_complex_path>` is the path to an existing protein-ligand complex, `<ligand_params>` is the params file
generated by Rosetta's molfile_to_params.py for the *exact* ligand in `<existing_complex_path>`, and `<ligand_ref>` is the
PDB file generated by Rosetta's molefile_to_params.py for the *exact* ligand in `<existing_complex_path>`.

For protein-ligand complexes generated from RosettaMatch:
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Use the params file and PDB used to generate the match. These are specified in the match PDB output filename and can be
found in the `<Compound_ID>/Inputs/Rosetta_Inputs` directory.

For all other protein-ligand complexes:
"""""""""""""""""""""""""""""""""""""""
You will need to use Rosetta's molfile_to_params.py to generate a params file and reference PDB for the *exact* ligand
in the protein-ligand complex you want to design. Use a program like Avogadro to extract the ligand from the PDB file
specified by `<existing_complex_path>`, save it as a Sybyl mol2 file, and use it as in input to molfile_to_params.py. Both
`<ligand_params>` and <ligand_ref> will be generated at the same time by molfile_to_params.py.

By default, the assembly protocol for an existing complex will attempt to clean the PDB file specified by
`<existing_complex_path>` to remove any ligands that do not share the same chemical component identifier as `<Compound_ID>`
and remove any residues that are missing backbone atoms. It will also attempt to replace MSE-> MET and SEC-> CYS. This
should not modify PDBs generated by RosettaMatch, and can be skipped using the `--skip_clean` option.

A contact pool will be generated for each ligand in the protein-ligand complex with the chemical component identifier as
indicated by `<Compound_ID>`.

Contact pools for existing complexes are generated in their own subdirectory under `Compounds/<Compound_ID>/Design` to
consolidate all files required for designing a protein-ligand complex with complementary RotamerSets as straightforward
as possible. The directory structure is as follows: ::

    +-- Compounds
            +-- <Compound_ID>
                +-- Inputs
                +-- Transformed_Aligned_PDBs
                +-- Fragment_PDB_Matches
                +-- Cluster_Results
                +-- Fuzzballs
                +-- Design
                    +-- Complex_name
                        +-- Fuzzballs
                        +-- Conformers
                        +-- Existing_complex-clean.pdb
                        +-- motif_residue_attributes.csv
                +-- PDB_search_results.json

The `Conformers` directory contains copies of <ligand_params> and <ligand_ref>. All information for the complex contact
pool are in the `Fuzzball` directory and `motif_residue_attributes.csv`. The cleaned version of the input protein-ligand
complex passed in as `<existing_complex_path>` is copied to the Design directory with `--clean` appended.

Designing a Protein-Ligand Complex with Complementary RotamerSets
-----------------------------------------------------------------

The following command will perform design on a specified protein-ligand complex: ::

    bsff design <ligand_conformer_path> <match_path> <match_residue_map> <params_path> [options]

Conveniently, all of these files are consolidated into one place by the assembly procedure for an existing complex.
This command can be run from anywhere, but should be run where you want you designs to be written.

Generating Complementary Rotamers for Design
--------------------------------------------

The `design` command uses the following function to generate complementary rotamers from a specified contact pool and
protein-ligand complex:

.. autofunction:: BindingSitesFromFragments.design.generate_fuzzball_contact_rotamersets

In this protocol, complementary rotamers are simply added to the Packer RotamerSets. The essential bits of code are
as follows: ::

    # Create a new RotamerSets
    rotamer_sets = rosetta.core.pack.rotamer_set.RotamerSetsFactory.create_rotamer_sets(match_pose)
    rotamer_sets.set_task(design_packer_task)
    rotamer_sets.initialize_pose_for_rotsets_creation(match_pose)
    rotamer_sets.build_rotamers(match_pose, sfxn, packer_neighbor_graph)

    # Add complementary rotamers to RotamerSets
    if use_complementary_rotsets:
        for position in viable_rotamers:
            if design_packer_task.design_residue(position):
                print(f"Adding complementary rotamers for position {position}")
                position_rotamer_set = rotamer_sets.rotamer_set_for_residue(position)

                # Add fuzzball rotamers to the appropriate rotamer_set in rotamer_sets
                if int(position_rotamer_set.resid()) == position:
                    for residue_type in viable_rotamers[position]:
                        print(f'Adding {len(viable_rotamers[position][residue_type])} {residue_type} rotamers at position {position}.')
                        for fuzz_rotamer in viable_rotamers[position][residue_type]:
                            position_rotamer_set.add_rotamer_into_existing_group(fuzz_rotamer)

    # Copy pose for design
    design_pose = match_pose.clone()
    design_path = os.path.join(designdir, f'{match_name}-{i}.pdb')

    # Perform design
    sfxn.setup_for_packing_with_rotsets(design_pose, rotamer_sets)
    rotamer_sets.prepare_sets_for_packing(design_pose, sfxn)
    ig = rosetta.core.pack.interaction_graph.InteractionGraphFactory.create_and_initialize_annealing_graph(design_packer_task, rotamer_sets, design_pose, sfxn, packer_neighbor_graph)
    rosetta.core.pack.pack_rotamers_run(design_pose, design_packer_task, rotamer_sets, ig)
    ig.clean_up_after_packing(design_pose)
    sfxn(design_pose)

We generate a new RotamerSets using the RotamerSets factory and then add rotamers to each position's RotamerSet based on
what is returned by `generate_fuzzball_contact_rotamersets()`. This augmented RotamerSets is then passed to the Packer
`pack_rotamers_run` method.
