*****
Align
*****

From the ``Compounds/`` directory, enter: ::

    bsff align <Compound_ID>

Where <Compound_ID> is the name of the project directory for your target molecule.

.. autofunction:: BindingSitesFromFragments.commands.bsff_align.align

It is **highly** recommended that you use a local copy of the PDB for this step, otherwise BSFF with need to download
each PDB file as it requires. Instructions to download a local copy of the PDB can be found `here
<http://www.wwpdb.org/ftp/pdb-ftp-sites>`_. Use the command under "Download coordinate files in PDB Format". Use
the ``--use_local_pdb_database`` option to tell BSFF where it should look for local PDB files.

For each PDB, the ``align`` function will identify all ligands that contain a fragment as a substructure and align each
substructure and all residues within 12A of the fragment onto the reference fragments in the ``Fragment_Inputs/``
directory. Protein-ligand complexes are only used if they pass certain filters:

    * All atoms in the ligand are represented in the PDB
    * No ligand atoms have alternate location records
    * User-defined RMSD cutoff for poor alignments of fragment substructure onto reference

Your project should now look like this: ::

    +-- Compounds
        +-- <Compound_ID>
            +-- Inputs
            +-- Transformed_Aligned_PDBs
                +-- Fragment_1
                +-- Fragment_2
                +-- ...
            +-- Fragment_PDB_Matches
            +-- PDB_search_results.json
