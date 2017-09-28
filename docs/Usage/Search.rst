***********************************
Search for Ligand-Protein Complexes
***********************************

From the ``Compounds/`` directory, enter: ::

    bsff search <Compound_ID>

Where <Compound_ID> is the name of the project directory for you target molecule. This function will search for proteins
bound to small molecules specificed in the PubChem search results for each fragment. PDBs will be downloaded and stored
locally so that your project should now look like this: ::

    +-- Compounds
    +-- <Compound_ID>
        +-- Inputs
        +-- Fragment_PDB_Matches
            +-- Fragment_1
            |   +-- <Ligand_1>
            |   |   +-- <PDB_1>
            |   |   +-- <PDB_2>
            |   +-- <Ligand_1>
            |       +-- <PDB_1>
            |       +-- <PDB_2>
            +-- Fragment_2
                +-- <Ligand_1>
                |   +-- <PDB_1>
                |   +-- <PDB_2>
                +-- <Ligand_1>
                    +-- <PDB_1>
                    +-- <PDB_2>

For each unique substructure-containing molecule, only protein-ligand complexes that pass a 70% sequence identity cutoff are
downloaded and stored locally to prevent over-representation of closely-related residue-ligand contacts.