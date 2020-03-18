***********************************
Search for Ligand-Protein Complexes
***********************************

From the ``Compounds/`` directory, enter: ::

    bsff search <Compound_ID>

Where <Compound_ID> is the name of the project directory for you target molecule. This function will search for proteins
bound to small molecules specificed in the PubChem search results for each fragment. A JSON file called
`PDB_search_results.json` keeps track of ligands and associated complexes in the PDB that possess each defined
fragment as a substructure. In addition, a directory called `Fragment_PDB_Matches` will be created that stores the
PubChem Compound IDs for fragment-containing compounds in CSV format.

Your project should now look like this: ::

    +-- Compounds
    +-- <Compound_ID>
        +-- Inputs
            +-- Fragment_Inputs
            |   +-- Fragment_Inputs.csv
            |   +-- * All fragment .pdb files *
            |   +-- * All fragment .csv search results from PubChem *
            +-- Rosetta_Inputs
            |   +-- * All conformer .pdb files *
            |   +-- * All conformer .param files *
            +-- User_Inputs
        +-- Fragment_PDB_Matches
        +-- PDB_search_results.json
