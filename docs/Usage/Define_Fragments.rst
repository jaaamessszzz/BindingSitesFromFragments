*******************************
Define Fragments
*******************************

Generate Fragment PDB Files
===========================
Use a molecular editor like Avogadro to generate fragments of your target molecule. It is important to use the PDB named
XXX_0001.pdb created by to generate your fragments! BSFF will be using this PDB as a reference throughout the protocol.
In addition, all atoms in this PDB have unique identifiers that are essential for mapping transformations in later steps.

General rules for defining fragments:

    1. **SHOULD BE RIGID** i.e. rings, aromatic, conjugated systems. An extra step is required for fragments that contain rotatable bonds, described below.
    2. Consist of 5-10 atoms
    3. Contain a maximum of three hydrogen bond donors/acceptors

Simply delete atoms from the molecule until you end up with your desired fragment. Save as a .pdb file, then rinse and repeat
until you have defined fragments for you entire target molecule. I find anywhere between 5-8 fragments is sufficient for
target molecules that are 600Da or less.

You should name each of your files ``Fragment_<number>.pdb``, where <number> enumerates each fragment.

Move all of your fragments to the Fragment_Inputs directory so that your project looks like this: ::

    +-- Compounds
    +-- <Compound_ID>
        +-- Inputs
            +-- Fragment_Inputs
            |   +-- Fragment_Inputs.csv
            |   +-- * All fragment .pdb files *
            +-- Rosetta_Inputs
            |   +-- * All conformer .pdb files *
            |   +-- * All conformer .param files *
            +-- User_Inputs

Fragments with Rotatable Bonds
==============================
You will need to tell BindingSitesFromFragments if your fragment contains rotatable bonds. Otherwise, cluster analysis
in later steps will return garbage. To do this, create a new directory under ``Fragment_Inputs/`` called ``Rigid_Fragment_Atoms``
and populate it with the rigid portions of your fragments. Rename the .pdb files to ``Fragment_<number>-rigid.pdb``, where
``<number>`` corresponds to the original fragment. BindingSitesFromFragments will use these atoms for alignments.

Finding Fragment-Containing Small Molecules
===========================================
Once you have defined your fragments, we will use the `PubChem Open Chemistry Database <https://pubchem.ncbi.nlm.nih.gov/search/>`_
to find small molecules that contain your defined fragments as substructures. We will specifically be using the Substructure
Search tool located under Compounds > Structure > Substructure. Click on the litte hexagon at the left of the textbox to bring up
the `<PubChem Sketcher https://pubchem.ncbi.nlm.nih.gov/sketch/sketchhelp.html>`_. This tool allows you to import your
fragment .pdb and specify the characteristics and connectivities of each atom of your fragment in the substructure search results.

First, import your fragment into the PubChem Skether using the Import panel located in the bottom left of the pop-up window.
You need to click import after selecting your fragment. Your fragment will now be imported into the sketcher, where atoms
that do not have full valency represent where covalent bonds were broken to create your fragment. You will need to tell
PubChem the allowed chemical environment for these atoms (e.g. what are these atoms allowed to be bonded to) using the Query tool
(**Qry** in the top row of the panel). Click the query tool and then click atoms in your fragment to specify attributes.

The attributes I find useful are:

    * **Query flags:** use this to specify where this fragment atom is allowed to be mapped onto superstructures. For instance, check aromatic if the selected atom should always be mapped to an aromatic ring, or if an atom should always be part of an aliphatic chain.
    * **Allowed Substituents:** use this to specify how many heavy atom neighbors an atom should have when the fragment is mapped to a superstructure. This is important to ensure the molecules returned by PubChem contain your fragment in the same chemical environment as your target molecule.
    * **Ring sizes:** use this to specify ring size if your fragment is a portion of a ring system.

After you finish specifying atom attributes, return to the main PubChem search page. To get all small molecules that we are most
likely to find bound to proteins in the PDB:

    1. Click "Search Options" and uncheck *Remove any explicit hydrogens before searching*. If your fragment contains
    stereochemistry, also enable exact sterochemistry.
    2. Click "Search All" where it displays the search result count
    3. Click on "Filters" and enter an upper limit of ~500Da. This removes large molecules we are unlikely to find in the PDB.
    In addition, we are interested in finding small molecules where our fragments contribute significantly to binding affinity,
    which is not necessarily the case in extremely large molecules.
    4. Sort by "Annotation Hit Count" and click the arrow pointing down to sort most annotated to least annotated. If all goes
    well the search results will have actual names instead of numeric IDs.
    5. Click "Detailed" to bring up detailed descriptions of molecules in the search results. This is required so that
    the .csv we are about to download will contain the InChiKey for each molecule in the search results.
    6. Click "Download". Save the .csv with the same name as the input fragment in the ``Fragment_Inputs`` directory, for
    example ``Fragment_<number>.csv``.

The .csv will contain the top 1,000,000 search results for small molecules that contain our fragment as a substructure with
the same connectivites and local chemical environments as our target molecule.

Add each fragment number and its SMILES string to ``Fragment_Inputs.csv``. BindingSitesFromFragments will read this file
to search the PDB for protein-ligand complexes where the ligand contains a fragment as a substructure. When you're
finished generating fragments, your ``Fragment_Inputs.csv`` file should look something like this:

============== ==================
Fragment       SMILES_fragment
============== ==================
<number>       <SMILES_string>
<number>       <SMILES_string>
<number>       <SMILES_string>
============== ==================

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

That's about it for inputs (despite an untouched directory called ``User_Inputs``...)! We can now proceed to searching for
PDBs that contain proteins bound to small molecules in the PubChem search results for each of your defined fragments.