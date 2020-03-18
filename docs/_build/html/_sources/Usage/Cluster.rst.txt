*********************************
Cluster Protein-Fragment Contacts
*********************************

From the ``Compounds/`` directory, enter: ::

    bsff cluster <Compound_ID>

Where <Compound_ID> is the name of the project directory for your target molecule.

Agglomerative hierarchical clustering is performed for each fragment’s ensemble of protein interactions to isolate
highly represented protein-fragment contacts. A feature vector is generated for residue in the ensemble consisting of:
spatial position relative to the fragment, contact distance, interaction chemistry, and residue identity. A filtering
step removes hydrophobic residues (ACFGILMPVW) >4A from any fragment atom or polar residues (DEHKNQRSTY) >3.5A from any
fragment atom. Residues are first clustered by protein-fragment interaction chemistry (Hamming distance, cutoff=2,
linkage=compete) followed by spatial distribution about the fragment (Cosine distance, cutoff =
:math:`1 - cos(20 * \pi / 180)`, linkage=average) to generate clusters of residues that mediate similar interaction
types with specific parts of the target fragment.

In practice, this is done to create small pools of similar contact types that expedite downstream analysis. Clusters
are stored for each fragment in a directory called ``Cluster_Results``, where fragment subdirectories contain all
clusters in .pdb and .ag.npz (Prody) format.

Your project should now look like this: ::

    +-- Compounds
        +-- <Compound_ID>
            +-- Inputs
            +-- Transformed_Aligned_PDBs
            +-- Fragment_PDB_Matches
            +-- Cluster_Results
                +-- Fragment_1
                +-- Fragment_2
                +-- ...
            +-- PDB_search_results.json
