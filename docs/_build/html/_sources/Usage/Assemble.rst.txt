**********************
Assemble Contact Pools
**********************

From the ``Compounds/`` directory, enter: ::

    bsff assemble <Compound_ID>

Where <Compound_ID> is the name of the project directory for your target molecule. The fuzzballs assembled at this step
are used to solve for composite binding sites. If you are using contact pools to directly redesign the context of an
existing protein-ligand complex, you can skip to :ref:`existingcomplexfuzzball`.

.. autofunction:: BindingSitesFromFragments.commands.bsff_assemble.assemble

Your project should now look like this: ::

    +-- Compounds
        +-- <Compound_ID>
            +-- Inputs
            +-- Transformed_Aligned_PDBs
            +-- Fragment_PDB_Matches
            +-- Cluster_Results
            +-- Fuzzballs
                +-- Iteration-0
                    +-- XXX_0001-iter_0-fuzz_0.ag.npz
                    +-- XXX_0001-iter_0-fuzz_0.pdb
                    +-- XXX_0001-iter_0-fuzz_0.csv
                    +-- XXX_0002-iter_0-fuzz_0.ag.npz
                    +-- XXX_0002-iter_0-fuzz_0.pdb
                    +-- XXX_0002-iter_0-fuzz_0.csv
                    +-- ..
            +-- PDB_search_results.json
