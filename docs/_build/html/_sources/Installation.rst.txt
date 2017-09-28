*******************************
Installation
*******************************

Dependencies
============
The following packages are required to use BindingSitesFromFragments:
    * Rosetta
    * RDKit

The following packages are highly recommended:
    * Gurobi
    * OpenEye Omega (or any other method to generate small molecule conformers)

Other Programs you may find useful:
    * Avogadro
    * Pymol

Rosetta
-------
Rosetta is a Macromolecular Design Suite that has ben used for all sorts of successful protein design applications.
BindingSitesFromFragments uses Rosetta to:

    1. Score pairwise interaction energies between sets of residues
    2. Match solved binding motifs into protein scaffolds

You will need access to the RosettaCommons Github repository to clone and compile Rosetta.

RDKit
-----
`RDKit <http://www.rdkit.org/docs/Overview.html>`_ is an open source chemoinformatics package that BindingSitesFromFragments
relies on to perform substructure mapping. There are several ways to install RDKit with varying degrees of pain and
suffering required. You may refer to the `installation documentaiton <http://www.rdkit.org/docs/Install.html>`_ provided
by RDKit, but I've consolidated the instructions for my recommended methods here:

    1. If you are using a Macbook, install `homebrew <https://brew.sh/>`_ and `brew install RDKIT <https://github.com/rdkit/homebrew-rdkit>`_
    2. If you are using Ubuntu, you can ``sudo apt-get install python-rdkit librdkit1 rdkit-data``
    3. Install `Anaconda <http://docs.continuum.io/anaconda/install.html>`_ and do the following::

        $ conda create -c rdkit -n my-rdkit-env rdkit
        $ source activate my-rdkit-env

    4. You can compile Boost and RDKit from source, but I've personally had a lot of trouble doing this using non-standard installations of Python/Boost.

Gubobi
------
`Gurobi <https://www.gurobi.com/index>`_ is an optimization package that BindingSitesFromFragments uses to solve for
binding sites with minimized cumulative two-body interaction energies between all residue-residue and residue-ligand
pairs. We have a UCSF license for using Gurobi on the QB3 cluster where everything is already set up, so you most likely
will not have to worry about this.

OpenEye Omega
-------------
OpenEye Omega generates conformers for you target small molecule. You do not necessarily need to use OpenEye Omega, any
method for generating small molecule conformers will do. However, we are using Omega as it comes highly recommended from
other labs.

Avogadro
--------
Avogadro is a molecular editor that makes it simple to generate and save both target molecules and fragments as Sybyl
.mol2 or .pdb files.

Pymol
------
Pymol is a molecular viewer that is extremely handly for visualizing proteins, residue clusters, small molecule fragments,
and basically anything that can be represented as a .pdb or sybyl .mol2 file.


Installation
============
Once all of the dependencies have been installed, clone the BindingSitesFromFragments repository somewhere convenient::

    mkdir path/to/project/directory
    cd path/to/project/directory
    git clone git@github.com:jaaamessszzz/BindingSitesFromFragments.git

Installation should be as simple as::

    pip install --editable ./BindingSitesFromFragments

You should now be able to use BindingSitesFromFragments through the command line. Type ``bsff --help`` to see possible
options. The accompanying repository *jaaamessszzz/BindingSiteFromFragments-Utilities* contains several scripts that you
may find to be helpful for automating menial tasks during the design process.
