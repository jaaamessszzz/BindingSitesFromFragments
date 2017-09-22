# Binding-Sites-From-Fragments

<p>This is a Python3 package developed by James Lucas and friends in the Kortemme lab at UCSF for generating composite
small molecule binding sites.</p>

## Requirements
Other than the requirements listed in requirements.txt, this package also requires:
* RDKit (installing through Homebrew is the easiest)
* MySQL or SQLite (MySQL preferred)

<b>Optional:</b>
* Gurobi (highly recommended)
* Openeye Omega (or some other way for generating small molecule conformers)

## Installation
First clone the repository:
<pre><code>git clone git@github.com:jaaamessszzz/BindingSitesFromFragments.git</code></pre>

Then installation (should be) as simple as:
<pre><code>pip3 install --editable ./Binding-Sites-From-Fragments</code></pre>

## Usage
<pre><code>bsff --help</code></pre>

This command lists all available functions.