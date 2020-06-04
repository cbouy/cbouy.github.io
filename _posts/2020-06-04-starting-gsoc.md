---
image: /assets/blog/2020/gsoc/rdkit_interop_banner.png
tags: gsoc mdanalysis rdkit python
title: Starting GSoC with MDAnalysis
---

This week I'm starting my Google Summer of Code (GSoC) project with MDAnalysis 😀 The goal of the project is to **make 
RDKit and MDAnalysis interoperable** (i.e. you can go back and forth between an RDKit molecule and an MDAnalysis 
universe). 

It will allow users to analyse MD trajectories with RDKit, compute 3D molecular descriptors from 
conformations of ligands in a trajectory, select groups of atoms using SMARTS patterns, display images of AtomGroups in 
Jupyter notebooks...etc. The ~~sky~~ `Segmentation fault (core dumped)` is the limit!

If you have any suggestions of cool RDKit features that should be wrapped in MDAnalysis, don't hesitate to tell me about
it in the comment section!

&nbsp;  
## What is GSoC?
---

![GSoC Logo](/assets/blog/2020/gsoc/GSoC-logo.png){:style="float: right; margin: 5px; height: 150px; width: 150px;"}GSoC is an opportunity for students to be involved in open source software development. 

Every year since 2005, Google hosts this 3 month event by matching students with organizations to write code while 
getting paid. It's a first-hand experience in real-world software development which should be extremely valuable for 
people who are interested in pursuing a career in this direction, like me.

Since its begining, GSoC has gathered 15,000+ students from 109 countries and 686 organizations! This year, we are 1191 students interacting with 199 organizations!

&nbsp;  
## MDAnalysis
---

![MDA Logo](/assets/blog/2020/gsoc/mdanalysis-logo.png){:style="float: right; margin: 5px; height: 133px; width: 300px;"} MDAnalysis is a Python library used to analyze MD trajectories in many popular formats. For the 2020 edition of GSoC, we are 3 students working on completely different projects: 
* [Hugo MacDermott-Opeskin](https://hmacdope.github.io/) implementing the much-needed TNG format (a collab' with GROMACS), 
* [Yuxuan Zhuang](https://yuxuanzhuang.github.io/) serializing MDAnalysis' Universe for parallelism, 
* and yours truly, mentored by [Richard Gowers](https://twitter.com/richardjgowers), [Irfan Alibay](https://twitter.com/HighSpeedMode) and [Fiona Naughton](https://twitter.com/ExplainedByCats)

You can learn more about all these projects in [this blog post](https://www.mdanalysis.org/2020/05/20/gsoc-students/).

&nbsp;  
## Project timeline
---

There are 3 main deliverables in the *RDKit interoperability project*:
1. RDKitReader 🔵 → 🟠
2. RDKitConverter 🟠 → 🔵
3. RDKitWrapper 🟠 + 🔵 = 💡

Here are the main lines of the project:

&nbsp;  
#### 1.1  Creating a Topology with the RDKitParser
---
*From June 1st to June 12th*

* Inherit `TopologyReaderBase` and override its `parse` method
* Add topology attributes (use guessers if N.A.):
  * Atoms (name, type, mass, charge, number)
  * Bonds (indices, bond order)
  * Protein-related attrs (residues, segments, chains…)
  * Others (angles, dihedrals…)
* Create new topology attributes:
  * Aromaticity
  * Partial charges
* Tests and docs

&nbsp;  
#### 1.2  RDKit to MDA with the RDKitReader
---
*From June 15th to June 26th*

* Inherit `MemoryReader`
* Add coordinates as a numpy array
* Tests and docs
* First tutorial

At this point, users will be able to do the following:
```python
import MDAnalysis as mda
from rdkit import Chem

# RDKit to Universe
mol = Chem.MolFromMol2File("molecule.mol2", removeHs=False)
u = mda.Universe(mol)
# do something...
```

&nbsp;  
#### 2.1  Creating a topology with the RDKitConverter
---
*From June 29th to July 10th*

* Create an editable molecule with `Chem.RWMol`
* Add `Chem.Atom` atoms to the molecule
  * Atomic number
  * Charges (partial and formal)
  * MonomerInfo (residues, segments...)
  * Store other properties with `atom.SetProp(name, value)`
* Add `Chem.Bond` bonds to the molecule
  * Indices
  * Bond type
* Tests and docs

&nbsp;  
#### 2.2  MDA to RDKit with the RDKitConverter
---
*From July 13th to July 24th*

* Guess bond types and formal charges
* Sanitization and aromaticity perception
* Add coordinates
  * each frame as a conformer
  * only store a single frame by default (to avoid exceeding the memory limit)
* Tests and docs
* Second tutorial

From the user's standpoint, it will look like this:
```python
import MDAnalysis as mda

# Universe to RDKit
u = mda.Universe('topology.gro', 'trajectory.xtc')
lig = u.select_atoms('resname LIG')
mol = lig.convert_to('RDKIT')
# do something...
```

&nbsp;  
#### 3.0  Wrapper to RDKit functionalities
---
*From July 27th to August 21st*

* New guessers:
  * Aromaticity
  * Partial and formal charges
* Molecular descriptors and fingerprints
  > A wrapper to commonly used descriptors and fingerprints:
  ```python
  lig = u.select_atoms('resname LIG')
  desc = get_descriptors(lig, ['RingCount', 'TPSA', 'NumRotatableBonds', 'GETAWAY'])
  fp = get_fingerprint(lig, 'Morgan', radius=2)
  ```
* Selecting AtomGroups with SMARTS patterns
  ```python
  group = u.select_atoms('[C,S](=O)-O', smarts=True)
  ```
* Displaying small molecules as images in Jupyter notebooks
  > Why print a string representation for molecules when you can just show an image instead ?  
    (or for meme connoisseurs, the virgin `__repr__` vs the Chad `_repr_png_`)
* Adding more Readers and Writers to MDAnalysis:
  * SDF
  * MOL
  * MAE (from Maestro)
  * TPL (from Catalyst)
  * HELM (from Pistoia Alliance)
  * Fasta
  * sequence of amino-acids

These will also come with short tutorials.

&nbsp;  
&nbsp;  

Hopefully you will find all these new features useful! I'll try to document my journey and post about cool RDKit or MDAnalysis stuff I've learned once or twice every 2 weeks.

Stay safe,

Cédric 🤖
