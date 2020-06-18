---
image: /assets/blog/2020/gsoc/rdkitreader.png
tags: gsoc mdanalysis rdkit python
title: From RDKit to the Universe
---

After almost 2 weeks since my last post, it's time for me to report on all the cool things I've learned and the features I've developed so far!


&nbsp;
### New features
---
I added 3 features to MDAnalysis since the begining of my GSoC project:

#### 1. Creating a Universe from an RDKit molecule

Here is a quick example on how to do it:
```python
from rdkit import Chem
import MDAnalysis as mda

mol = Chem.MolFromMol2File("docking_poses.mol2", removeHs=False)
u = mda.Universe(mol)
u.trajectory
# <RDKitReader with 10 frames of 42 atoms>
```

Pretty simple, right? Now obviously MDAnalysis already supports MOL2 files, and there are other formats that both libraries support (like PDB files), so the benefit of doing this might not be apparent at first. There are actually two uses for this:
  * **Making new formats indirectly available to MDAnalysis**  
    This concerns mostly SDF and MOL files.
  * **Covering all the subtle differences when parsing files**  
    There are probably as many PDB flavors as there are Linux distributions, and chances are they will be parsed differently by RDKit and MDAnalysis. Now thanks to this feature, if the default PDB parser from MDAnalysis wasn't doing what you expected, you can ~~blame the developers~~ use RDKit's parser which will hopefully give you what you expected, and then convert the RDKit molecule to a Universe.

#### 2. Constructing a Universe from a SMILES string

To create a Universe with 3 different conformations of ethanol (don't ask me why I chose this molecule):
```python
u = mda.Universe.from_smiles("CCO", numConfs=3)

u.trajectory
# <RDKitReader with 3 frames of 9 atoms>
```

By default, this method will add hydrogens and generate coordinates for 1 conformer. Internally, this is equivalent to:
```python
mol = Chem.MolFromSmiles("CCO")
mol = Chem.AddHs(mol)
confids = AllChem.EmbedMultipleConfs(mol, numConfs=1)
u = mda.Universe(mol)
```

The `EmbedMultipleConfs` function has a lot of interesting arguments available, and they are still accessible in the `from_smiles` method through the `rdkit_kwargs` parameter. For example, to use a different conformer generation algorithm like ETKDGv3:
```python
u = mda.Universe.from_smiles("CCO", rdkit_kwargs=dict(params=AllChem.ETKDGv3()))
```
To fix the seed for conformer generation:
```python
u = mda.Universe.from_smiles("CCO", rdkit_kwargs=dict(randomSeed=42))
```

#### 3. Aromatic atoms selection

Something else that MDAnalysis wasn't able to do previously: selecting aromatic atoms. There's a very simple reason behind this: an aromaticity attribute did not exist and there was no perception of the aromatic nature of bonds and atoms in MDAnalysis. This is now a thing of the past as we can leverage RDKit's aromaticity perception ability, although currently it's only available when the Universe was created from an RDKit molecule.
```python
u = mda.Universe.from_smiles("Cc1cNcc1")

u.select_atoms("aromatic")
# <AtomGroup with 5 atoms>

u.select_atoms("not aromatic")
# <AtomGroup with 1 atom>

u.select_atoms("type N and not aromatic")
# <AtomGroup with 0 atoms>
```

Support for all formats will come in the late stage of the project, as it depends on over things (mainly the RDKit converter).

&nbsp;

All these features are part of [PR#2707](https://github.com/MDAnalysis/mdanalysis/pull/2707) and will be available on the [develop branch](https://github.com/MDAnalysis/mdanalysis) of MDAnalysis very soon.


&nbsp;
### What I've learned
---
This first step was a good way for me to get acquainted with the MDAnalysis library as I had never used it before GSoC. I now have a good overview of how the core parts of the Universe work (topology, topology attributes and trajectory).

I also got more familiar with writing **unit tests** with pytest, which I find much more intuitive and convenient to use than Python's default library "unittest". It's really important to make sure that the code you write returns what you'd expect, but also that you don't break the previous code when adding new features. This means that a lot of time is spent writing tests that cover almost all possible scenarios: in my case I'd say I spend as much time on tests than on the actual code.

![I have no idea what I'm doing](/assets/blog/2020/gsoc/noideawhatimdoing.gif){:style="float: right; margin: 5px; width:270px; height: 166px"}Another important aspect of software development, and which I'm quite new to, is writing a proper **documentation**. 
I guess you know how frustrating it can be when you find a library that does what you're looking for but there's barely any documentation on how to use it, so you either try to guess how it works, or you just search for another solution.  
The most popular choice for documenting Python code seems to be Sphinx which makes it very easy. The most helpful features I've seen so far are the automatic inclusion of docstrings as actual documentation, and the possibility of embedding an ipython shell in the sphinx document. I've never used Sphinx before so I'm learning by making mistakes for now, but I'll make sure to spend some time following tutorials on that! 

And lastly the most important part: **code review**. Richard and Irfan have been super helpful and reactive at every commit so I'm really thankful to them for mentoring me :pray:
They are much more experienced with the library and software development in general so they quickly see things that I missed, especially for tests and docs.


&nbsp;
### What's next
---
I'll now start working on the most exciting part: converting an MDAnalysis universe to an RDKit molecule! This will open MDAnalysis to all the chemoinformatics analysis and awesome features available in RDKit.

I'm a bit in advance on the schedule that I proposed, which means that I should be able to spend more time on wrapping RDKit functionalities inside MDAnalysis during the month of August.


Until next time,

Cédric