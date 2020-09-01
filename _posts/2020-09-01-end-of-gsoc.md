---
image: /assets/blog/2020/gsoc/gsoc-end.png
tags: gsoc mdanalysis rdkit python
title: Conclusion of my GSoC adventure
---

After 3 months of hard work, I've reached the end of my journey as a Google Summer of Code student with MDAnalysis.

It has been an extremely valuable experience for me and I hope that it will inspire other people to jump on the GSoC wagon with [MDAnalysis](https://www.mdanalysis.org/) or others (like [Open Chemistry](https://www.openchemistry.org/gsoc/)), or to contribute to the community by publishing open-source code.

You can find a recap of my whole project on this MDAnalysis [blog post](https://www.mdanalysis.org/2020/08/29/gsoc-report-cbouy/).

&nbsp;
### Preview of upcoming features
---

I don't have a lot of new things to present, I've been mainly wrapping everything up before the end of the project. Hopefully you will find these features usefull (`<buzzfeed>you won't believe what's in feature #2</buzzfeed>`) !

&nbsp;
#### 1/ Drawing code

You can now directly visualize small AtomGroups RDKit-style in Jupyter notebooks. You can also generate images (svg or png) and even gif of trajectories, similarly to what I was showing in my previous post. The rendering can be customized through the `rdkit.Chem.Draw.MolDrawOptions` class (documentation [here](http://rdkit.org/docs/cppapi/structRDKit_1_1MolDrawOptions.html))

<video controls loop src="/assets/blog/2020/gsoc/rdkit-drawer.webm"></video>

&nbsp;
#### 2/ Descriptors

Descriptors can now be computed on an AtomGroup for each frames:

```python
>>> from MDAnalysis.analysis.RDKit import RDKitDescriptors
>>> u = mda.Universe.from_smiles("CCO", numConfs=3)
>>> descriptors = ["MolWt", "RadiusOfGyration"]
>>> calc = RDKitDescriptors(u.atoms, *descriptors).run()
>>> calc.results
array([[46.06900000000002, 1.161278342193013],
       [46.06900000000002, 1.175492972121405],
       [46.06900000000002, 1.173230936577319]],
      dtype=object)
```
There are 305 descriptors (not a currated list) available this way, taken from 8 RDKit modules:

* `rdkit.Chem.Descriptors`: Molecular descriptors
* `rdkit.Chem.Descriptors3D`: Descriptors derived from a molecule's 3D structure
* `rdkit.Chem.EState.EState`: Basic EState definitions
* `rdkit.Chem.EState.EState_VSA`: Hybrid EState-VSA descriptors
* `rdkit.Chem.GraphDescriptors`: Topological/topochemical descriptors
* `rdkit.Chem.Lipinski`: Lipinski parameters for molecules
* `rdkit.Chem.MolSurf`: Approximate molecular surface area descriptors
* `rdkit.Chem.rdMolDescriptors`: Molecular descriptors (some redundancies with Chem.Descriptors)

More interestingly, you can also pass your own function as long as it takes an RDKit molecule as argument:

```python
def num_conjugated_bonds(mol):
    return sum(bond.GetIsConjugated() for bond in mol.GetBonds())
>>> u = mda.Universe.from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
>>> calc = RDKitDescriptors(u.atoms, num_conjugated_bonds).run()
>>> calc.results[0,0]
12
```

This means that you can use other packages that compute descriptors from RDKit molecules, like [Mordred](https://github.com/mordred-descriptor/mordred):

```python
>>> from mordred import Calculator, descriptors
>>> u = mda.Universe.from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
>>> calc = RDKitDescriptors(u.atoms, Calculator(descriptors)).run()
>>> len(calc.results[0,0])} # number of descriptors calculated
1826
>>> from mordred import SLogP
>>> calc = RDKitDescriptors(u.atoms, SLogP.SLogP()).run()
>>> calc.results[0,0]
-1.0293
```

Interoperability saves the day: you can now mix MDAnalysis with RDKit-friendly packages 👍

&nbsp;
#### 3/ Fingerprints

Fingerprints are also accessible now, and the following are available: AtomPair, Morgan, TopologicalTorsion, RDKit (all of these in classical or hashed version) and MACCS.

```python
>>> from MDAnalysis.analysis.RDKit import get_fingerprint
>>> from rdkit import DataStructs
>>> from nglview.datafiles import PDB, XTC
>>> # load MD
>>> u = mda.Universe(PDB, XTC)
>>> elements = mda.topology.guessers.guess_types(u.atoms.names)
>>> u.add_TopologyAttr('elements', elements)
>>> ag = u.select_atoms("resname LRT")
>>> # load retinal
>>> u2 = mda.Universe.from_smiles("CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C=O)/C)/C")
>>> # tanimoto
>>> fp1 = get_fingerprint(ag, "Morgan", radius=2, hashed=True, nBits=1024)
>>> fp2 = get_fingerprint(u2.atoms, "Morgan", radius=2, hashed=True, nBits=1024)
>>> DataStructs.TanimotoSimilarity(fp1, fp2)
0.6785714285714286
```

&nbsp;
### Lessons learnt
---

Aside from learning how to use two great packages in depth, MDAnalysis and RDKit, GSoC has taught me three valuable lessons:

&nbsp;
#### A/ Untested code is or will be broken

While some people will hate to see that the new feature they implemented broke something else in the code, they should be thankful that someone wrote tests to make sure this doesn't go unnoticed. MDAnalysis uses [pytest](https://docs.pytest.org/en/stable/) for that matter and I must say that it's ten times better than the default `unittest` package. They also have a variety of continuous integration pipelines to ensure support for different OS and dependencies versions...etc.

&nbsp;
#### B/ Undocumented code will neither be used nor maintained

Whether you're coming back to your code after 3 months or you're reviewing the code of a colleague or you're a user trying to understand what a function does and how it works, having a proper documentation can go a long way. [Sphinx](https://www.sphinx-doc.org/en/master/) makes this super easy by allowing you to automatically build documentation from functions and classes docstrings. You can even add images, code examples and links to other packages Sphinx documentation.

&nbsp;
#### C/ No software development project is a one-man project

I've kept the most important one for the conclusion! Whether you've just started programming or you're an expert, having someone review your code at the end of the day is immensely valuable to identify flaws or shortcomings or to provide some help.

With that being said, this project wouldn't have been possible without the help of the following people: my amazing mentors (Irfan Alibay, Fiona Naughton and Richard Gowers) who provided extremely helpful feedback on a daily basis together with the rest of the MDAnalysis team, more specifically Oliver Beckstein, Lily Wang, Jonathan Barnoud and Tyler Reddy.  
I'm also thankful to the RDKit community that contributed to all the tutorials on the [Cookbook](https://www.rdkit.org/docs/Cookbook.html) as well as the awesome blog posts (more specifically [iwatobipen](https://iwatobipen.wordpress.com/) and [Pat Walters](https://practicalcheminformatics.blogspot.com/)), those have been a source of inspiration for me!  
Finally, I want to thank [my lab](http://chemosim.unice.fr/) for supporting me continuously during this project.

![Thank You](/assets/blog/2020/gsoc/keanu.gif)

I hope you enjoyed following my adventure so far! I'll try my best to keep the blog posts coming, although it won't be as regular.

Best,

Cédric 🤖