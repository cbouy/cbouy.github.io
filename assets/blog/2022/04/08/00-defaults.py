"""Configure your IPython shells and Jupyter notebooks!

All the functions and variables imported or defined in this file will be available
in any ipython shell/notebook

Place this file in `~/.ipython/profile_default/startup/` and you're good to go.
If the `~/.ipython/profile_default/` folder doesn't exist, run `ipython profile create` in a shell.

---
License: Apache License 2.0
A permissive license whose main conditions require preservation of copyright and license notices.
Contributors provide an express grant of patent rights.
Licensed works, modifications, and larger works may be distributed under different terms and without source code.

Copyright 2022 Cédric Bouysset <cedric@bouysset.net>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
from contextlib import suppress


# add space as thousand separator for ints and floats
try:
    import numpy as np
except ImportError:
    _numeric_types = [int, float]
else:
    _numeric_types =  [
        int, float,
        np.uint16, np.uint32, np.uint64,
        np.int16, np.int32, np.int64,
        np.float16, np.float32, np.float64
    ]

def numbers_with_thousand_separator(self, pretty_printer, cycle):
    pretty_printer.text(f"{self:,}".replace(",", " "))

for numeric_type in _numeric_types:
    (get_ipython()
     .display_formatter.formatters['text/plain']
     .for_type(numeric_type, numbers_with_thousand_separator))

del np, numeric_type, _numeric_types, numbers_with_thousand_separator


# RDKit
with suppress(ImportError):
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Draw import IPythonConsole
    from rdkit.Chem import rdDepictor

    # configure defaults
    IPythonConsole.ipython_useSVG = True
    rdDepictor.SetPreferCoordGen(True)

    # SMILES string as RDKit mols
    def smiles_to_svg(smiles):
        """Display str types as SVG depictions of molecules"""
        SMILES_CHARS = "AaBbCcDdEeFfGgHhIiKkLlMmNnOoPpRrSsTtUuVvWXYyZ0123456789()[].-=#$:/\\+-@"
        if all(char in SMILES_CHARS for char in smiles):
            RDLogger.DisableLog("rdApp.*")
            try:
                mol = Chem.MolFromSmiles(smiles)
                svg = IPythonConsole._toSVG(mol)
            except Exception:
                svg = None
            RDLogger.EnableLog("rdApp.*")
            return svg

    (get_ipython()
     .display_formatter.formatters['image/svg+xml']
     .for_type(str, smiles_to_svg))

    del smiles_to_svg, rdDepictor


del suppress