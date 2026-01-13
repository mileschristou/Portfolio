#!/usr/bin/env python3
"""
Generate caffeine molecular structure for the decaf coffee blog post.
"""

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

output_dir = "/Users/mileschristou/Desktop/portfolio/public"

# Caffeine SMILES: 1,3,7-trimethylxanthine
caffeine_smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

print("Generating caffeine structure...")
mol = Chem.MolFromSmiles(caffeine_smiles)

if mol:
    AllChem.Compute2DCoords(mol)

    # Generate larger SVG for better visibility
    drawer = Draw.MolDraw2DSVG(500, 400)
    draw_opts = drawer.drawOptions()
    draw_opts.fontFile = "/System/Library/Fonts/Supplemental/Times New Roman.ttf"
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # Save
    filename = f"{output_dir}/caffeine-molecule.svg"
    with open(filename, 'w') as f:
        f.write(svg)
    print(f"  Created: caffeine-molecule.svg")
    print("\n✓ Caffeine structure generated successfully!")
else:
    print("Error: Could not parse caffeine SMILES")
