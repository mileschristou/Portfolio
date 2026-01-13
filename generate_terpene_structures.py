#!/usr/bin/env python3
"""
Generate molecular structure diagrams for the terpenes blog post.
Outputs SVG files to the public folder.
"""

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import os

# Create output directory
output_dir = "/Users/mileschristou/Desktop/portfolio/public"
os.makedirs(output_dir, exist_ok=True)

# SMILES strings for terpenes
terpenes = {
    "isoprene": {
        "smiles": "CC(=C)C=C",
        "name": "Isoprene (C₅H₈)"
    },
    "limonene": {
        "smiles": "CC1=CCC(CC1)C(=C)C",
        "name": "Limonene"
    },
    "alpha-pinene": {
        "smiles": "CC1=CCC2CC1C2(C)C",
        "name": "α-Pinene"
    },
    "linalool": {
        "smiles": "CC(=CCCC(C)(C=C)O)C",
        "name": "Linalool"
    },
    "myrcene": {
        "smiles": "CC(=CCCC(=C)C=C)C",
        "name": "Myrcene"
    },
    "beta-caryophyllene": {
        "smiles": "CC1=CCCC(=C)C2CC(C2CC1)(C)C",
        "name": "β-Caryophyllene"
    }
}

# Generate individual structure images
print("Generating individual terpene structures...")
for key, data in terpenes.items():
    mol = Chem.MolFromSmiles(data["smiles"])
    if mol:
        AllChem.Compute2DCoords(mol)

        # Generate SVG
        drawer = Draw.MolDraw2DSVG(400, 300)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()

        # Save individual file
        filename = f"{output_dir}/terpene-{key}.svg"
        with open(filename, 'w') as f:
            f.write(svg)
        print(f"  Created: terpene-{key}.svg")

# Generate combined structure comparison image
print("\nGenerating combined terpene structures image...")
mols = []
labels = []
for key, data in terpenes.items():
    mol = Chem.MolFromSmiles(data["smiles"])
    if mol:
        AllChem.Compute2DCoords(mol)
        mols.append(mol)
        labels.append(data["name"])

# Create grid image (2x3 layout)
img = Draw.MolsToGridImage(
    mols,
    molsPerRow=3,
    subImgSize=(400, 300),
    legends=labels,
    useSVG=True
)

# Save combined image
combined_filename = f"{output_dir}/terpene-structures-combined.svg"
with open(combined_filename, 'w') as f:
    f.write(img)
print(f"  Created: terpene-structures-combined.svg")

print("\n✓ All terpene structure images generated successfully!")
print(f"  Output directory: {output_dir}")
