#!/usr/bin/env python3
"""
Generate improved molecular structure diagrams for terpenes with stereochemistry.
"""

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont
import io
import os

output_dir = "/Users/mileschristou/Desktop/portfolio/public"

# SMILES strings with proper stereochemistry
terpenes = [
    {
        "smiles": "CC(=C)C=C",
        "name": "Isoprene (C₅H₈)",
        "filename": "terpene-structures-improved.png"
    },
    {
        "smiles": "C[C@H]1CC=C(CC1)C(=C)C",  # R-limonene (d-limonene)
        "name": "Limonene (R)",
        "filename": None
    },
    {
        "smiles": "CC1=CC[C@H]2C[C@@H]1C2(C)C",  # (1S,5S)-α-pinene
        "name": "α-Pinene",
        "filename": None
    },
    {
        "smiles": "CC(=CCC[C@@](C)(C=C)O)C",  # (S)-linalool
        "name": "(S)-Linalool",
        "filename": None
    },
    {
        "smiles": "CC(=CCCC(=C)C=C)C",
        "name": "Myrcene",
        "filename": None
    },
    {
        "smiles": "CC1=CCC[C@]2(C)C3CC[C@@]([C@@H]3CC=C12)(C)C",  # Natural β-caryophyllene
        "name": "β-Caryophyllene",
        "filename": None
    }
]

print("Generating improved terpene structures with stereochemistry...")

# Parse molecules
mols = []
legends = []
for terp in terpenes:
    mol = Chem.MolFromSmiles(terp["smiles"])
    if mol:
        AllChem.Compute2DCoords(mol)
        mols.append(mol)
        legends.append(terp["name"])

# Create high-resolution image grid (3x2 layout)
# Increased size for 300 DPI web display
img_width = 500  # per molecule
img_height = 400  # per molecule

# Draw options for better rendering
drawer = rdMolDraw2D.MolDraw2DCairo(
    img_width * 3,  # 3 columns
    img_height * 2,  # 2 rows
    img_width,
    img_height
)

# Set drawing options
draw_opts = drawer.drawOptions()
draw_opts.bondLineWidth = 2.5  # Bolder lines for readability
draw_opts.fixedBondLength = 30  # Larger molecules
draw_opts.addStereoAnnotation = True  # Show stereochemistry
draw_opts.includeChiralFlagLabel = False
draw_opts.backgroundColour = (1, 1, 1, 1)  # White background
draw_opts.fontFile = "/System/Library/Fonts/Supplemental/Times New Roman.ttf"  # Times New Roman font

# Draw grid (no highlighting)
drawer.DrawMolecules(mols, legends=legends)

drawer.FinishDrawing()

# Save as PNG
png_data = drawer.GetDrawingText()
png_filename = f"{output_dir}/terpene-structures-improved.png"

with open(png_filename, 'wb') as f:
    f.write(png_data)

print(f"  ✓ Created: terpene-structures-improved.png")
print(f"    Resolution: {img_width * 3}x{img_height * 2} pixels")
print(f"    Layout: 3x2 grid")
print(f"    Features: Stereochemistry shown with wedge/dash bonds, bolder lines")

# Also create individual high-res structures
print("\nGenerating individual high-resolution structures...")

individual_terpenes = {
    "isoprene": {"smiles": "CC(=C)C=C", "name": "Isoprene (C₅H₈)"},
    "limonene": {"smiles": "C[C@H]1CC=C(CC1)C(=C)C", "name": "(R)-Limonene"},
    "alpha-pinene": {"smiles": "CC1=CC[C@H]2C[C@@H]1C2(C)C", "name": "α-Pinene"},
    "linalool": {"smiles": "CC(=CCC[C@@](C)(C=C)O)C", "name": "(S)-Linalool"},
    "myrcene": {"smiles": "CC(=CCCC(=C)C=C)C", "name": "Myrcene"},
    "beta-caryophyllene": {"smiles": "CC1=CCC[C@]2(C)C3CC[C@@]([C@@H]3CC=C12)(C)C", "name": "β-Caryophyllene"}
}

for key, data in individual_terpenes.items():
    mol = Chem.MolFromSmiles(data["smiles"])
    if mol:
        AllChem.Compute2DCoords(mol)

        drawer = rdMolDraw2D.MolDraw2DCairo(600, 500)
        draw_opts = drawer.drawOptions()
        draw_opts.bondLineWidth = 2.5
        draw_opts.fixedBondLength = 35
        draw_opts.addStereoAnnotation = True
        draw_opts.backgroundColour = (1, 1, 1, 1)
        draw_opts.fontFile = "/System/Library/Fonts/Supplemental/Times New Roman.ttf"

        # Draw molecule without highlighting
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        png_data = drawer.GetDrawingText()
        filename = f"{output_dir}/terpene-{key}-hires.png"
        with open(filename, 'wb') as f:
            f.write(png_data)
        print(f"  ✓ Created: terpene-{key}-hires.png")

print("\n✓ All high-resolution terpene structure images generated successfully!")
print(f"  Output directory: {output_dir}")
print("\nNotes:")
print("  - Stereochemistry indicated with wedge/dash bonds")
print("  - Greek letters (α, β) included in labels")
print("  - High resolution suitable for web display (300+ DPI equivalent)")
