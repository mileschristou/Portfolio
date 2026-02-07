#!/usr/bin/env python3
"""
Generate PFAS-related diagrams for the blog post:
1. PFOA molecular structure with C-F bonds highlighted
2. TFE to PTFE polymerisation diagram
3. Environmental pathway diagram
"""

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont
import os

output_dir = "public"

def generate_pfoa_structure():
    """Generate clean PFOA structure without highlighting"""
    # PFOA: Perfluorooctanoic acid - C8HF15O2
    pfoa_smiles = "C(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    mol = Chem.MolFromSmiles(pfoa_smiles)

    # Generate 2D coordinates for cleaner layout
    AllChem.Compute2DCoords(mol)

    # Draw without highlights - clean structure
    drawer = rdMolDraw2D.MolDraw2DCairo(900, 400)
    drawer.drawOptions().addStereoAnnotation = False
    drawer.drawOptions().addAtomIndices = False
    drawer.drawOptions().bondLineWidth = 2

    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    # Save
    output_path = os.path.join(output_dir, "pfoa-structure.png")
    with open(output_path, 'wb') as f:
        f.write(drawer.GetDrawingText())

    # Add labels using PIL
    img = Image.open(output_path)
    img = img.convert('RGBA')

    # Create a new larger image with space for labels
    new_img = Image.new('RGBA', (1000, 500), (255, 255, 255, 255))
    new_img.paste(img, (50, 50))

    draw = ImageDraw.Draw(new_img)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 22)
        font_small = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 16)
    except:
        font = ImageFont.load_default()
        font_small = font

    # Add title and labels
    draw.text((500, 15), "PFOA (Perfluorooctanoic acid)", fill=(50, 50, 50), font=font, anchor="mt")
    draw.text((500, 470), "15 C-F bonds per molecule - bond strength ~485 kJ/mol", fill=(100, 100, 100), font=font_small, anchor="mt")

    new_img.save(output_path)
    print(f"Saved: {output_path}")


def generate_ptfe_polymerisation():
    """Generate TFE to PTFE polymerisation diagram with brackets around structures"""
    # TFE monomer: CF2=CF2
    tfe_smiles = "FC(F)=C(F)F"
    tfe = Chem.MolFromSmiles(tfe_smiles)

    # Short PTFE chain (3 units)
    ptfe_smiles = "FC(F)(C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)F"
    ptfe = Chem.MolFromSmiles(ptfe_smiles)

    # Create image with PIL
    img_width, img_height = 1200, 520
    img = Image.new('RGBA', (img_width, img_height), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img)

    try:
        font_large = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 28)
        font_medium = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 22)
        font_small = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 16)
        font_formula = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 24)
        font_subscript = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 14)
        font_n = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 20)
    except:
        font_large = ImageFont.load_default()
        font_medium = font_large
        font_small = font_large
        font_formula = font_large
        font_subscript = font_small
        font_n = font_medium

    # Draw TFE monomer
    drawer_tfe = rdMolDraw2D.MolDraw2DCairo(220, 220)
    drawer_tfe.drawOptions().addStereoAnnotation = False
    drawer_tfe.DrawMolecule(tfe)
    drawer_tfe.FinishDrawing()

    tfe_path = os.path.join(output_dir, "temp_tfe.png")
    with open(tfe_path, 'wb') as f:
        f.write(drawer_tfe.GetDrawingText())
    tfe_img = Image.open(tfe_path)

    # Position for TFE (with space for "n" prefix)
    tfe_x, tfe_y = 80, 150
    img.paste(tfe_img, (tfe_x, tfe_y))

    # Draw "n" before TFE monomer
    draw.text((tfe_x - 25, tfe_y + 100), "n", fill=(50, 50, 50), font=font_n)

    # Draw PTFE polymer
    drawer_ptfe = rdMolDraw2D.MolDraw2DCairo(480, 280)
    drawer_ptfe.drawOptions().addStereoAnnotation = False
    drawer_ptfe.DrawMolecule(ptfe)
    drawer_ptfe.FinishDrawing()

    ptfe_path = os.path.join(output_dir, "temp_ptfe.png")
    with open(ptfe_path, 'wb') as f:
        f.write(drawer_ptfe.GetDrawingText())
    ptfe_img = Image.open(ptfe_path)

    # Position for PTFE (with space for brackets)
    ptfe_x, ptfe_y = 620, 120
    img.paste(ptfe_img, (ptfe_x, ptfe_y))

    # Draw brackets around PTFE structure
    bracket_color = (50, 50, 50)
    bracket_width = 3
    # Left bracket [
    bracket_left_x = ptfe_x - 15
    bracket_top = ptfe_y + 20
    bracket_bottom = ptfe_y + 260
    draw.line([(bracket_left_x, bracket_top), (bracket_left_x, bracket_bottom)], fill=bracket_color, width=bracket_width)
    draw.line([(bracket_left_x, bracket_top), (bracket_left_x + 12, bracket_top)], fill=bracket_color, width=bracket_width)
    draw.line([(bracket_left_x, bracket_bottom), (bracket_left_x + 12, bracket_bottom)], fill=bracket_color, width=bracket_width)

    # Right bracket ]
    bracket_right_x = ptfe_x + 490
    draw.line([(bracket_right_x, bracket_top), (bracket_right_x, bracket_bottom)], fill=bracket_color, width=bracket_width)
    draw.line([(bracket_right_x, bracket_top), (bracket_right_x - 12, bracket_top)], fill=bracket_color, width=bracket_width)
    draw.line([(bracket_right_x, bracket_bottom), (bracket_right_x - 12, bracket_bottom)], fill=bracket_color, width=bracket_width)

    # Draw subscript "n" after right bracket
    draw.text((bracket_right_x + 8, bracket_bottom - 20), "n", fill=(50, 50, 50), font=font_n)

    # Add labels
    draw.text((tfe_x + 110, 50), "TFE Monomer", fill=(50, 50, 50), font=font_medium, anchor="mt")
    draw.text((tfe_x + 110, 80), "(Tetrafluoroethylene)", fill=(100, 100, 100), font=font_small, anchor="mt")

    # TFE formula with balanced subscripts
    # CF₂=CF₂
    x_tfe = 130
    y_tfe = 390
    draw.text((x_tfe, y_tfe), "CF", fill=(50, 50, 50), font=font_formula)
    draw.text((x_tfe + 31, y_tfe + 8), "2", fill=(50, 50, 50), font=font_subscript)
    draw.text((x_tfe + 43, y_tfe), "=CF", fill=(50, 50, 50), font=font_formula)
    draw.text((x_tfe + 90, y_tfe + 8), "2", fill=(50, 50, 50), font=font_subscript)

    # Draw arrow
    arrow_start = (350, 260)
    arrow_end = (570, 260)
    draw.line([arrow_start, arrow_end], fill=(50, 50, 50), width=3)
    draw.polygon([(570, 260), (550, 250), (550, 270)], fill=(50, 50, 50))

    # Arrow label
    draw.text((460, 225), "Addition", fill=(100, 100, 100), font=font_small, anchor="mt")
    draw.text((460, 285), "polymerisation", fill=(100, 100, 100), font=font_small, anchor="mt")

    # PTFE labels
    draw.text((ptfe_x + 245, 50), "PTFE Polymer Chain", fill=(50, 50, 50), font=font_medium, anchor="mt")
    draw.text((ptfe_x + 245, 80), "(Polytetrafluoroethylene / Teflon)", fill=(100, 100, 100), font=font_small, anchor="mt")

    # PTFE formula with balanced subscripts
    # −(CF₂−CF₂)n−
    x_ptfe = 780
    y_ptfe = 420
    draw.text((x_ptfe, y_ptfe), "−(CF", fill=(50, 50, 50), font=font_formula)
    draw.text((x_ptfe + 53, y_ptfe + 8), "2", fill=(50, 50, 50), font=font_subscript)
    draw.text((x_ptfe + 65, y_ptfe), "−CF", fill=(50, 50, 50), font=font_formula)
    draw.text((x_ptfe + 110, y_ptfe + 8), "2", fill=(50, 50, 50), font=font_subscript)
    draw.text((x_ptfe + 122, y_ptfe), ")", fill=(50, 50, 50), font=font_formula)
    draw.text((x_ptfe + 132, y_ptfe + 8), "n", fill=(50, 50, 50), font=font_subscript)
    draw.text((x_ptfe + 148, y_ptfe), "−", fill=(50, 50, 50), font=font_formula)

    # Bottom note
    draw.text((600, 495), "Fluorine atoms (green) form a protective sheath around the carbon backbone", fill=(100, 100, 100), font=font_small, anchor="mt")

    output_path = os.path.join(output_dir, "ptfe-polymerisation.png")
    img.save(output_path)

    # Clean up temp files
    os.remove(tfe_path)
    os.remove(ptfe_path)

    print(f"Saved: {output_path}")


def generate_environmental_pathway():
    """Generate PFAS environmental pathway diagram with cleaner arrow routing"""
    img_width, img_height = 1200, 750
    img = Image.new('RGBA', (img_width, img_height), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img)

    try:
        font_title = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 24)
        font_box = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 16)
        font_small = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 13)
    except:
        font_title = ImageFont.load_default()
        font_box = font_title
        font_small = font_title

    # Colors
    source_color = (220, 80, 80)      # Red - sources
    product_color = (80, 130, 200)    # Blue - products
    env_color = (100, 160, 100)       # Green - environment
    human_color = (180, 100, 180)     # Purple - human exposure
    box_outline = (80, 80, 80)

    def draw_box(x, y, w, h, text, color, subtext=None):
        draw.rectangle([x, y, x+w, y+h], fill=color, outline=box_outline, width=2)
        draw.text((x + w//2, y + h//2 - (8 if subtext else 0)), text, fill=(255, 255, 255), font=font_box, anchor="mm")
        if subtext:
            draw.text((x + w//2, y + h//2 + 12), subtext, fill=(240, 240, 240), font=font_small, anchor="mm")

    def draw_arrow(start, end, color=(80, 80, 80)):
        draw.line([start, end], fill=color, width=2)
        import math
        angle = math.atan2(end[1] - start[1], end[0] - start[0])
        arrow_len = 10
        arrow_angle = math.pi / 6
        x1 = end[0] - arrow_len * math.cos(angle - arrow_angle)
        y1 = end[1] - arrow_len * math.sin(angle - arrow_angle)
        x2 = end[0] - arrow_len * math.cos(angle + arrow_angle)
        y2 = end[1] - arrow_len * math.sin(angle + arrow_angle)
        draw.polygon([end, (x1, y1), (x2, y2)], fill=color)

    def draw_arrow_path(points, color=(80, 80, 80)):
        """Draw arrow with multiple waypoints (for routing around boxes)"""
        for i in range(len(points) - 1):
            if i == len(points) - 2:
                draw_arrow(points[i], points[i+1], color)
            else:
                draw.line([points[i], points[i+1]], fill=color, width=2)

    # Title
    draw.text((400, 30), "PFAS Environmental Pathways", fill=(50, 50, 50), font=font_title, anchor="mt")

    # Layout - more spread out to avoid arrow crossing
    # Row 1: Sources (y=80)
    draw_box(80, 80, 160, 50, "Manufacturing", source_color)
    draw_box(500, 80, 180, 50, "Firefighting Foam", source_color, "(AFFF)")

    # Row 2: Products & Direct contamination (y=180)
    draw_box(80, 180, 160, 50, "Consumer Products", product_color)
    draw_box(500, 180, 180, 50, "Soil Contamination", env_color, "airports, military bases")

    # Row 3: Waste paths (y=280)
    draw_box(50, 280, 130, 50, "Wastewater", env_color)
    draw_box(210, 280, 130, 50, "Landfills", env_color)
    draw_box(500, 280, 180, 50, "Groundwater", env_color)

    # Row 4: Water system (y=380)
    draw_box(50, 380, 180, 55, "Water Treatment", env_color, "(PFAS not removed)")
    draw_box(500, 380, 180, 50, "Rivers & Lakes", env_color)

    # Row 5: Human exposure routes (y=500)
    draw_box(100, 500, 180, 55, "Drinking Water", human_color)
    draw_box(500, 500, 180, 55, "Food Chain", human_color, "fish, crops, dairy")

    # Final: Humans (y=620)
    draw_box(300, 620, 200, 50, "Human Exposure", (150, 70, 150))

    # === ARROWS - Left side (manufacturing pathway) ===
    draw_arrow((160, 130), (160, 180))  # Manufacturing -> Products
    draw_arrow((120, 230), (100, 280))  # Products -> Wastewater
    draw_arrow((200, 230), (260, 280))  # Products -> Landfills
    draw_arrow((115, 330), (130, 380))  # Wastewater -> Treatment
    draw_arrow((140, 435), (180, 500))  # Treatment -> Drinking water

    # Landfills leaching - route around
    draw_arrow_path([(275, 330), (275, 360), (500, 360), (500, 380)])  # Landfills -> Groundwater (via Rivers level)

    # === ARROWS - Right side (AFFF pathway) ===
    draw_arrow((590, 130), (590, 180))  # AFFF -> Soil
    draw_arrow((590, 230), (590, 280))  # Soil -> Groundwater
    draw_arrow((590, 330), (590, 380))  # Groundwater -> Rivers
    draw_arrow((590, 430), (590, 500))  # Rivers -> Food chain

    # Groundwater to drinking water - route down right side then across
    draw_arrow_path([(500, 305), (420, 305), (420, 520), (280, 520)])

    # === ARROWS - To humans ===
    draw_arrow((190, 555), (350, 620))  # Drinking water -> Humans
    draw_arrow((590, 555), (450, 620))  # Food chain -> Humans

    # Legend
    draw.text((850, 100), "Legend:", fill=(50, 50, 50), font=font_box)
    draw.rectangle([850, 130, 880, 150], fill=source_color, outline=box_outline)
    draw.text((890, 140), "Sources", fill=(80, 80, 80), font=font_small, anchor="lm")
    draw.rectangle([850, 160, 880, 180], fill=product_color, outline=box_outline)
    draw.text((890, 170), "Products", fill=(80, 80, 80), font=font_small, anchor="lm")
    draw.rectangle([850, 190, 880, 210], fill=env_color, outline=box_outline)
    draw.text((890, 200), "Environment", fill=(80, 80, 80), font=font_small, anchor="lm")
    draw.rectangle([850, 220, 880, 240], fill=human_color, outline=box_outline)
    draw.text((890, 230), "Exposure routes", fill=(80, 80, 80), font=font_small, anchor="lm")

    # Note box
    draw.rectangle([830, 280, 1050, 370], fill=(255, 245, 230), outline=(200, 150, 100), width=2)
    draw.text((940, 295), "Key insight:", fill=(50, 50, 50), font=font_box, anchor="mt")
    draw.text((940, 320), "Standard water treatment", fill=(80, 80, 80), font=font_small, anchor="mt")
    draw.text((940, 337), "does NOT remove PFAS", fill=(150, 50, 50), font=font_small, anchor="mt")

    output_path = os.path.join(output_dir, "pfas-environmental-pathway.png")
    img.save(output_path)
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    print("Generating PFAS diagrams...")
    print()

    print("1. Generating PFOA structure...")
    generate_pfoa_structure()

    print("\n2. Generating PTFE polymerisation diagram...")
    generate_ptfe_polymerisation()

    print("\n3. Generating environmental pathway diagram...")
    generate_environmental_pathway()

    print("\nDone! Images saved to public/ directory")
