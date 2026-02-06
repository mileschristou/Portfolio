"""
Generate a simple micelle diagram showing surfactant molecules surrounding an oil droplet.
Surfactant molecules have hydrophobic tails (zigzag lines) pointing inward toward oil,
and hydrophilic heads (circles) pointing outward toward water.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib.patches import FancyBboxPatch
from matplotlib.path import Path
import matplotlib.patches as mpatches

# Set up the figure
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
ax.set_aspect('equal')
ax.axis('off')

# Central oil droplet (smaller)
oil_radius = 0.6
oil_circle = plt.Circle((0, 0), oil_radius, color='#FFD700', alpha=0.6, zorder=1)
ax.add_patch(oil_circle)

# Add "Oil droplet" label with arrow (positioned below to avoid overlap)
ax.annotate('Oil\ndroplet', xy=(0.4, -0.4), xytext=(1.8, -2.5),
            fontsize=14, ha='center', va='center',
            arrowprops=dict(arrowstyle='->', lw=1.5, color='black'),
            fontweight='bold')

# Number of surfactant molecules around the micelle
n_molecules = 16

# Generate positions for surfactant molecules
angles = np.linspace(0, 2*np.pi, n_molecules, endpoint=False)

def draw_surfactant(ax, angle, show_na=True, oil_radius=0.6):
    """
    Draw a single surfactant molecule at the given angle.
    Hydrophobic tail points toward center, hydrophilic head points outward.
    """
    # Head position (outer)
    head_r = 3.0
    head_radius = 0.25
    head_x = head_r * np.cos(angle)
    head_y = head_r * np.sin(angle)

    # Tail starts from oil droplet edge
    tail_start_r = oil_radius

    # Tail ends at the edge of the head (connect properly)
    tail_end_r = head_r - head_radius

    # Create zigzag tail (hydrophobic chain)
    # Reduced segments for cleaner, less chaotic look
    n_segments = 5

    # Points along the radial line for the tail
    r_points = np.linspace(tail_start_r, tail_end_r, n_segments)

    # Create zigzag by alternating perpendicular offset
    # Slightly reduced amplitude for smoother appearance
    zigzag_amplitude = 0.12
    tail_x = []
    tail_y = []

    for i, r in enumerate(r_points):
        # Base point on radial line
        base_x = r * np.cos(angle)
        base_y = r * np.sin(angle)

        # Perpendicular direction (rotate by 90 degrees)
        perp_angle = angle + np.pi/2

        # Alternate offset
        offset = zigzag_amplitude * (1 if i % 2 == 0 else -1)

        tail_x.append(base_x + offset * np.cos(perp_angle))
        tail_y.append(base_y + offset * np.sin(perp_angle))

    # Add final point that connects to the center of the head circle
    connection_r = head_r  # Connect to center of head
    tail_x.append(connection_r * np.cos(angle))
    tail_y.append(connection_r * np.sin(angle))

    # Draw the zigzag tail
    ax.plot(tail_x, tail_y, 'k-', linewidth=2, zorder=2)

    # Draw hydrophilic head (circle)
    head = plt.Circle((head_x, head_y), head_radius, color='#9370DB',
                     edgecolor='black', linewidth=1.5, zorder=3)
    ax.add_patch(head)

    # Add Na+ label to some heads (fewer to avoid clutter, positioned further out)
    if show_na:
        na_x = head_x * 1.25
        na_y = head_y * 1.25
        ax.text(na_x, na_y, 'Na⁺', fontsize=11, ha='center', va='center',
               bbox=dict(boxstyle='circle', facecolor='#E6E6FA',
                        edgecolor='none', pad=0.3), zorder=4)

# Draw all surfactant molecules
# Show Na+ labels only at cardinal directions to minimize clutter
na_positions = [0, 4, 8, 12]  # Indices to show Na+ labels (top, right, bottom, left)

for i, angle in enumerate(angles):
    show_na = i in na_positions
    draw_surfactant(ax, angle, show_na=show_na, oil_radius=oil_radius)

# Add title
ax.text(0, 4.5, 'Micelle Structure: Soap Molecules Surrounding Oil',
        fontsize=16, ha='center', va='center', fontweight='bold')

# Add legend/explanation
legend_y = -4.2
# Draw purple circle for head legend (matching the actual head color)
legend_head = plt.Circle((-4.7, legend_y), 0.08, color='#9370DB',
                        edgecolor='black', linewidth=1, zorder=5)
ax.add_patch(legend_head)
ax.text(-4.5, legend_y, 'Hydrophilic head (COO⁻)', fontsize=11, ha='left', va='center')

# Draw a small zigzag for the tail legend
legend_tail_x = [-4.7, -4.65, -4.75, -4.65, -4.7]
legend_tail_y = [legend_y - 0.4, legend_y - 0.475, legend_y - 0.55, legend_y - 0.625, legend_y - 0.7]
ax.plot(legend_tail_x, legend_tail_y, 'k-', linewidth=2)
ax.text(-4.5, legend_y - 0.55, '⚊ Hydrophobic tail (C₁₂₋₁₈)', fontsize=11, ha='left', va='center')

# Add note about water environment (moved down for more spacing)
ax.text(0, -5.1, 'Water surrounds the micelle (hydrophilic heads interact with H₂O)',
        fontsize=10, ha='center', va='center', style='italic', color='#444444')

plt.tight_layout()

# Save the diagram
output_path = 'public/micelle-structure.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Micelle diagram saved to {output_path}")

# Also save a higher resolution version
output_path_hires = 'public/micelle-structure-hires.png'
plt.savefig(output_path_hires, dpi=600, bbox_inches='tight', facecolor='white')
print(f"High-res micelle diagram saved to {output_path_hires}")

plt.close()

print("\nDiagram generation complete!")
print("The diagram shows:")
print("  - Central oil droplet (yellow/gold)")
print("  - Surfactant molecules arranged radially")
print("  - Hydrophobic tails (zigzag lines) pointing toward oil")
print("  - Hydrophilic heads (purple circles) pointing toward water")
print("  - Na⁺ counterions near some heads")
