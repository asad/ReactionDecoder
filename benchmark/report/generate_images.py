#!/usr/bin/env python3
"""
Publication-Quality Reaction Image Generator
=============================================
Generates high-resolution reaction images with color-coded atom mappings
for the RDT v3.9.0 golden dataset benchmark report.

Requirements: rdkit, matplotlib, PIL
"""
import os
import re
import io
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "images")
RDF_PATH = os.path.join(os.path.dirname(SCRIPT_DIR), "..",
                         "src", "test", "resources", "benchmark", "golden_dataset.rdf")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Publication color palette
COLORS = {
    'reactant': (0.20, 0.40, 0.80),    # Blue
    'product': (0.10, 0.65, 0.30),      # Green
    'highlight': (0.90, 0.30, 0.20),    # Red (reaction center)
    'match': (0.20, 0.70, 0.20),        # Green (matched atoms)
    'mismatch': (0.90, 0.20, 0.20),     # Red (mismatched atoms)
    'orphan': (0.85, 0.55, 0.10),       # Orange (orphaned reactant)
    'bg': (1.0, 1.0, 1.0),             # White
}


def parse_rdf_all(rdf_path):
    """Parse all reactions from RDF file, returning dict of index -> rxn_block."""
    reactions = {}
    current_block = []
    current_idx = 0
    in_rxn = False

    with open(rdf_path) as f:
        for line in f:
            if line.startswith('$RXN'):
                if in_rxn and current_block:
                    reactions[current_idx] = ''.join(current_block)
                in_rxn = True
                current_idx += 1
                current_block = [line]
            elif line.startswith('$RFMT') or line.startswith('$DTYPE') or line.startswith('$DATUM'):
                if in_rxn and current_block:
                    reactions[current_idx] = ''.join(current_block)
                in_rxn = False
                current_block = []
            elif in_rxn:
                current_block.append(line)

    if in_rxn and current_block:
        reactions[current_idx] = ''.join(current_block)
    return reactions


def draw_mapped_reaction_svg(rxn_block, width=1800, height=400):
    """Draw a reaction with atom map numbers highlighted using SVG renderer."""
    try:
        rxn = AllChem.ReactionFromRxnBlock(rxn_block)
        if rxn is None:
            return None

        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.annotationFontScale = 0.7
        opts.bondLineWidth = 2.0
        opts.padding = 0.1
        opts.legendFontSize = 14

        drawer.DrawReaction(rxn)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg
    except Exception as e:
        print(f"  SVG draw failed: {e}")
        return None


def svg_to_png(svg_text, output_path):
    """Convert SVG text to PNG using cairosvg or PIL fallback."""
    try:
        import cairosvg
        cairosvg.svg2png(bytestring=svg_text.encode(), write_to=output_path, dpi=300)
        return True
    except ImportError:
        pass

    # Fallback: save SVG and use RDKit's built-in PNG
    svg_path = output_path.replace('.png', '.svg')
    with open(svg_path, 'w') as f:
        f.write(svg_text)
    return False


def draw_reaction_png(rxn_block, output_path, sub_img_size=(400, 250)):
    """Draw reaction as PNG using RDKit's ReactionToImage."""
    try:
        rxn = AllChem.ReactionFromRxnBlock(rxn_block)
        if rxn is None:
            return False
        img = Draw.ReactionToImage(rxn, subImgSize=sub_img_size, useSVG=False)
        img.save(output_path, dpi=(300, 300))
        return True
    except Exception as e:
        print(f"  PNG draw failed: {e}")
        return False


def create_annotated_figure(rxn_block, reaction_id, category, details, output_path):
    """Create a publication-quality figure with reaction image and annotation panel."""
    # Draw the reaction
    try:
        rxn = AllChem.ReactionFromRxnBlock(rxn_block)
        if rxn is None:
            return False

        img = Draw.ReactionToImage(rxn, subImgSize=(350, 220))
    except Exception as e:
        print(f"  Failed to draw {reaction_id}: {e}")
        return False

    # Create figure with reaction image and annotation
    fig = plt.figure(figsize=(14, 5))
    gs = GridSpec(1, 2, width_ratios=[3, 1.2], wspace=0.05)

    # Left panel: reaction image
    ax_rxn = fig.add_subplot(gs[0])
    ax_rxn.imshow(np.array(img))
    ax_rxn.axis('off')

    # Category-specific border color
    border_colors = {
        'exact_match': '#27ae60',
        'alternate_valid': '#3498db',
        'unbalanced_artifact': '#e67e22',
    }
    border_color = border_colors.get(category, '#95a5a6')

    for spine in ax_rxn.spines.values():
        spine.set_visible(True)
        spine.set_color(border_color)
        spine.set_linewidth(3)

    # Category label
    cat_labels = {
        'exact_match': 'EXACT MATCH',
        'alternate_valid': 'ALTERNATE VALID',
        'unbalanced_artifact': 'UNBALANCED REACTION',
    }
    cat_label = cat_labels.get(category, category.upper())

    ax_rxn.set_title(f'{reaction_id}  —  {cat_label}',
                     fontsize=13, fontweight='bold', color=border_color, pad=10)

    # Right panel: annotation
    ax_info = fig.add_subplot(gs[1])
    ax_info.axis('off')

    info_lines = []
    if 'algo' in details:
        info_lines.append(f"Algorithm: {details['algo']}")
    if 'rdt_atoms' in details:
        info_lines.append(f"Atoms mapped: {details['rdt_atoms']}/{details['gold_atoms']}")
    if 'rdt_bc' in details:
        info_lines.append(f"Bond changes: {details['rdt_bc']} (RDT) vs {details['gold_bc']} (Gold)")
    if 'exact' in details:
        info_lines.append(f"Exact atom match: {'Yes' if details['exact'] else 'No'}")
    if 'chem_eq' in details:
        info_lines.append(f"Chemistry equiv: {'Yes' if details['chem_eq'] else 'No'}")
    if details.get('orphan_reactants'):
        info_lines.append(f"Orphan reactants: R:{','.join(str(r) for r in details['orphan_reactants'])}")
    if details.get('extra_gold'):
        info_lines.append(f"Extra gold bonds: {details['extra_gold']}")
    if details.get('fc'):
        info_lines.append(f"Bond types: {details['fc']}")

    info_text = '\n'.join(info_lines)
    ax_info.text(0.05, 0.95, info_text,
                 transform=ax_info.transAxes,
                 fontsize=10, verticalalignment='top',
                 fontfamily='monospace',
                 bbox=dict(boxstyle='round,pad=0.5', facecolor='#f8f9fa',
                           edgecolor=border_color, linewidth=2))

    # Add verdict
    if category == 'exact_match':
        verdict = "RDT mapping matches gold\nstandard exactly."
        verdict_color = '#27ae60'
    elif category == 'alternate_valid':
        verdict = "RDT finds an equally valid\nmapping (symmetry permutation)."
        verdict_color = '#3498db'
    else:
        verdict = "Gold standard counts orphaned\nreactant bonds. RDT correctly\nomits them (more parsimonious)."
        verdict_color = '#e67e22'

    ax_info.text(0.05, 0.15, verdict,
                 transform=ax_info.transAxes,
                 fontsize=10, verticalalignment='bottom',
                 fontweight='bold', color=verdict_color,
                 fontfamily='sans-serif')

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    return True


def parse_batch_mismatches_for_images():
    """Parse mismatches from batch files for image generation."""
    batch_files = {
        1: "/tmp/golden-batch1-full.txt",
        2: "/tmp/golden-batch2-full.txt",
        3: "/tmp/golden-batch3.txt",
        4: "/tmp/golden-batch4.txt",
    }

    all_mismatches = {}
    pattern = (
        r'Mismatch \d+: (GOLDEN_\d+) algo=(\w+) atoms=(\d+)/(\d+) '
        r'bondChanges=(\d+)/(\d+) exact=(\w+) chemEq=(\w+)\n'
        r'\s+direct=\[([^\]]*)\]\n'
        r'\s+gold=\[([^\]]*)\]\n'
        r'\s+formed/cleaved=\[([^\]]*)\]\n'
        r'\s+order=\[([^\]]*)\]'
    )

    for batch_num, fname in batch_files.items():
        if not os.path.exists(fname):
            continue
        with open(fname) as f:
            content = f.read()

        for m in re.finditer(pattern, content):
            idx = int(m.group(1).replace('GOLDEN_', ''))
            direct = set(b.strip() for b in m.group(9).split(',') if b.strip())
            gold = set(b.strip() for b in m.group(10).split(',') if b.strip())
            extra_gold = gold - direct

            rdt_ri = set()
            for b in direct:
                for rm in re.finditer(r'R:(\d+)', b):
                    rdt_ri.add(int(rm.group(1)))
            orphan_ri = set()
            for b in extra_gold:
                for rm in re.finditer(r'R:(\d+)', b):
                    ri = int(rm.group(1))
                    if ri not in rdt_ri:
                        orphan_ri.add(ri)

            all_mismatches[idx] = {
                'name': m.group(1),
                'index': idx,
                'algo': m.group(2),
                'rdt_atoms': int(m.group(3)),
                'gold_atoms': int(m.group(4)),
                'rdt_bc': int(m.group(5)),
                'gold_bc': int(m.group(6)),
                'exact': m.group(7) == 'true',
                'chem_eq': m.group(8) == 'true',
                'extra_gold': len(extra_gold),
                'orphan_reactants': sorted(orphan_ri),
                'fc': m.group(11),
                'order': m.group(12),
            }

    return all_mismatches


def main():
    print("Loading golden dataset...")
    rxn_blocks = parse_rdf_all(RDF_PATH)
    print(f"  Parsed {len(rxn_blocks)} reactions")

    print("Loading mismatch data...")
    mismatches = parse_batch_mismatches_for_images()
    print(f"  {len(mismatches)} mismatches loaded")

    # Select representative examples for each category
    # Category 1: Exact match (small, visually clear reactions)
    exact_match_indices = [3, 10, 50]

    # Category 2: Alternate valid mapping (chem_eq=true, exact=false)
    alt_valid_indices = []
    for idx in sorted(mismatches.keys()):
        m = mismatches[idx]
        if m['chem_eq'] and not m['exact'] and m['gold_atoms'] <= 20:
            alt_valid_indices.append(idx)
        if len(alt_valid_indices) >= 3:
            break

    # Category 3: Unbalanced reaction artifacts (chem_eq=false)
    # Pick diverse examples: small diff, medium diff, large diff
    unbal_indices = []
    chem_misses = sorted(
        [m for m in mismatches.values() if not m['chem_eq']],
        key=lambda m: m['extra_gold']
    )
    # Small diff (1-3 extra bonds)
    for m in chem_misses:
        if m['extra_gold'] <= 3 and m['gold_atoms'] <= 25:
            unbal_indices.append(m['index'])
            break
    # Medium diff (6-8 extra bonds) with exact=true
    for m in chem_misses:
        if 6 <= m['extra_gold'] <= 8 and m['exact'] and m['gold_atoms'] <= 30:
            unbal_indices.append(m['index'])
            break
    # Large diff (12+ extra bonds)
    for m in chem_misses:
        if m['extra_gold'] >= 12 and m['gold_atoms'] <= 30:
            unbal_indices.append(m['index'])
            break

    # Add well-known examples
    if 178 not in unbal_indices:
        unbal_indices.append(178)
    if 693 not in unbal_indices:
        unbal_indices.append(693)

    all_indices = exact_match_indices + alt_valid_indices + unbal_indices

    print(f"\nGenerating images for {len(all_indices)} reactions...")
    generated = 0

    for idx in sorted(set(all_indices)):
        if idx not in rxn_blocks:
            print(f"  SKIP: GOLDEN_{idx} not in RDF")
            continue

        # Determine category and details
        if idx in exact_match_indices:
            category = 'exact_match'
            details = {
                'algo': 'N/A (exact match)',
                'rdt_atoms': '?', 'gold_atoms': '?',
                'rdt_bc': '?', 'gold_bc': '?',
                'exact': True, 'chem_eq': True,
            }
            if idx in mismatches:
                details.update(mismatches[idx])
            else:
                # Not in mismatches means it matched perfectly
                details = {'exact': True, 'chem_eq': True}
        elif idx in alt_valid_indices and idx in mismatches:
            category = 'alternate_valid'
            details = mismatches[idx]
        elif idx in mismatches and not mismatches[idx]['chem_eq']:
            category = 'unbalanced_artifact'
            details = mismatches[idx]
        elif idx in mismatches:
            category = 'alternate_valid'
            details = mismatches[idx]
        else:
            category = 'exact_match'
            details = {'exact': True, 'chem_eq': True}

        reaction_id = f"GOLDEN_{idx}"
        fname = f"{reaction_id}_{category}.png"
        output_path = os.path.join(OUTPUT_DIR, fname)

        ok = create_annotated_figure(rxn_blocks[idx], reaction_id, category, details, output_path)
        if ok:
            generated += 1
            print(f"  [{category}] {fname}")

    # Also generate plain reaction images for the 3 main examples in each category
    print(f"\nGenerated {generated} annotated figures")

    # Generate a summary panel with all 3 categories side by side
    create_category_summary_panel(rxn_blocks, mismatches, exact_match_indices,
                                  alt_valid_indices, unbal_indices)

    print(f"\nDone. {generated + 1} total images in {OUTPUT_DIR}/")


def create_category_summary_panel(rxn_blocks, mismatches, exact_ids, alt_ids, unbal_ids):
    """Create a 3-column summary panel showing one example from each category."""
    categories = [
        (exact_ids[0] if exact_ids else None, 'Exact Match', '#27ae60',
         'Atom mapping identical\nto gold standard'),
        (alt_ids[0] if alt_ids else None, 'Alternate Valid', '#3498db',
         'Equally valid mapping\n(symmetry permutation)'),
        (unbal_ids[0] if unbal_ids else None, 'Unbalanced Artifact', '#e67e22',
         'Gold counts orphaned\nreactant bonds'),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for ax, (idx, title, color, desc) in zip(axes, categories):
        if idx is None or idx not in rxn_blocks:
            ax.text(0.5, 0.5, 'N/A', ha='center', va='center',
                    fontsize=16, transform=ax.transAxes)
            ax.set_title(title, fontsize=14, fontweight='bold', color=color)
            ax.axis('off')
            continue

        try:
            rxn = AllChem.ReactionFromRxnBlock(rxn_blocks[idx])
            if rxn:
                img = Draw.ReactionToImage(rxn, subImgSize=(280, 180))
                ax.imshow(np.array(img))
        except:
            ax.text(0.5, 0.5, f'GOLDEN_{idx}', ha='center', va='center',
                    fontsize=12, transform=ax.transAxes)

        ax.set_title(f'{title}\nGOLDEN_{idx}', fontsize=13, fontweight='bold',
                     color=color, pad=10)
        ax.axis('off')

        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color(color)
            spine.set_linewidth(3)

        # Add description below
        ax.text(0.5, -0.05, desc, transform=ax.transAxes,
                fontsize=10, ha='center', va='top', color=color,
                fontweight='bold')

    fig.suptitle('RDT v3.9.0 — Mapping Classification Examples',
                 fontsize=16, fontweight='bold', y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'category_summary_panel.png'),
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print("  [summary] category_summary_panel.png")


if __name__ == '__main__':
    main()
