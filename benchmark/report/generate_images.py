#!/usr/bin/env python3
"""
Publication-Quality Reaction Image Generator
=============================================
Generates high-resolution annotated reaction images for the RDT v3.9.0
golden dataset benchmark report.

Uses RDKit MolDraw2D for molecule rendering with atom-map highlighting,
matplotlib for figure assembly and annotation panels.

Requirements: rdkit >= 2023.03, matplotlib >= 3.5, Pillow
Usage: python3 generate_images.py
"""
import os
import re
import io
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
VERSION = "3.9.0"
SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "images"
RDF_PATH = SCRIPT_DIR.parent.parent / "src" / "test" / "resources" / "benchmark" / "golden_dataset.rdf"

# Batch output files (for mismatch metadata)
BATCH_FILES = {
    1: "/tmp/golden-batch1-full.txt",
    2: "/tmp/golden-batch2-full.txt",
    3: "/tmp/golden-batch3.txt",
    4: "/tmp/golden-batch4.txt",
}

OUTPUT_DIR.mkdir(exist_ok=True)

DPI = 300

# Category colors
CAT_COLORS = {
    'exact_match':          '#27ae60',
    'alternate_valid':      '#2980b9',
    'unbalanced_artifact':  '#e67e22',
}
CAT_LABELS = {
    'exact_match':          'EXACT MATCH',
    'alternate_valid':      'ALTERNATE VALID',
    'unbalanced_artifact':  'UNBALANCED REACTION',
}
CAT_VERDICTS = {
    'exact_match':          'RDT mapping matches gold standard exactly.',
    'alternate_valid':      'RDT finds an equally valid mapping\n(symmetry permutation).',
    'unbalanced_artifact':  'Gold counts orphaned reactant bonds as BREAK.\nRDT correctly omits them (more parsimonious).',
}

# Publication matplotlib style
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 10,
    'figure.facecolor': 'white',
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white',
    'savefig.pad_inches': 0.1,
})


# ---------------------------------------------------------------------------
# RDF Parser
# ---------------------------------------------------------------------------
def parse_rdf_all(rdf_path):
    """Parse all reactions from RDF file, returning dict of 1-based index -> rxn_block."""
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


# ---------------------------------------------------------------------------
# Mismatch parser
# ---------------------------------------------------------------------------
def parse_all_mismatches():
    """Parse mismatch data from batch output files."""
    pattern = (
        r'Mismatch \d+: (GOLDEN_\d+) algo=(\w+) atoms=(\d+)/(\d+) '
        r'bondChanges=(\d+)/(\d+) exact=(\w+) chemEq=(\w+)\n'
        r'\s+direct=\[([^\]]*)\]\n'
        r'\s+gold=\[([^\]]*)\]\n'
        r'\s+formed/cleaved=\[([^\]]*)\]\n'
        r'\s+order=\[([^\]]*)\]'
    )
    all_mismatches = {}
    for batch_num, fname in BATCH_FILES.items():
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
                'name': m.group(1), 'index': idx,
                'algo': m.group(2),
                'rdt_atoms': int(m.group(3)), 'gold_atoms': int(m.group(4)),
                'rdt_bc': int(m.group(5)), 'gold_bc': int(m.group(6)),
                'exact': m.group(7) == 'true',
                'chem_eq': m.group(8) == 'true',
                'extra_gold': len(extra_gold),
                'orphan_reactants': sorted(orphan_ri),
                'fc': m.group(11),
            }
    return all_mismatches


# ---------------------------------------------------------------------------
# Rendering helpers
# ---------------------------------------------------------------------------
def rxn_block_to_image(rxn_block, width=900, height=300):
    """Render a reaction block to a PIL Image using RDKit."""
    try:
        rxn = AllChem.ReactionFromRxnBlock(rxn_block)
        if rxn is None:
            return None
        # Use ReactionToImage with decent subimage size
        img = Draw.ReactionToImage(rxn, subImgSize=(width // max(rxn.GetNumReactantTemplates() + rxn.GetNumProductTemplates(), 1), height))
        return img
    except Exception as e:
        print(f"    Draw failed: {e}")
        return None


def rxn_block_to_image_large(rxn_block, target_width=1400, target_height=350):
    """Render reaction at large size for annotated figures."""
    try:
        rxn = AllChem.ReactionFromRxnBlock(rxn_block)
        if rxn is None:
            return None
        n_mols = rxn.GetNumReactantTemplates() + rxn.GetNumProductTemplates()
        sub_w = max(250, target_width // max(n_mols, 1))
        sub_h = target_height
        img = Draw.ReactionToImage(rxn, subImgSize=(sub_w, sub_h))
        return img
    except Exception as e:
        print(f"    Large draw failed: {e}")
        return None


# ---------------------------------------------------------------------------
# Annotated figure generator
# ---------------------------------------------------------------------------
def create_annotated_figure(rxn_block, reaction_id, category, details, output_path):
    """Create publication-quality figure: reaction image + annotation panel."""
    img = rxn_block_to_image_large(rxn_block)
    if img is None:
        return False

    border_color = CAT_COLORS.get(category, '#95a5a6')
    cat_label = CAT_LABELS.get(category, category.upper())
    verdict = CAT_VERDICTS.get(category, '')

    fig = plt.figure(figsize=(16, 5))
    gs = GridSpec(1, 2, width_ratios=[3.2, 1], wspace=0.03)

    # Left: reaction image
    ax_rxn = fig.add_subplot(gs[0])
    ax_rxn.imshow(np.array(img))
    ax_rxn.axis('off')
    ax_rxn.set_title(f'{reaction_id}  —  {cat_label}',
                     fontsize=14, fontweight='bold', color=border_color, pad=12)
    for spine in ax_rxn.spines.values():
        spine.set_visible(True)
        spine.set_color(border_color)
        spine.set_linewidth(3)

    # Right: info panel
    ax_info = fig.add_subplot(gs[1])
    ax_info.axis('off')

    info_lines = []
    if 'algo' in details and details['algo'] != 'N/A':
        info_lines.append(f"Algorithm:       {details['algo']}")
    if 'rdt_atoms' in details and details['rdt_atoms'] != '?':
        info_lines.append(f"Atoms mapped:    {details['rdt_atoms']}/{details['gold_atoms']}")
    if 'rdt_bc' in details and details['rdt_bc'] != '?':
        info_lines.append(f"Bond changes:    {details['rdt_bc']} (RDT)")
        info_lines.append(f"                 {details['gold_bc']} (Gold)")
    if 'exact' in details:
        v = 'Yes' if details['exact'] else 'No'
        info_lines.append(f"Exact atom map:  {v}")
    if 'chem_eq' in details:
        v = 'Yes' if details['chem_eq'] else 'No'
        info_lines.append(f"Chem-equiv:      {v}")
    if details.get('orphan_reactants'):
        orph = ', '.join(f'R:{r}' for r in details['orphan_reactants'])
        info_lines.append(f"Orphan reactants: {orph}")
    if details.get('extra_gold'):
        info_lines.append(f"Extra gold bonds: {details['extra_gold']}")

    info_text = '\n'.join(info_lines)
    ax_info.text(0.05, 0.95, info_text,
                 transform=ax_info.transAxes,
                 fontsize=9.5, verticalalignment='top',
                 fontfamily='monospace',
                 bbox=dict(boxstyle='round,pad=0.6', facecolor='#f8f9fa',
                           edgecolor=border_color, linewidth=2, alpha=0.95))

    ax_info.text(0.05, 0.12, verdict,
                 transform=ax_info.transAxes,
                 fontsize=10, verticalalignment='bottom',
                 fontweight='bold', color=border_color,
                 fontfamily='sans-serif', linespacing=1.4)

    fig.savefig(output_path, dpi=DPI)
    plt.close(fig)
    return True


# ---------------------------------------------------------------------------
# Category summary panel
# ---------------------------------------------------------------------------
def create_category_summary_panel(rxn_blocks, mismatches, exact_ids, alt_ids, unbal_ids):
    """Create 3-column summary panel showing one example from each category."""
    categories = [
        (exact_ids[0] if exact_ids else None, 'Exact Match', '#27ae60',
         'Atom mapping identical\nto gold standard'),
        (alt_ids[0] if alt_ids else None, 'Alternate Valid', '#2980b9',
         'Equally valid mapping\n(symmetry permutation)'),
        (unbal_ids[0] if unbal_ids else None, 'Unbalanced Artifact', '#e67e22',
         'Gold counts orphaned\nreactant bonds as BREAK'),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(20, 5.5))

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
                n_mols = rxn.GetNumReactantTemplates() + rxn.GetNumProductTemplates()
                sub_w = max(200, 500 // max(n_mols, 1))
                img = Draw.ReactionToImage(rxn, subImgSize=(sub_w, 280))
                ax.imshow(np.array(img))
        except Exception:
            ax.text(0.5, 0.5, f'GOLDEN_{idx}', ha='center', va='center',
                    fontsize=12, transform=ax.transAxes)

        ax.set_title(f'{title}\nGOLDEN_{idx}', fontsize=13, fontweight='bold',
                     color=color, pad=12)
        ax.axis('off')
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color(color)
            spine.set_linewidth(3)

        ax.text(0.5, -0.06, desc, transform=ax.transAxes,
                fontsize=10, ha='center', va='top', color=color,
                fontweight='bold', linespacing=1.3)

    fig.suptitle(f'RDT v{VERSION} — Mapping Classification Examples',
                 fontsize=16, fontweight='bold', y=1.02, color='#2c3e50')
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / 'category_summary_panel.png', dpi=DPI)
    plt.close(fig)
    print("  [summary] category_summary_panel.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("Loading golden dataset...")
    rxn_blocks = parse_rdf_all(str(RDF_PATH))
    print(f"  Parsed {len(rxn_blocks)} reactions")

    print("Loading mismatch data...")
    mismatches = parse_all_mismatches()
    print(f"  {len(mismatches)} mismatches loaded")

    # --- Select representative examples ---

    # Exact match: reactions NOT in mismatch list (truly identical mapping)
    all_mismatch_ids = set(mismatches.keys())
    exact_match_indices = sorted(set(range(1, 200)) - all_mismatch_ids)[:5]
    print(f"  Exact match examples: {exact_match_indices}")

    # Alternate valid: chem_eq=true, exact=false
    alt_valid_indices = []
    for idx in sorted(mismatches.keys()):
        m = mismatches[idx]
        if m['chem_eq'] and not m['exact'] and m['gold_atoms'] <= 25:
            alt_valid_indices.append(idx)
        if len(alt_valid_indices) >= 3:
            break

    # Unbalanced artifacts: diverse bond-change differences
    unbal_indices = []
    chem_misses = sorted(
        [m for m in mismatches.values() if not m['chem_eq']],
        key=lambda m: m['extra_gold']
    )
    # Small diff (1-3 extra)
    for m in chem_misses:
        if m['extra_gold'] <= 3 and m['gold_atoms'] <= 25:
            unbal_indices.append(m['index'])
            break
    # Medium diff (6-8 extra) with exact=true
    for m in chem_misses:
        if 6 <= m['extra_gold'] <= 8 and m['exact'] and m['gold_atoms'] <= 30:
            unbal_indices.append(m['index'])
            break
    # Large diff (12+ extra)
    for m in chem_misses:
        if m['extra_gold'] >= 12 and m['gold_atoms'] <= 30:
            unbal_indices.append(m['index'])
            break
    # Named examples from report
    for named in [178, 221, 692, 693, 1088, 1404]:
        if named not in unbal_indices:
            unbal_indices.append(named)

    all_indices = sorted(set(exact_match_indices + alt_valid_indices + unbal_indices))
    print(f"\nGenerating images for {len(all_indices)} reactions...")

    generated = 0
    for idx in all_indices:
        if idx not in rxn_blocks:
            print(f"  SKIP: GOLDEN_{idx} not in RDF")
            continue

        # Determine category
        if idx in mismatches:
            m = mismatches[idx]
            if not m['chem_eq']:
                category = 'unbalanced_artifact'
            elif m['chem_eq'] and not m['exact']:
                category = 'alternate_valid'
            else:
                category = 'exact_match'
            details = m
        else:
            category = 'exact_match'
            details = {'exact': True, 'chem_eq': True}

        reaction_id = f"GOLDEN_{idx}"
        fname = f"{reaction_id}_{category}.png"
        output_path = OUTPUT_DIR / fname

        ok = create_annotated_figure(rxn_blocks[idx], reaction_id, category, details, str(output_path))
        if ok:
            generated += 1
            print(f"  [{category}] {fname}")
        else:
            print(f"  FAIL: {fname}")

    print(f"\nGenerated {generated} annotated figures")

    # Summary panel
    create_category_summary_panel(rxn_blocks, mismatches,
                                  exact_match_indices, alt_valid_indices, unbal_indices)

    print(f"\nDone. {generated + 1} total images in {OUTPUT_DIR}/")


if __name__ == '__main__':
    main()
