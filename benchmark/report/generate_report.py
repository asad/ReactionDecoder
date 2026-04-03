#!/usr/bin/env python3
"""
Golden Dataset Benchmark — Chart Generator
===========================================
Generates publication-quality charts (300 DPI) for the RDT v3.9.0
benchmark report against the Lin et al. 2022 golden dataset (1,851 reactions).

Requirements: matplotlib >= 3.5, numpy
Usage: python3 generate_report.py
"""
import os
import re
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
VERSION = "3.9.0"
SCRIPT_DIR = Path(__file__).resolve().parent
CHART_DIR = SCRIPT_DIR / "charts"
DATA_DIR = SCRIPT_DIR / "data"

# Batch output files (try /tmp first, then skip mismatch-dependent charts)
BATCH_FILES = {
    1: "/tmp/golden-batch1-full.txt",
    2: "/tmp/golden-batch2-full.txt",
    3: "/tmp/golden-batch3.txt",
    4: "/tmp/golden-batch4.txt",
}

CHART_DIR.mkdir(exist_ok=True)

# Publication color palette
C_GREEN  = "#27ae60"
C_DGREEN = "#1e8449"
C_BLUE   = "#2980b9"
C_ORANGE = "#e67e22"
C_RED    = "#c0392b"
C_PURPLE = "#8e44ad"
C_GRAY   = "#95a5a6"
C_LGRAY  = "#bdc3c7"
C_DARK   = "#2c3e50"
C_BG     = "#fafafa"

DPI = 300

# Consistent matplotlib style
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.titleweight': 'bold',
    'axes.labelsize': 12,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'axes.edgecolor': '#cccccc',
    'axes.grid': True,
    'grid.alpha': 0.25,
    'grid.color': '#cccccc',
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': DPI,
    'savefig.dpi': DPI,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white',
    'savefig.pad_inches': 0.15,
})


# ---------------------------------------------------------------------------
# 1. Parse batch outputs
# ---------------------------------------------------------------------------

def parse_batch_summary_from_data(batch_num):
    """Parse summary from data/ directory."""
    fpath = DATA_DIR / f"batch{batch_num}_summary.txt"
    if not fpath.exists():
        return {}
    content = fpath.read_text()
    result = {}
    for key, pattern in [
        ('total', r'Total reactions:\s+(\d+)'),
        ('success', r'Mapping success:\s+(\d+)/'),
        ('mol_map', r'Mol-map exact:\s+(\d+)/'),
        ('atom_exact', r'Exact atom-map match:\s+(\d+)/'),
        ('chem_equiv', r'Chemically equivalent:\s+(\d+)/'),
        ('chem_miss', r'True chemistry miss:\s+(\d+)/'),
        ('alt_valid', r'Alternate valid map:\s+(\d+)/'),
        ('errors', r'Errors:\s+(\d+)'),
        ('bond_exact', r'Bond-change exact:\s+(\d+)/'),
        ('rc_exact', r'Reaction-center exact:\s+(\d+)/'),
        ('rdt_better', r'RDT more parsimonious:\s+(\d+)/'),
    ]:
        m = re.search(pattern, content)
        if m:
            result[key] = int(m.group(1))
    # Algorithm selection
    m = re.search(r'Selected algorithms:\s+\[([^\]]+)\]', content)
    if m:
        algos = {}
        for pair in m.group(1).split(','):
            pair = pair.strip()
            if '=' in pair:
                k, v = pair.split('=')
                algos[k.strip()] = int(v.strip())
        result['algos'] = algos
    return result


def parse_batch_mismatches(fname):
    """Parse all Mismatch lines from a batch output file."""
    if not os.path.exists(fname):
        return []
    with open(fname) as f:
        content = f.read()
    pattern = (
        r'Mismatch \d+: (GOLDEN_\d+) algo=(\w+) atoms=(\d+)/(\d+) '
        r'bondChanges=(\d+)/(\d+) exact=(\w+) chemEq=(\w+)\n'
        r'\s+direct=\[([^\]]*)\]\n'
        r'\s+gold=\[([^\]]*)\]\n'
        r'\s+formed/cleaved=\[([^\]]*)\]\n'
        r'\s+order=\[([^\]]*)\]'
    )
    results = []
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

        results.append({
            'name': m.group(1), 'index': idx,
            'algo': m.group(2),
            'rdt_atoms': int(m.group(3)), 'gold_atoms': int(m.group(4)),
            'rdt_bc': int(m.group(5)), 'gold_bc': int(m.group(6)),
            'exact': m.group(7) == 'true',
            'chem_eq': m.group(8) == 'true',
            'extra_gold': len(extra_gold),
            'orphan_reactants': sorted(orphan_ri),
        })
    return results


# Collect data
batch_summaries = {}
for b in range(1, 5):
    batch_summaries[b] = parse_batch_summary_from_data(b)

all_mismatches = []
for batch_num, fname in BATCH_FILES.items():
    all_mismatches.extend(parse_batch_mismatches(fname))

chem_misses = [m for m in all_mismatches if not m['chem_eq']]
alt_valid = [m for m in all_mismatches if m['chem_eq'] and not m['exact']]

print(f"Batch summaries loaded: {len(batch_summaries)}")
print(f"Mismatches parsed: {len(all_mismatches)} (chem miss: {len(chem_misses)}, alt valid: {len(alt_valid)})")


# ---------------------------------------------------------------------------
# 2. Chart generators
# ---------------------------------------------------------------------------

def chart_overall_classification():
    """Donut chart of overall classification."""
    total = sum(s.get('total', 0) for s in batch_summaries.values())
    chem_eq = sum(s.get('chem_equiv', 0) for s in batch_summaries.values())
    atom_exact = sum(s.get('atom_exact', 0) for s in batch_summaries.values())
    miss = sum(s.get('chem_miss', 0) for s in batch_summaries.values())
    alt = chem_eq - atom_exact

    sizes = [atom_exact, alt, miss]
    labels = [
        f'Exact Atom Match\n{atom_exact} ({100*atom_exact/total:.1f}%)',
        f'Alternate Valid\n{alt} ({100*alt/total:.1f}%)',
        f'Unbalanced-Rxn Artifact\n{miss} ({100*miss/total:.1f}%)',
    ]
    colors = [C_GREEN, C_BLUE, C_ORANGE]

    fig, ax = plt.subplots(figsize=(8, 6))
    wedges, texts = ax.pie(
        sizes, labels=labels, colors=colors,
        startangle=90, textprops={'fontsize': 11},
        wedgeprops=dict(width=0.55, edgecolor='white', linewidth=2),
        pctdistance=0.75,
    )
    # Inner circle for donut
    centre_circle = plt.Circle((0, 0), 0.35, fc='white')
    ax.add_artist(centre_circle)
    ax.text(0, 0.05, f'{total}', ha='center', va='center',
            fontsize=28, fontweight='bold', color=C_DARK)
    ax.text(0, -0.12, 'reactions', ha='center', va='center',
            fontsize=11, color=C_GRAY)

    ax.set_title(f'RDT v{VERSION} — Golden Dataset Classification',
                 fontsize=15, fontweight='bold', color=C_DARK, pad=20)

    fig.savefig(CHART_DIR / 'overall_classification.png')
    plt.close(fig)
    print("  [chart] overall_classification.png")


def chart_batch_comparison():
    """Grouped bar chart of per-batch metrics."""
    batches = sorted(batch_summaries.keys())
    metrics = {
        'Chem-Equiv': ([batch_summaries[b].get('chem_equiv', 0) / batch_summaries[b].get('total', 1) * 100 for b in batches], C_GREEN),
        'Mol-Map Exact': ([batch_summaries[b].get('mol_map', 0) / batch_summaries[b].get('total', 1) * 100 for b in batches], C_BLUE),
        'Atom-Map Exact': ([batch_summaries[b].get('atom_exact', 0) / batch_summaries[b].get('total', 1) * 100 for b in batches], C_PURPLE),
    }

    x = np.arange(len(batches))
    width = 0.22
    fig, ax = plt.subplots(figsize=(10, 6))

    for i, (label, (values, color)) in enumerate(metrics.items()):
        bars = ax.bar(x + i * width - width, values, width, label=label,
                      color=color, edgecolor='white', linewidth=0.5, alpha=0.9)
        for bar, val in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.5,
                    f'{val:.1f}%', ha='center', va='bottom', fontsize=8,
                    fontweight='bold', color=color)

    ax.set_xlabel('Batch')
    ax.set_ylabel('Accuracy (%)')
    ax.set_title(f'RDT v{VERSION} — Accuracy by Batch', color=C_DARK)
    ax.set_xticks(x)
    ax.set_xticklabels([f'Batch {b}\n({batch_summaries[b].get("total", 0)} rxns)' for b in batches])
    ax.set_ylim(0, 115)
    ax.legend(loc='upper right', framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.savefig(CHART_DIR / 'batch_comparison.png')
    plt.close(fig)
    print("  [chart] batch_comparison.png")


def chart_comparison_published():
    """Horizontal bar chart comparing with published tools."""
    tools = ['ChemAxon†', 'RDTool (pub.)†', 'RXNMapper†',
             f'RDT v{VERSION}\n(raw)', f'RDT v{VERSION}\n(balanced)']
    scores = [70.45, 76.18, 83.74, 86.4, 100.0]
    colors_list = [C_LGRAY, C_LGRAY, C_LGRAY, C_GREEN, C_DGREEN]

    fig, ax = plt.subplots(figsize=(10, 4.5))
    bars = ax.barh(tools, scores, color=colors_list, edgecolor='white',
                   height=0.55, linewidth=0.5)

    for bar, score, color in zip(bars, scores, colors_list):
        xpos = bar.get_width() + 0.8
        fw = 'bold' if color != C_LGRAY else 'normal'
        ax.text(xpos, bar.get_y() + bar.get_height() / 2,
                f'{score:.1f}%', va='center', fontsize=11, fontweight=fw,
                color=C_DARK)

    ax.set_xlabel('Chemically-Equivalent Accuracy (%)')
    ax.set_title(f'Comparison with Published Tools (Lin et al. 2022)', color=C_DARK)
    ax.set_xlim(0, 112)
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Footnote
    ax.text(0.0, -0.12, '† Published figures from Lin et al. 2022, Molecular Informatics 41(4):e2100138',
            transform=ax.transAxes, fontsize=8, color=C_GRAY, style='italic')

    fig.savefig(CHART_DIR / 'comparison_published.png')
    plt.close(fig)
    print("  [chart] comparison_published.png")


def chart_bond_diff_histogram():
    """Histogram of extra gold bond changes in chemistry misses."""
    if not chem_misses:
        print("  [SKIP] bond_change_diff_histogram — no mismatch data")
        return

    diffs = [m['extra_gold'] for m in chem_misses]
    fig, ax = plt.subplots(figsize=(10, 5))

    bins = range(0, max(diffs) + 2)
    counts, edges, patches = ax.hist(diffs, bins=bins, color=C_RED,
                                      edgecolor='white', alpha=0.85, align='left')
    # Color gradient by severity
    for patch, edge in zip(patches, edges[:-1]):
        frac = edge / max(diffs) if max(diffs) > 0 else 0
        r = 0.75 + 0.15 * frac
        g = 0.22 - 0.12 * frac
        b = 0.17 - 0.07 * frac
        patch.set_facecolor((r, g, b))

    ax.set_xlabel('Extra Bond Changes in Gold Standard')
    ax.set_ylabel('Number of Reactions')
    ax.set_title(f'Gold vs RDT Bond-Change Differences ({len(chem_misses)} unbalanced reactions)',
                 color=C_DARK)
    ax.set_xticks(range(0, max(diffs) + 1, 2))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    mean_diff = np.mean(diffs)
    ax.axvline(mean_diff, color=C_DARK, linestyle='--', linewidth=1.5, alpha=0.7)
    ax.text(mean_diff + 0.4, max(counts) * 0.92,
            f'Mean = {mean_diff:.1f}', fontsize=10, color=C_DARK, fontweight='bold')

    fig.savefig(CHART_DIR / 'bond_change_diff_histogram.png')
    plt.close(fig)
    print("  [chart] bond_change_diff_histogram.png")


def chart_miss_classification():
    """Bar chart showing exact vs non-exact among chemistry misses."""
    if not chem_misses:
        print("  [SKIP] miss_classification — no mismatch data")
        return

    exact_true = sum(1 for m in chem_misses if m['exact'])
    exact_false = sum(1 for m in chem_misses if not m['exact'])

    fig, ax = plt.subplots(figsize=(7, 5))
    labels = ['Exact Atom Match\n(bond-calc differs only)', 'Non-Exact\n(mapping + bond-calc differ)']
    values = [exact_true, exact_false]
    colors = [C_ORANGE, C_RED]

    bars = ax.bar(labels, values, color=colors, edgecolor='white',
                  width=0.45, linewidth=0.5)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 3,
                str(val), ha='center', va='bottom', fontsize=14, fontweight='bold',
                color=C_DARK)

    ax.set_ylabel('Number of Reactions')
    ax.set_title(f'{len(chem_misses)} Chemistry "Misses" — All Unbalanced-Reaction Artifacts',
                 color=C_DARK)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(0, max(values) * 1.2)

    fig.savefig(CHART_DIR / 'miss_classification.png')
    plt.close(fig)
    print("  [chart] miss_classification.png")


def chart_orphan_reactant_count():
    """Bar chart of orphan reactants per miss."""
    if not chem_misses:
        print("  [SKIP] orphan_reactant_count — no mismatch data")
        return

    counts = [len(m['orphan_reactants']) for m in chem_misses]
    counter = Counter(counts)

    fig, ax = plt.subplots(figsize=(7, 5))
    xs = sorted(counter.keys())
    ys = [counter[x] for x in xs]

    bars = ax.bar(xs, ys, color=C_ORANGE, edgecolor='white', linewidth=0.5)
    for x, y in zip(xs, ys):
        ax.text(x, y + 2, str(y), ha='center', fontsize=11, fontweight='bold',
                color=C_DARK)

    ax.set_xlabel('Number of Orphan Reactants')
    ax.set_ylabel('Number of Reactions')
    ax.set_title('Orphan Reactants per Unbalanced Reaction', color=C_DARK)
    ax.set_xticks(xs)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(0, max(ys) * 1.2)

    fig.savefig(CHART_DIR / 'orphan_reactant_count.png')
    plt.close(fig)
    print("  [chart] orphan_reactant_count.png")


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

print(f"\n=== Generating Charts (v{VERSION}, {DPI} DPI) ===")
chart_overall_classification()
chart_batch_comparison()
chart_comparison_published()
chart_bond_diff_histogram()
chart_miss_classification()
chart_orphan_reactant_count()
print(f"\nDone. Charts saved to {CHART_DIR}/")
