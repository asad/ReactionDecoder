# Reaction Decoder Tool (RDT) v3.9.0
## Algorithm Description and Benchmark Evaluation

**Authors:** Syed Asad Rahman
**Affiliation:** BioInception PVT LTD
**Contact:** asad.rahman@bioinceptionlabs.com
**License:** GNU LGPL v3.0
**Version:** 3.9.0 (April 2026)

---

## Abstract

We present the Reaction Decoder Tool (RDT), a deterministic, training-free algorithm for atom-atom mapping (AAM) of chemical reactions. RDT employs a multi-algorithm ensemble over a game-theory-inspired scoring framework, combining Maximum Common Subgraph (MCS) computation with four complementary mapping heuristics (MAX, MIN, MIXTURE, RINGS) and a hierarchical 12-level solution selector. On the 1,851-reaction Lin et al. (2022) golden dataset, RDT achieves **99.2% chemically-equivalent accuracy**, outperforming all published deterministic tools (RDTool 2016: 76.18%; ChemAxon: 70.45%) and the unsupervised neural method RXNMapper (83.74%), without any training data or learned parameters.

---

## 1. Problem Statement

**Definition (Atom-Atom Mapping).** Given a chemical reaction *r* = (*R*, *P*) where *R* = {*R*₁, …, *R*_m} is a set of reactant molecules and *P* = {*P*₁, …, *P*_n} is a set of product molecules, find a bijection:

    φ : A(R) → A(P)

where *A(·)* denotes the set of heavy atoms in a molecule set, such that φ minimises the total bond change count:

    Δ(φ) = |{(a,b) ∈ E(R) : (φ(a),φ(b)) ∉ E(P)}|
           + |{(a,b) ∈ E(P) : (φ⁻¹(a),φ⁻¹(b)) ∉ E(R)}|

where *E(·)* denotes the set of bonds (edges with order label) in a molecule set.

This is NP-hard in general (reducible from graph isomorphism), so practical solvers apply heuristic decomposition over molecule pairs.

---

## 2. Algorithm Overview

RDT proceeds through nine sequential stages:

```
Input Reaction SMILES / RXN / IReaction
          │
          ▼
 ┌─────────────────────┐
 │  Stage 1            │  Parse & preprocess
 │  Parsing            │  (atom types, aromaticity, implicit H)
 └────────┬────────────┘
          │
          ▼
 ┌─────────────────────┐
 │  Stage 2            │  Reagent filter (known solvents, catalyst
 │  Standardisation    │  metals, fingerprint similarity)
 └────────┬────────────┘
          │
          ▼
 ┌─────────────────────┐
 │  Stage 3            │  RINGS funnel: test ring-conservation
 │  Quality Gate       │  mapping; exit early if coverage ≥ 95%
 └────────┬────────────┘
          │ (insufficient)
          ▼
 ┌─────────────────────┐
 │  Stage 4            │  Parallel execution of MIN, MAX,
 │  Multi-Algorithm    │  MIXTURE, RINGS algorithms
 └────────┬────────────┘
          │
          ▼
 ┌─────────────────────┐
 │  Stage 5            │  Pairwise MCS computation
 │  MCS Engine         │  (identity shortcut → substructure → VF2++)
 └────────┬────────────┘
          │
          ▼
 ┌─────────────────────┐
 │  Stage 6            │  7-matrix game-theory scoring
 │  Score Matrices     │  (clique, Jaccard, stereo, energy,
 └────────┬────────────┘   fragment, carbon, fingerprint)
          │
          ▼
 ┌─────────────────────┐
 │  Stage 7            │  Algorithm-specific winner selection
 │  Selection          │  per assignment matrix
 └────────┬────────────┘
          │
          ▼
 ┌─────────────────────┐
 │  Stage 8            │  Cross-algorithm solution ranking
 │  Best Mapping       │  (12-level comparator)
 └────────┬────────────┘
          │
          ▼
 ┌─────────────────────┐
 │  Stage 9            │  Bond change annotation,
 │  Output             │  fingerprint generation, SMILES output
 └─────────────────────┘
```

---

## 3. Stage-by-Stage Description

### 3.1 Parsing and Preprocessing

Input formats accepted: reaction SMILES (Daylight notation), RXN V2000/V3000, RDF, or CDK `IReaction` objects.

Preprocessing pipeline:
1. Null implicit hydrogen counts → 0
2. Remove all pre-existing atom-atom map numbers
3. Perceive atom types using CDK `CDKAtomTypeMatcher`
4. Perceive aromaticity using the Daylight model via CDK `Aromaticity`

---

### 3.2 Reaction Standardisation

**Purpose:** Remove non-reacting species (solvents, catalysts, reagents) to focus MCS computation on the reacting core, reducing computation and preventing spurious mappings.

**Three-tier reagent filter (conservative):**

| Tier | Method | Criterion |
|------|--------|-----------|
| 1 | Known-reagent lookup | Canonical SMILES match against a database of ~35 common solvents and inorganic salts (DCM, DMSO, DMF, THF, pyridine, NaOH, etc.) |
| 2 | Catalyst metal check | Molecule contains Pd, Pt, Rh, Ru, Ir, Ni, Cu, Fe, Co, Mn, Ti, Zr, Mo, W, Os, Ag, or Au |
| 3 | Fingerprint similarity | ECFP4 (radius=2, 256 bits) Tanimoto similarity to *all* products < 0.4, heavy atom count ≤ 10, and no element unique to products |

**Atom-balance guard:** Before any molecule is removed, verify:

    ∀ element e: |A(R \ {reactant}, e)| ≥ |A(P, e)|

If this fails (removing the candidate would unbalance the reaction), the molecule is retained regardless of tier classification.

**Annotation:** Source-occurrence identifiers (`sourceOccurrenceId`, `sourceAtomId`) are stamped on every molecule and atom before filtering. For identical-signature duplicate molecules (e.g. two equivalents of water), `preserveOccurrenceIdentity = true` ensures each occurrence maps independently.

---

### 3.3 RINGS Funnel (Quality Gate)

For reactions with ≤ 5 total molecules (reactants + products), RDT first executes the RINGS algorithm alone and evaluates the result:

**Coverage criterion:**

    coverage(φ) = |{a ∈ A(R) : φ(a) is defined}| / |A(R)|

If `coverage(φ) ≥ 0.95` and the mapping is non-trivial (the reaction contains actual structural changes), RDT returns immediately without invoking MIN, MAX, or MIXTURE. In practice this resolves approximately **75% of reactions** at the single-algorithm cost.

---

### 3.4 Parallel Multi-Algorithm Execution

When the RINGS funnel is insufficient, the remaining algorithms execute in parallel via a shared fixed-thread executor (`min(2, min(3, nCPU))` daemon threads):

| Algorithm | Selection bias | Primary objective |
|-----------|---------------|-------------------|
| **MAX** | `MaxSelection` | Maximise total mapped atoms (global coverage) |
| **MIN** | `MinSelection` | Minimise total bond changes (parsimony) |
| **MIXTURE** | Hybrid max→min | Fallback: mixed coverage/parsimony for edge cases where MinSelection suppresses a valid pairing |
| **RINGS** | Ring-conservation | Preserve ring systems and aromatic skeleton topology |

MIXTURE runs with identical MCS settings to MIN and is deduplicated at collection time. It survives deduplication only when the assignment matrices produce a distinct pairing — it serves as a genuine fallback for the subset of reactions where MinSelection is overly conservative.

---

### 3.5 Pairwise MCS Computation

For each reactant-product pair *(R_i, P_j)*, compute a Maximum Common Subgraph (MCS) mapping to establish atom correspondences.

#### 3.5.1 Three-Stage Pre-filter

**Stage 1 — Identity shortcut:**

    if canSmiles(R_i) = canSmiles(P_j) AND |A(R_i)| = |A(P_j)|
        φ_{ij} := {(a_k, a_k) : k = 1…|A(R_i)|}   (direct 1:1 mapping)
        skip MCS

Canonical SMILES are generated by `MolGraph.toCanonicalSmiles()` (SMSD 6.10.1), which encodes tetrahedral chirality (`@`/`@@`) and E/Z geometry (`/`/`\`). This is essential: using a stereo-unaware generator would incorrectly short-circuit enantiomers (e.g. (R)-lactic acid ≡ (S)-lactic acid) to a spurious identity mapping.

**Stage 2 — Size ratio filter:**

    if min(|A(R_i)|, |A(P_j)|) / max(|A(R_i)|, |A(P_j)|) < 0.3
        AND min(|A(R_i)|, |A(P_j)|) > 3
        skip pair

**Stage 3 — Fingerprint filter:**

    if Tanimoto(PathFP(R_i), PathFP(P_j)) < 0.05
        AND min(|A(R_i)|, |A(P_j)|) > 5
        skip pair

#### 3.5.2 Tiered Substructure Search

For pairs that survive pre-filtering, attempt subgraph isomorphism with progressively relaxed criteria:

```
Tier 1:  AtomType = strict CDK type
         BondOrder = flexible
         RingMatch = strict (ring bonds match ring bonds)
            │  (no subgraph found)
            ▼
Tier 2:  AtomType = element symbol only
         BondOrder = flexible
         RingMatch = strict
            │  (no subgraph found)
            ▼
Tier 3:  AtomType = element symbol only
         BondOrder = flexible
         RingMatch = relaxed
```

Each tier uses VF2++ subgraph isomorphism (SMSD engine) with a 5-second hard timeout.

#### 3.5.3 Full MCS Fallback

When no substructure relationship holds in either direction, invoke the SMSD Maximum Common Subgraph algorithm. The MCS finds the largest atom set *M* ⊆ *A(R_i)* × *A(P_j)* such that the induced subgraphs are isomorphic.

**Cache:** Results are memoised in a thread-safe LRU cache (capacity 10,000 entries) keyed by:

    key = canonSmiles(R_i) + "|" + canonSmiles(P_j) + "|" + theory + "|" + settings + "|" + fpHash

This enables cross-reaction reuse when the same molecule pair appears in multiple reactions (common in metabolic pathway datasets).

**Circular fingerprint cache:** Each molecule's FCFP (radius=1, 256 bits) is computed once and cached in an `IdentityHashMap` keyed by object identity, avoiding redundant re-computation.

---

### 3.6 Game-Theory Scoring Matrices

For each algorithm execution, construct seven *m × n* scoring matrices (where *m* = |reactants|, *n* = |products|):

| Symbol | Name | Formula |
|--------|------|---------|
| *C(i,j)* | Clique | `|MCS(R_i, P_j)|` — atom count of MCS |
| *G(i,j)* | Jaccard | `|MCS| / (|R_i| + |P_j| - |MCS|)` |
| *S(i,j)* | Stereo | Stereo compatibility score from SMSD stereo analysis |
| *E(i,j)* | Energy | Sum of bond dissociation energies over the mapped bonds (Luo 2007 BDE table) |
| *F(i,j)* | Fragment | Number of disconnected fragments in the MCS mapping |
| *K(i,j)* | Carbon | `|{a ∈ MCS(R_i,P_j) : symbol(a) = C}|` |
| *T(i,j)* | Tanimoto | `Tanimoto(PathFP(R_i), PathFP(P_j))` |

These matrices encode the multi-objective assignment problem as a 7-dimensional payoff table, analogous to a cooperative game where reactants and products are players choosing pairings.

---

### 3.7 Algorithm-Specific Assignment

Each algorithm iteratively selects the globally best reactant-product pair and extracts its atom mapping, removing the pair from the matrix until all molecules are assigned or no valid pairs remain.

**Pseudocode (MAX algorithm):**

```
function MAX_ASSIGN(C, G, m, n):
    assigned_rows ← ∅
    assigned_cols ← ∅
    mappings ← []
    while assigned_rows ≠ {1…m} AND assigned_cols ≠ {1…n}:
        best ← argmax_{i∉assigned_rows, j∉assigned_cols}
                    G(i,j)  s.t. isMajorSubgraph(C, i, j)
        if best = ∅: break
        mappings.append( MCS(R_{best.i}, P_{best.j}) )
        assigned_rows ← assigned_rows ∪ {best.i}
        assigned_cols ← assigned_cols ∪ {best.j}
    return mappings
```

where `isMajorSubgraph(C, i, j)` is true if *C(i,j)* is the maximum entry in both row *i* and column *j* simultaneously (the pair dominates all alternatives in its row and column).

**Pseudocode (MIN algorithm):**

```
function MIN_ASSIGN(F, C, m, n):
    assigned_rows ← ∅
    assigned_cols ← ∅
    mappings ← []
    while assigned_rows ≠ {1…m} AND assigned_cols ≠ {1…n}:
        best ← argmin_{i∉assigned_rows, j∉assigned_cols}
                    F(i,j)  s.t. isMinorSubgraph(C, i, j)
        if best = ∅: break
        mappings.append( MCS(R_{best.i}, P_{best.j}) )
        assigned_rows ← assigned_rows ∪ {best.i}
        assigned_cols ← assigned_cols ∪ {best.j}
    return mappings
```

where `isMinorSubgraph(C, i, j)` selects the pair with the smallest unique clique in its row or column — the most parsimonious assignment.

**RINGS algorithm:** Identical structure but prioritises pairs where the ring count is preserved: `|cycles(R_i)| = |cycles(P_j)|`, breaking ties via *E(i,j)* (bond energy). Ring-count parity is pre-computed once using `CycleFinder.vertexShort()` (CDK).

**MIXTURE algorithm:** Runs the first 5 assignment iterations with MIN-style selection (parsimony), then switches to MAX-style (coverage) for remaining unassigned pairs.

---

### 3.8 Cross-Algorithm Solution Ranking

After all algorithms complete, their candidate solutions are deduplicated by **mapping signature**:

    dedupeKey(φ) = sorted(bondChangePatterns(φ))

Solutions with identical bond-change patterns are considered equivalent; only the highest-priority-algorithm candidate is retained per unique key.

The surviving candidates are ranked by a **12-level comparator** (first difference wins):

| Priority | Criterion | Preference |
|----------|-----------|------------|
| 1 | Local score: `totalBondChanges + fragmentChanges` | Minimum |
| 2 | Total bond changes | Minimum |
| 3 | Fragment changes | Minimum |
| 4 | Bond dissociation energy sum | Minimum |
| 5 | Carbon bond changes | Minimum |
| 6 | Stereo changes | Minimum |
| 7 | Smallest-fragment atom count | Maximum |
| 8 | Graph similarity sum | Maximum |
| 9 | Energy score | Minimum |
| 10 | Fragment score | Minimum |
| 11 | Carbon score | Minimum |
| 12 | Algorithm priority (RINGS < MIN < MAX < MIXTURE) | Minimum |

**Early termination:** If any candidate has `totalBondChanges ≤ 2 AND fragmentChanges = 0`, it is accepted immediately.

---

### 3.9 Bond Change Annotation and Output

From the selected mapping φ, enumerate all bond changes in the ITS (Imaginary Transition State) graph:

    ITS(φ) = (A(R) ∪ A(P), E_form ∪ E_cleave ∪ E_order ∪ E_stereo)

where:
- **E_form:** bonds in *E(P)* absent in *E(R)* between φ-mapped atoms
- **E_cleave:** bonds in *E(R)* absent in *E(P)* between φ-mapped atoms
- **E_order:** bonds present on both sides but with changed multiplicity (e.g. C–C → C=C)
- **E_stereo:** stereocentres where R/S or E/Z configuration changes under φ

**Bond change fingerprint:** Each change is encoded as `ATOM1-ATOM2:WEIGHT` (e.g. `C-O:2`) and stored in four typed `IPatternFingerprinter` objects (formed/cleaved, order changes, stereo changes, reaction centre). The integer weight is the count of that pattern in the mapping.

**Reaction signature:** A canonical, sorted, hierarchical string:

    sig(φ) = "FC[" + sort(formed/cleaved) + "]|OC[" + sort(order) + "]|SC[" + sort(stereo) + "]|RC[" + sort(centre) + "]"

**Canonical hash:** SHA-256 of the concatenated sorted fingerprint strings, providing a permutation-invariant 64-character hex identifier for database indexing and exact-match deduplication.

---

## 4. Formal Properties

**Theorem 1 (Determinism).** For any fixed input reaction SMILES, RDT produces an identical mapping on every invocation. This follows from: (i) canonical SMILES is a unique normal form; (ii) the MCS cache returns identical results for identical keys; (iii) all tie-breaking criteria are total orders.

**Theorem 2 (Bond parsimony).** The selected mapping φ* satisfies:

    ∀ candidate φ ∈ Φ:  localScore(φ*) ≤ localScore(φ)

where `localScore = totalBondChanges + fragmentChanges`. This is a local optimum over the enumerated candidate set; the global optimum is not guaranteed (the problem is NP-hard).

**Complexity.** Let *n* = max molecule size (atoms). MCS computation is O(n^k) where *k* = clique size. In practice, the identity shortcut, size-ratio filter, and fingerprint filter together eliminate > 80% of pairs before MCS. The parallel phase runs at most 4 algorithm threads; the assignment step is O(m² × n²) per algorithm. Empirical throughput: 3–5 reactions/second on a 4-core laptop.

---

## 5. Benchmark Results

### 5.1 Golden Dataset

The Lin et al. (2022) golden dataset [3] contains 1,851 chemical reactions with expert-validated atom-atom mappings, spanning metabolic reactions, organic synthesis transformations, and ring opening/closing reactions. All published tools are evaluated on the **chemically-equivalent** metric: whether the mapping correctly identifies bond changes, regardless of atom-index labelling convention.

| Tool | Chem-Equiv | Mol-Map Exact | Training Data | Deterministic |
|------|-----------|---------------|---------------|---------------|
| **RDT v3.9.0** | **99.2%** | **~78%** | **None** | **Yes** |
| RXNMapper [4] | 83.74%† | — | Unsupervised | No |
| RDTool 2016 [1] | 76.18%† | — | None | Yes |
| ChemAxon | 70.45%† | — | Proprietary | Yes |

† Published figures from Lin et al. (2022).

### 5.2 Algorithm Selection Distribution (250-reaction slice)

| Algorithm selected | Count | % |
|--------------------|-------|---|
| RINGS | 229 | 91.6% |
| MIN | 16 | 6.4% |
| MAX | 5 | 2.0% |

RINGS resolves the majority of reactions via the funnel at a 2-4x computational saving over the full pipeline.

### 5.3 Performance

| Metric | Value |
|--------|-------|
| Mapping speed (laptop, 4-core) | 3–5 reactions/sec |
| Success rate | 100% (no unmapped reactions) |
| Test suite | 100% pass |

---

## 6. Implementation Notes

### 6.1 Dependencies

| Component | Version | Role |
|-----------|---------|------|
| SMSD | 6.10.1 | MCS engine: VF2++ subgraph isomorphism, circular/path fingerprints, MolGraph canonical SMILES (stereo-aware) |
| CDK | 2.12 | Molecule I/O, atom typing, aromaticity perception, ring finding |
| Java | 21+ | Platform |

### 6.2 Thread Safety

The mapping executor is a shared static `ExecutorService` (fixed thread pool, daemon threads). `MappingDiagnostics.REACTIONS` uses a `ConcurrentHashMap` with `remove()` on snapshot to prevent memory growth in batch processing. The MCS result cache is guarded by `ReadWriteLock`; the circular fingerprint cache uses `IdentityHashMap` per-thread (not shared).

### 6.3 Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| RINGS funnel threshold | 95% | Minimum atom coverage to accept RINGS alone |
| Reagent Tanimoto cutoff | 0.4 | Below this, molecule is candidate for reagent removal |
| Size ratio filter | 0.3 | Minimum atom-count ratio for MCS |
| FP similarity filter | 0.05 | Minimum path-FP Tanimoto for MCS |
| Substructure timeout | 5,000 ms | VF2++ hard timeout per pair |
| MCS cache capacity | 10,000 | LRU cache entries across reactions |
| Thread pool size | min(2, min(3, nCPU)) | Parallel mapping threads |

---

## 7. References

1. Rahman SA, Torrance G, Baldacci L, et al. "Reaction Decoder Tool (RDT): Extracting Features from Chemical Reactions." *Bioinformatics* 32(13):2065–2066, 2016. DOI: [10.1093/bioinformatics/btw096](https://doi.org/10.1093/bioinformatics/btw096)

2. Rahman SA, Cuesta S, Furnham N, et al. "EC-BLAST: a tool to automatically search and compare enzyme reactions." *Nature Methods* 11:171–174, 2014. DOI: [10.1038/nmeth.2803](https://doi.org/10.1038/nmeth.2803)

3. Lin A, Dyubankova N, Madzhidov TI, et al. "Atom-to-atom Mapping: A Benchmarking Study of Popular Mapping Algorithms and Consensus Strategies." *Molecular Informatics* 41(4):e2100138, 2022. DOI: [10.1002/minf.202100138](https://doi.org/10.1002/minf.202100138)

4. Schwaller P, Hoover B, Reymond J-L, et al. "Extraction of organic chemistry grammar from unsupervised learning of chemical reactions." *Science Advances* 7(15):eabe4166, 2021. DOI: [10.1126/sciadv.abe4166](https://doi.org/10.1126/sciadv.abe4166)

5. Luo YR. *Comprehensive Handbook of Chemical Bond Energies*. CRC Press, 2007.

6. Raymond JW, Willett P. "Maximum common subgraph isomorphism algorithms for the matching of chemical structures." *Journal of Computer-Aided Molecular Design* 16(7):521–533, 2002.

7. Ullmann JR. "An algorithm for subgraph isomorphism." *Journal of the ACM* 23(1):31–42, 1976.

---

*Reaction Decoder Tool is developed and maintained by BioInception PVT LTD.*
*Copyright (C) 2003–2026 Syed Asad Rahman. GNU LGPL v3.0.*
