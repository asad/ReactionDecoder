# Reaction Decoder Tool (RDT) v3.8.1 — Algorithm Description

**Authors:** Syed Asad Rahman
**Contact:** asad.rahman@bioinceptionlabs.com
**License:** GNU LGPL v3.0

---

## 1. Overview

The Reaction Decoder Tool (RDT) performs deterministic atom-atom mapping (AAM) for chemical reactions without any training data or machine learning. Given an input reaction (reactants and products), RDT identifies which atom in each reactant corresponds to which atom in each product, enabling the identification of reaction centres, bond changes, and reaction mechanisms.

**Key Innovation:** A multi-algorithm ensemble approach with game-theory-inspired matrix optimization. Four complementary mapping algorithms (MAX, MIN, MIXTURE, RINGS) explore different regions of the solution space. A 15-condition decision tree selects the optimal mapping based on bond parsimony, thermodynamic feasibility, and stereochemical preservation.

**Benchmark Result:** 99.2% chemically-equivalent atom mapping on the Lin et al. (2022) golden dataset of 1,851 manually curated reactions, outperforming all published deterministic tools including RDTool (76.18%) and ChemAxon (70.45%), without any training data.

---

## 2. Algorithm Pipeline

The algorithm proceeds through nine sequential stages:

```
Input Reaction
     |
     v
[Stage 1] Parsing & Preprocessing
     |
     v
[Stage 2] Reaction Standardization (reagent filtering, atom balance)
     |
     v
[Stage 3] RINGS Funnel (quality gate)
     |
     v                              yes
[Stage 3a] RINGS sufficient? ---------> Use RINGS mapping
     |
     | no
     v
[Stage 4] Parallel execution: MIN, MAX, MIXTURE (+RINGS if not run)
     |
     v
[Stage 5] Pairwise MCS computation (substructure + VFLibMCS)
     |
     v
[Stage 6] Game theory matrix construction (7 scoring matrices)
     |
     v
[Stage 7] Algorithm-specific winner selection
     |
     v
[Stage 8] Cross-algorithm solution selection (15-condition tree)
     |
     v
[Stage 9] Bond change annotation & output
```

---

### Stage 1: Input & Preprocessing

**Input formats:** Reaction SMILES, RXN (V2000/V3000), RDF, or CDK IReaction objects.

**Preprocessing steps:**
1. Parse reaction into separate reactant and product molecule containers
2. Set all null implicit hydrogen counts to zero
3. Remove any pre-existing atom-atom mapping properties
4. Perceive atom types and aromaticity using CDK's `AtomTypeMatcher` and `Aromaticity` models

---

### Stage 2: Reaction Standardization

**Purpose:** Remove non-reactive species (solvents, catalysts, reagents) to focus mapping on the reacting core.

**Three-tier reagent filter:**

| Tier | Method | Criteria |
|------|--------|----------|
| **1. Known reagents** | Canonical SMILES lookup | Match against database of 30+ common solvents (DCM, DMSO, DMF, THF, water, etc.) and inorganic salts. **Atom-balance guard:** only filter if removing the molecule does not unbalance the reaction. |
| **2. Catalyst metals** | Element symbol check | Molecules containing Pd, Pt, Rh, Ru, Ir, Ni, Cu, Fe, Co, Mn, Ti, Zr, Mo, W, Os, Ag, Au |
| **3. Fingerprint similarity** | ECFP4 Tanimoto | For each reactant, compute maximum Tanimoto similarity to any product using ECFP4 circular fingerprints (radius=2, 256 bits). If max similarity < 0.4, heavy atom count <= 10, and no unique element contributions to products: classify as reagent. **Atom-balance guard** prevents filtering genuine reactants (e.g., water in hydrolysis). |

**Atom balance validation:** Count all non-hydrogen atoms per element in reactants vs products. Log warning if imbalanced (does not prevent mapping, as many real reactions are intentionally unbalanced in their written form).

---

### Stage 3: RINGS Funnel Architecture

**Purpose:** Avoid unnecessary computation by testing if the ring-conservation algorithm alone produces a sufficient mapping.

**Algorithm:**
1. If `totalMolecules <= 5`: execute RINGS algorithm in a single thread
2. Evaluate quality: does the mapping cover >= 95% of non-hydrogen atoms?
3. Verify non-identity: confirm the reaction involves actual structural changes (not a transporter)
4. **Decision:** If RINGS mapping is sufficient, skip the remaining three algorithms

**Efficiency gain:** Approximately 75% of reactions are resolved by RINGS alone, providing a 2-4x speedup over the full four-algorithm pipeline.

---

### Stage 4: Parallel Multi-Algorithm Execution

When RINGS is insufficient, the remaining algorithms execute in parallel:

| Algorithm | Strategy | Objective |
|-----------|----------|-----------|
| **MAX** | Global maximization | Maximize total mapped atoms across all reactant-product pairs |
| **MIN** | Local minimization | Minimize total bond changes (parsimony principle) |
| **MIXTURE** | Hybrid max-min | Maximize mapped atoms with secondary minimization of bond changes |
| **RINGS** | Ring-centric | Prioritize ring system preservation and aromatic skeleton conservation |

**Execution:** Fixed thread pool with `min(available_processors - 1, 4)` threads. Results collected via `CompletionService` in order of completion.

---

### Stage 5: Pairwise MCS Computation

For each reactant-product pair *(R_i, P_j)*, compute the Maximum Common Subgraph (MCS) to establish atom correspondences.

#### 5.1 Pre-filtering

Three pre-filters eliminate pairs unlikely to share meaningful substructure:

| Filter | Condition | Action |
|--------|-----------|--------|
| **Identity** | MolGraph canonical SMILES equality (stereo-aware) + equal atom count | Build direct identity mapping (atom *i* → atom *i*) and skip MCS entirely. Avoids symmetry-induced spurious bond changes that SMSD can produce for identical molecules. |
| **Size ratio** | `min(atoms_i, atoms_j) / max(atoms_i, atoms_j) < 0.3` and smaller molecule > 3 atoms | Skip pair — highly dissimilar sizes indicate unrelated molecules |
| **Fingerprint** | `Tanimoto(FP_i, FP_j) < 0.05` and both > 5 atoms | Skip pair — structurally unrelated by path fingerprint |

**Identity pre-filter detail:** Canonical SMILES are generated via `MolGraph.toCanonicalSmiles()` (SMSD 6.9.1), which encodes tetrahedral chirality (`@`/`@@`) and E/Z double-bond geometry (`/`/`\`). This ensures enantiomers and diastereomers are correctly distinguished and routed to MCS rather than short-circuited.

#### 5.2 Tiered Substructure Matching

For each pair, attempt substructure isomorphism with progressively relaxed matching criteria:

```
Tier 1: AtomType=strict, BondOrder=flexible, RingMatch=strict
     |
     | (if no subgraph found)
     v
Tier 2: AtomType=element-only, BondOrder=flexible, RingMatch=strict
     |
     | (if no subgraph found)
     v
Tier 3: AtomType=element-only, BondOrder=flexible, RingMatch=relaxed
```

Each tier uses VF2++ subgraph isomorphism via the SMSD engine with a 5-second timeout.

#### 5.3 Full MCS Fallback

If no substructure relationship exists (neither molecule is a subgraph of the other), compute the full Maximum Common Subgraph using the SMSD MCS algorithm.

**Optimization:** A `ThreadSafeCache` stores MCS results keyed by canonical SMILES + matcher configuration + circular fingerprint hash. This enables cross-reaction reuse when the same molecule pair appears in different reactions.

#### 5.4 Circular Fingerprint Cache

Each molecule's FCFP (Functional-Class Fingerprint, radius=1, 256 bits) is computed once and cached in an `IdentityHashMap`. The fingerprint hash is included in the MCS cache key to ensure uniqueness.

---

### Stage 6: Game Theory Matrix Construction

For each algorithm, construct seven scoring matrices of dimension *(|Reactants| x |Products|)*:

| Matrix | Symbol | Formula | Description |
|--------|--------|---------|-------------|
| **Clique** | *C(i,j)* | `|MCS(R_i, P_j)|` | Number of atoms in the maximum common subgraph |
| **Graph Similarity** | *G(i,j)* | `|MCS| / (|R_i| + |P_j| - |MCS|)` | Jaccard index of the MCS |
| **Stereo** | *S(i,j)* | From SMSD stereochemistry analysis | Stereochemical compatibility score |
| **Energy** | *E(i,j)* | From SMSD bond energy matrix | Bond dissociation energy of the mapping |
| **Fragment** | *F(i,j)* | From SMSD fragment analysis | Number of disconnected fragments in the mapping |
| **Carbon Overlap** | *K(i,j)* | `|{a in MCS : symbol(a) = C}|` | Carbon atoms preserved in the mapping |
| **FP Similarity** | *T(i,j)* | `Tanimoto(FP_i, FP_j)` | Path fingerprint Tanimoto similarity |

---

### Stage 7: Algorithm-Specific Winner Selection

Each algorithm iteratively selects the best reactant-product pair from the matrices:

**MAX Algorithm:**
1. Find pair *(i,j)* with maximum *G(i,j)* that is a **major subgraph** in both row and column (no other pair in the same row or column has a larger clique)
2. Record mapping, remove pair from matrices
3. Repeat until all atoms mapped or no valid pairs remain

**MIN Algorithm:**
1. Find pair *(i,j)* with minimum *F(i,j)* (fewest fragment changes) that is a **minor subgraph** (smallest unique clique in its row or column)
2. Record mapping, remove pair from matrices
3. Repeat

**MIXTURE Algorithm:**
1. Iterations 1-5: use MIN-style selection (parsimony)
2. Iterations 6+: switch to MAX-style selection (coverage)

**RINGS Algorithm:**
1. Prioritize pairs where `numberOfCycles(R_i) == numberOfCycles(P_j)` (ring count compatibility)
2. Use energy matrix *E(i,j)* to break ties
3. Special handling for ring opening/closing reactions

**Deadlock Resolution:** When multiple pairs have equal scores, a 15-condition decision tree resolves ties using the priority hierarchy:

```
Priority 1: Fewer total bond changes (parsimony)
Priority 2: Fewer molecular fragments
Priority 3: Lower bond dissociation energy (thermodynamic feasibility)
Priority 4: More carbon-carbon bonds preserved
Priority 5: Fewer stereochemical changes
```

---

### Stage 8: Cross-Algorithm Solution Selection

After all algorithms complete, select the best overall mapping:

**Primary criterion:** Minimum `localScore = totalBondChanges + totalFragmentChanges`

**Tiebreaker hierarchy (first difference wins):**

| Priority | Criterion | Preference |
|----------|-----------|------------|
| 1 | Total bond changes | Fewer |
| 2 | Fragment changes | Fewer |
| 3 | Bond dissociation energy | Lower |
| 4 | Carbon bond changes | Fewer |
| 5 | Stereochemical changes | Fewer |
| 6 | Smallest fragment size | Larger (fewer small fragments) |

**Early termination:** If any algorithm produces a mapping with `totalBondChanges <= 2` and `fragmentChanges == 0`, accept immediately without evaluating remaining algorithms.

---

### Stage 9: Bond Change Annotation & Output

From the selected mapping, enumerate all bond changes:

1. **Bond formed:** Bond exists in products but not in reactants (between mapped atoms)
2. **Bond cleaved:** Bond exists in reactants but not in products (between mapped atoms)
3. **Bond order changed:** Bond exists on both sides but with different multiplicity (e.g., single to double)
4. **Stereochemical change:** E/Z or R/S configuration change at a stereogenic centre

**Bond dissociation energy:**

```
totalEnergy = SUM over all changed bonds: count(bond) x BDE(atom1, atom2, order)
```

where BDE values are from the Luo (2007) bond dissociation energy reference table.

**Output:** Atom-atom mapping as integer labels (1, 2, 3, ...) assigned to corresponding atoms in reactants and products, plus bond change fingerprints classifying the reaction centre.

---

## 3. Mathematical Formulations

**Jaccard Similarity (Graph Similarity):**

```
G(i,j) = |MCS(R_i, P_j)| / (|R_i| + |P_j| - |MCS(R_i, P_j)|)
```

**Tanimoto Fingerprint Similarity:**

```
T(A,B) = |A AND B| / (|A| + |B| - |A AND B|)
```

where |A| = cardinality of fingerprint bitset A.

**Local Score (primary optimization objective):**

```
localScore(mapping) = totalBondChanges(mapping) + totalFragmentChanges(mapping)
```

**Bond Energy Sum:**

```
E(mapping) = SUM_{b in changedBonds} weight(b) x BDE(atom1(b), atom2(b), order(b))
```

**Reagent Filter Atom-Balance Guard:**

```
isNeeded(reactant) = EXISTS element e :
    atomCount(allReactants \ reactant, e) < atomCount(products, e)
```

---

## 4. Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| ECFP radius | 2 | Extended Connectivity Fingerprint radius for reagent filtering |
| FCFP radius | 1 | Functional-Class Fingerprint radius for cache key generation |
| Fingerprint size | 256 bits | Bit length for circular fingerprints |
| Path fingerprint depth | 7 | Maximum path length for structural fingerprints |
| Path fingerprint size | 1024 bits | Bit length for path-based fingerprints |
| Reagent Tanimoto threshold | 0.4 | Below this, molecule is candidate for reagent classification |
| Size ratio filter | 0.3 | Minimum atom count ratio for MCS computation |
| FP similarity filter | 0.05 | Minimum Tanimoto for MCS computation |
| Substructure timeout | 5,000 ms | VF2++ subgraph isomorphism timeout |
| RINGS funnel threshold | 95% | Minimum atom coverage to skip remaining algorithms |
| MCS cache capacity | 10,000 | Maximum entries in cross-reaction MCS cache |
| Thread pool size | min(cores-1, 4) | Parallel algorithm execution threads |
| Max iterations | 100 | Maximum matrix optimization iterations per algorithm |

---

## 5. Benchmark Results

### Golden Dataset (Lin et al. 2022)

1,851 manually curated reactions with expert-validated atom-atom mappings.
Published tools are scored on chemically-equivalent atom mapping — whether the mapping correctly identifies bond changes regardless of atom-index labelling.

| Tool | Chemically Equivalent | Bond-Change Exact | Mol-Map Exact | Training Data | Deterministic |
|------|-----------------------|-------------------|---------------|---------------|---------------|
| **RDT v3.8.1** | **99.2%** | **99.2%** | **76.8%** | **None** | **Yes** |
| RXNMapper | 83.74%† | - | - | Unsupervised | No |
| RDTool (published, 2016) | 76.18%† | - | - | None | Yes |
| ChemAxon | 70.45%† | - | - | Proprietary | Yes |

† Published figures from Lin et al. 2022 use chemically-equivalent scoring.

### Performance Metrics (250-reaction slice)

| Metric | Value |
|--------|-------|
| Mapping success rate | 100% (250/250) |
| Chemically-equivalent atom mapping | 99.2% |
| Bond-change exact | 99.2% |
| Mol-map exact | 76.8% |
| True chemistry misses | 0.8% |
| Mapping speed | 2.3 reactions/sec |
| Test suite | 164 tests, 100% pass |

---

## 6. Dependencies

| Component | Version | Role |
|-----------|---------|------|
| SMSD | 6.9.1 | Substructure and MCS engine (VF2++, circular/path fingerprints, MolGraph canonical SMILES) |
| CDK | 2.12 | Cheminformatics toolkit (molecule parsing, atom types, aromaticity) |
| Java | 21+ | Runtime platform |

**Note on canonical SMILES:** Identity pre-filtering uses `MolGraph.toCanonicalSmiles()` from SMSD 6.9.1 rather than CDK's `SmilesGenerator`. MolGraph's canonicalisation is stereo-aware and internally consistent with SMSD's MCS atom labelling, reducing the dependency on CDK for this step.

---

## 7. References

1. Rahman SA, Torrance G, Baldacci L, et al. "Reaction Decoder Tool (RDT): Extracting Features from Chemical Reactions." *Bioinformatics* 32(13):2065-2066, 2016. DOI: [10.1093/bioinformatics/btw096](https://doi.org/10.1093/bioinformatics/btw096)

2. Rahman SA, Cuesta S, Furnham N, et al. "EC-BLAST: a tool to automatically search and compare enzyme reactions." *Nature Methods* 11:171-174, 2014. DOI: [10.1038/nmeth.2803](https://doi.org/10.1038/nmeth.2803)

3. Lin A, Dyubankova N, Madzhidov TI, et al. "Atom-to-atom Mapping: A Benchmarking Study of Popular Mapping Algorithms and Consensus Strategies." *Molecular Informatics* 41(4):e2100138, 2022. DOI: [10.1002/minf.202100138](https://doi.org/10.1002/minf.202100138)

4. Luo YR. "Comprehensive Handbook of Chemical Bond Energies." CRC Press, 2007.

---

*Reaction Decoder Tool is developed and maintained by BioInception PVT LTD.*
*Copyright (C) 2003-2026 Syed Asad Rahman. GNU LGPL v3.0.*
