# Golden Benchmark Report

Release: RDT v3.8.0

Date: 2026-04-02

Dataset: Lin et al. 2022 golden dataset

Commands:

```bash
mvn -q -Dtest=GoldenDatasetBenchmarkTest -Dgolden.max=20 test
mvn -q -Dtest=GoldenDatasetBenchmarkTest -Dgolden.max=100 test
mvn -q -Dtest=GoldenDatasetBenchmarkTest -Dgolden.max=250 test
```

## Metric definitions

- `Mapping success`: mapper returned a solution without hard failure
- `Mol-map exact`: exact equality of induced reactant-molecule to product-molecule relation set
- `Atom-map exact`: exact reactant-atom to product-atom match against the gold file
- `Atom-map chemically equivalent`: same bond-change set as gold, even if atom numbering differs
- `Bond-change exact`: exact equality of the full bond-change set
- `Bond-change count exact`: exact equality of total bond-change count
- `Bond-change type exact`: exact equality of `FORM`/`BREAK`/`ORDER` counts
- `Reaction-center exact`: exact equality of the changed-atom set
- `Reaction-center atoms`: atom-level reaction-center accuracy
- `True chemistry miss`: bond-change set differs from gold
- `Speed`: reactions per second for the measured slice

## Results

| Slice | Mapping success | Mol-map exact | Atom-map exact | Atom-map chemically equivalent | Bond-change exact | Bond-change count exact | Bond-change type exact | Reaction-center exact | Reaction-center atoms | True chemistry miss | Speed |
|------|-----------------|---------------|----------------|-------------------------------|-------------------|-------------------------|------------------------|-----------------------|----------------------|---------------------|-------|
| `20` | `20/20 (100.0%)` | `11/20 (55.0%)` | `2/20 (10.0%)` | `20/20 (100.0%)` | `20/20 (100.0%)` | `20/20 (100.0%)` | `20/20 (100.0%)` | `20/20 (100.0%)` | `870/870 (100.0%)` | `0/20 (0.0%)` | `1.7 rxn/sec` |
| `100` | `100/100 (100.0%)` | `71/100 (71.0%)` | `27/100 (27.0%)` | `100/100 (100.0%)` | `100/100 (100.0%)` | `100/100 (100.0%)` | `100/100 (100.0%)` | `100/100 (100.0%)` | `4509/4509 (100.0%)` | `0/100 (0.0%)` | `2.6 rxn/sec` |
| `250` | `250/250 (100.0%)` | `189/250 (75.6%)` | `58/250 (23.2%)` | `248/250 (99.2%)` | `248/250 (99.2%)` | `248/250 (99.2%)` | `248/250 (99.2%)` | `248/250 (99.2%)` | `11747/11769 (99.8%)` | `2/250 (0.8%)` | `2.0 rxn/sec` |

## Interpretation

- The current branch is strong on chemistry correctness.
- The main benchmark penalty is strict atom numbering, not wrong reaction chemistry.
- On the first `100` reactions there were `0` true chemistry misses.
- On the `250` slice there were `2` true chemistry misses and `190` alternate valid maps.
- Reaction-center quality is effectively saturated on the measured slices.
- Mol-map exact is much higher than atom-map exact, which is consistent with symmetry-equivalent atom labeling inside otherwise correct component mappings.

## Practical conclusion

The current benchmark should be read as:

- high chemistry correctness
- moderate molecule-level exactness
- low strict atom-number exactness
- low throughput relative to the long-term target

The next optimization target should be strict atom-map canonicalization under symmetry, not bond-change chemistry.
