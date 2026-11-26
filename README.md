# `src/popgenml/preprocessing/` — detailed data QC

Five modules provide reproducible preprocessing for the 250-genotype,
5-location, 3-replication, 300-marker synthetic experiment and compatible
real-data tables.

## Modules
| Module | Main functionality |
|---|---|
| `genotype_qc.py` | Sample/marker call rates, MAF, heterozygosity, HWE diagnostic, genotype-mode imputation |
| `phenotype_qc.py` | Design completeness, duplicate detection, trait missingness, MAD outlier flags, means |
| `expression_qc.py` | Integer-count checks, sample/gene filtering, library sizes, log2(CPM+1) |
| `metadata_qc.py` | Genotype/location/replication validation and cross-table ID checks |
| `data_integration.py` | Merge genotype markers and phenotype means without simulation-truth leakage |

## Dependencies
`pip install numpy pandas`

From repository root: `export PYTHONPATH="$PWD/src:$PYTHONPATH"` (Linux/WSL).

## Usage
```python
from popgenml.preprocessing import (
    GenotypeQC, PhenotypeQC, ExpressionQC, MetadataQC, DataIntegrator
)

filtered_geno, marker_qc, sample_qc, stats = GenotypeQC(
    min_marker_call_rate=0.90, min_sample_call_rate=0.90,
    min_maf=0.01
).filter(genotype_df, output_dir="results/qc/genotype")

flagged, trait_qc, means = PhenotypeQC(
    expected_locations=5, expected_reps=3
).inspect(phenotype_df, traits=["boll_weight_g", "fiber_strength_g_tex"],
          output_dir="results/qc/phenotype")

design_summary = MetadataQC(250,5,3).validate_design(
    phenotype_df, output_dir="results/qc/metadata")

ml_table = DataIntegrator.join_genotype_phenotype(
    filtered_geno, means, trait="fiber_strength_g_tex",
    output_dir="results/qc/integration")

# For actual or simulated raw integer counts only:
filtered_counts, log_cpm, sample_info, gene_info = ExpressionQC(
    min_total_count=10, min_samples_expressed=3
).process(expression_counts_df, output_dir="results/qc/expression")
```

## Input requirements
- Genotypes: `genotype_id` plus `SNP_...` columns coded 0/1/2/NaN.
- Phenotypes: one row per genotype-location-replication plot, with trait
  columns; duplicate plot keys are not permitted.
- Expression: one row per sample, `sample_id`, optional `group`, gene count
  columns; counts must be nonnegative integers.
- Trait means: `genotype_id`, `location_id` and numeric trait columns.

## Scientific and analytical limitations
- **Cotton is allotetraploid.** The diploid 0/1/2 QC routines are only
  appropriate for a validated diploid-like marker encoding, not arbitrary
  polyploid dosage calls.
- HWE chi-square diagnostics can be misleading in structured, selected,
  inbred, or small populations; they are reported, never auto-filtered.
- Mode imputation is a simple baseline. More advanced imputation methods
  should be evaluated against missingness and population structure.
- Outliers are **flagged, not deleted**. Transgressive segregants and genuine
  biological extremes must not be removed without review.
- CPM is an exploratory normalization, not a DESeq2 normalization substitute.
  Use untransformed raw counts with DESeq2 or edgeR for differential expression.
- When genotype-location rows repeat genotypes, group cross-validation by
  genotype (or use environment-aware validation) to prevent information leakage.
- Synthetic truth labels must never be included as ML features.
