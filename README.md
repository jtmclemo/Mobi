# Möbi v1.0.0

## Overview

Möbi is a MATLAB-based analytical workflow for transforming NetMHCpan-style peptide-HLA prediction outputs into sample-level mathematical representations, identifying cohort structure through topological grouping in weighted feature space, and selecting constrained neopeptide panels for the resulting cohorts.

The workflow supports multiple HLA-selection modes, including broad reference-panel filtering and patient-specific OptiType-based filtering. It explicitly preserves both the direct topological grouping and a downstream practical grouping for optimization.

Möbi is intended for exploratory neopeptide-panel design, cohort comparison, and reproducible computational analysis of peptide/HLA prediction outputs. It is **not** clinical decision software and should not be interpreted as a substitute for experimental validation or patient-specific clinical decision-making.

---

## Scientific and Mathematical Rationale

The scientific motivation for Möbi arises from the structure of peptide-HLA prediction data itself. NetMHCpan-style outputs are inherently multivariate, heterogeneous, and variable in length: each sample may contain many peptide-HLA prediction rows, and two samples may be biologically similar even when they do not share identical peptides row-for-row.

Similarity may instead appear in broader statistical and structural properties, such as peptide diversity, HLA breadth, binding-rank behavior, peptide promiscuity, peptide-family organization, and the overall distribution of prediction strengths across peptide-HLA combinations. Möbi therefore parses each sample into an interpretable feature vector that summarizes the peptide/HLA landscape of that sample, including quantities such as the number of unique peptides, the number of distinct HLA alleles represented, rank-based summary statistics, peptide promiscuity measures, HLA entropy, and peptide-family structure based on contiguous sequence overlap.

In this way, a variable-length biological prediction output is converted into a fixed-dimensional descriptor that can support geometry, topology, and downstream optimization.

Rather than imposing a single clustering threshold *a priori*, Möbi uses a persistent-homology approach to study how sample connectivity changes as the similarity scale is gradually relaxed. After feature standardization and weighting, each sample is treated as a point in a weighted Euclidean feature space. Pairwise distances between samples are then used to induce a Vietoris-Rips filtration, and the $H_0$ barcode records the merger of connected components as the filtration parameter increases.

In this application, the most immediate structural question is which samples become connected as one progressively relaxes the notion of proximity. $H_0$ persistence is therefore a natural object to study because it gives a multiscale description of cohort connectivity. Instead of defining groups only at one manually chosen threshold, Möbi first examines the merger structure across the filtration and then extracts a clustering threshold from that empirical topological behavior. This provides a more principled basis for grouping than relying entirely on an arbitrary cutoff chosen in advance.

Formally, let each sample $i$ be represented by a parser-derived feature vector $x_i \in \mathbb{R}^p$. Möbi standardizes each feature across samples, producing standardized values $z_{ik}$ for feature $k$ in sample $i$. It then computes a weighted Euclidean distance between samples $i$ and $j$:

$$
d_w(i,j)=\sqrt{\sum_{k=1}^{p} w_k^2 \left(z_{ik}-z_{jk}\right)^2}
$$

where $w_k$ is the configured weight assigned to feature $k$.

These pairwise distances define the order in which samples become connected as the filtration scale increases. Möbi uses this distance geometry to build a Vietoris-Rips filtration and track the $H_0$ merger structure; that is, the way connected components merge as the distance threshold grows. The resulting $H_0$ barcode provides a multiscale description of cohort connectivity. Rather than fixing one arbitrary clustering threshold in advance, Möbi derives a clustering scale from the empirical merger structure and reports the corresponding raw topological groups.

The representation step converts each variable-length NetMHCpan prediction table into a fixed-dimensional vector, the weighted distance defines the sample geometry, and the $H_0$ filtration captures how cohort structure changes across scales. This mapping can be interpreted as a feature embedding

$$
\Phi: \text{NetMHCpan outputs} \to \mathbb{R}^p
$$

which enables the use of geometric and topological methods. In that sense, the grouping stage is derived directly from the geometry induced by the parsed feature space and is not imposed independently of the data representation.

However, groups that are mathematically valid are not always convenient for downstream peptide optimization. Some groups may be too small to support practical cohort-level panel design, while others may be too large relative to the allowed peptide budget. Möbi therefore also produces a **practical grouping**, in which undersized groups may be merged and oversized groups may be split, while preserving an explicit audit trail showing how the practical grouping was derived from the raw topological result.

Möbi also distinguishes between patient-specific HLA information and population-informed approximation. When direct patient typing is available, OptiType-based filtering provides the most specific HLA support. When it is unavailable, Möbi can instead operate in global, region-based, or ancestry-group modes using internal reference logic derived from the Allele Frequency Net Database. This allows the workflow to remain self-contained and usable across multiple analysis settings while preserving a clear conceptual distinction between population-informed filtering and true individual HLA typing.

Another important design principle is explicit sensitivity analysis. Because Möbi combines feature engineering, weighted distance geometry, threshold-based grouping, and constrained peptide selection, modest parameter changes can propagate through the workflow in nontrivial ways. Möbi therefore includes optional robustness analysis to assess how stable the resulting cluster assignments and peptide outputs are under perturbations of clustering thresholds, peptide-count limits, and optimization weights. This makes the workflow more transparent and provides users with a direct way to evaluate whether a given result appears stable or highly parameter-sensitive.

---

## Disclaimer

Möbi is intended as an analytical and methodological framework for downstream study of peptide-HLA prediction outputs. It is designed to support exploratory panel design, cohort stratification, comparative computational analysis, and reproducible methodological experimentation. It is **not** a clinical pipeline, does not infer immunogenicity directly, and does not replace experimental validation, immunological assay data, or patient-specific clinical interpretation.

Several limitations should be kept in mind when interpreting results:

- The representation stage depends on the selected feature map, and any feature representation inevitably compresses biological complexity.
- The grouping stage depends on normalization, feature weights, and the threshold-selection rule used to extract clusters from the $H_0$ merger structure.
- The reference-panel HLA modes are only approximations when direct patient-specific typing is unavailable.
- NetMHCpan-derived quantities remain predictive model outputs rather than experimentally confirmed biological outcomes.
- The optimizer uses a transparent greedy heuristic rather than an exact combinatorial solver. While computationally efficient and interpretable, it does not guarantee global optimality.

---

## Software Architecture and File-Level Design

Möbi is implemented as a modular MATLAB suite because the underlying scientific workflow is itself modular. The pipeline must coordinate user input, parameter standardization, sample-level feature construction, topological grouping, constrained peptide selection, reference-panel support, optional robustness analysis, and baseline validation. These tasks are separated into dedicated files because they represent conceptually distinct layers of the method.

### `Mobi_frontend.m`

`Mobi_frontend.m` is the user-facing orchestration layer. It is the file that the user runs directly, and its responsibility is to coordinate the full workflow from input collection to final reporting. It handles HLA-source selection, OptiType file loading when patient-specific typing is used, NetMHCpan folder selection, filtering of candidate CSV files, iterative parsing of all accepted files, construction of the feature matrix, initiation of the TDA stage, reporting of merge events and group summaries, optional transition into optimization, optional robustness analysis, and run logging through a detailed text file.

### `Mobi_config.m`

`Mobi_config.m` defines the default state of the suite. It contains program metadata, expected NetMHCpan column names, parsing defaults, peptide-family overlap settings, TDA weights and threshold rules, optimization defaults and scoring weights, reference-panel defaults, logging/export conventions, and dormant future switches. If alternate non-user-prompted weights are desirable, they should be changed in this file.

### `Mobi_parsing.m`

`Mobi_parsing.m` is the representation layer. It reads a NetMHCpan-style CSV, verifies required columns, sorts EL-rank values, normalizes HLA nomenclature, optionally filters rows to a specified set of typed HLAs, and computes a sample-level feature representation. It organizes this information into a matrix `X`, with each row representing an input CSV file.

Möbi uses the following 13 parser-derived features for TDA:

- `num_unique_peptides`: number of distinct peptides in the sample
- `num_unique_HLAs_hit`: number of distinct HLA alleles represented after filtering
- `mean_rank`: mean `%Rank_EL` value across retained peptide-HLA rows
- `best_rank`: best (lowest) `%Rank_EL` value observed in the sample
- `rank_std`: standard deviation of `%Rank_EL` values
- `peptide_promiscuity_mean`: mean number of HLA alleles associated with each peptide
- `mean_best_rank_per_peptide`: mean of each peptide's best observed `%Rank_EL`
- `hla_entropy`: entropy of the HLA distribution across retained rows
- `rows_per_unique_peptide`: retained row count divided by unique peptide count
- `num_peptide_families`: number of peptide families defined by contiguous peptide overlap
- `largest_peptide_family_size`: size of the largest peptide-overlap family
- `frac_peptides_in_family_size_gt_1`: fraction of peptides belonging to multi-peptide families
- `mean_best_rank_per_family`: mean best rank across peptide families

The HLA entropy feature is computed using Shannon entropy:

$$
H_{\text{HLA}}=-\sum_h p_h \log p_h
$$

where $p_h$ is the empirical proportion of retained peptide-HLA rows assigned to HLA allele $h$. This quantity measures the breadth and balance of HLA representation within a sample.

### `Mobi_tda.m`

`Mobi_tda.m` implements the grouping stage. It performs conservative handling of non-finite values, z-score normalization, feature weighting, pairwise distance calculation, and weighted edge-list construction. It then builds an $H_0$ barcode by tracking the merger of connected components as the distance threshold grows, selects a threshold from the finite death values according to the configured rule, and extracts the raw cluster labels at that threshold.

When requested, it additionally constructs a practical optimization grouping by merging undersized groups and/or splitting oversized groups. The result therefore preserves both the direct topological grouping and a downstream optimization-oriented grouping, including an audit trail that records how raw clusters became practical groups.

### `Mobi_optimization.m`

`Mobi_optimization.m` implements constrained peptide-panel selection for a specified group of CSV files. It reads and concatenates the relevant NetMHCpan outputs, verifies required columns, normalizes HLA strings, filters rows by EL-rank threshold, constructs peptide families, computes peptide-level coverage and binding descriptors, and applies a greedy selection rule up to the requested panel size. The scoring function balances incremental patient coverage, incremental HLA coverage, transformed binding quality, prevalence across source files, family novelty, and a redundancy penalty based on contiguous peptide overlap. It returns the selected panel, the peptide feature table, and a stepwise selection history.

After samples are grouped, one must still decide which peptides best serve the cohort under a fixed size constraint while avoiding gratuitous redundancy. That is an optimization problem, not a topological one.

The optimization weights act on the following scoring features:

- `wPatientCoverage`: rewards peptides that newly cover previously uncovered samples or source files
- `wHLACoverage`: rewards peptides that newly cover previously uncovered HLA alleles
- `wBinding`: rewards peptides with stronger transformed EL-rank performance after EL-rank thresholding
- `wPrevalence`: rewards peptides that appear across more samples within the optimization group
- `wFamilyNovelty`: rewards peptides from peptide families not yet represented in the selected panel
- `wRedundancy`: penalizes peptides that overlap strongly with already selected peptides and are therefore likely to be redundant

Formally, at each greedy selection step, each candidate peptide $p$ is assigned a score

$$
S(p)=
w_1 \Delta C_{\text{patient}}(p)
+w_2 \Delta C_{\text{HLA}}(p)
+w_3 B(p)
+w_4 P(p)
+w_5 N(p)
-w_6 R(p)
$$

where:

- $\Delta C_{\text{patient}}(p)$ denotes the fraction of previously uncovered samples newly covered by peptide $p$
- $\Delta C_{\text{HLA}}(p)$ denotes the fraction of previously uncovered HLA alleles newly covered by peptide $p$
- $B(p)$ is transformed binding quality
- $P(p)$ is within-group prevalence
- $N(p)$ is peptide-family novelty
- $R(p)$ is redundancy penalty

The peptide with maximal score is selected at each iteration until the requested panel size is reached.

Binding quality is derived from EL-rank values by the transformation

$$
B(p) = \max\left(0, 1 - \frac{\mathrm{Rank}_{\mathrm{EL}}(p)}{\tau}\right)
$$

where $\tau$ is the configured EL-rank threshold. This converts lower, stronger EL-rank values into larger bounded scores.

In practice, the optimizer favors peptides that expand sample coverage, broaden HLA representation, retain favorable binding quality, and avoid repeatedly selecting near-duplicate peptide-family content. When weight normalization is enabled, the user-provided weights are rescaled so that their relative importance is preserved while the total weight mass remains standardized.

### `Mobi_run_optimization_suite.m`

`Mobi_run_optimization_suite.m` is the execution path for the peptide optimization workflow. It runs the global cohort optimization, honest/raw TDA-group optimization, and practical-group optimization as requested by the selected run mode, then returns a standardized suite structure with explicit peptide result fields.

These fields distinguish:

- `globalSelectedPeptides`
- `honestSelectedPeptidesByGroup`
- `practicalSelectedPeptidesByGroup`
- `honestSelectedPeptideUnion`
- `practicalSelectedPeptideUnion`
- `selectedPeptideUnionAcrossModes`

This file prevents the frontend and robustness analysis from developing separate interpretations of what optimization means. The frontend uses this suite for the main analysis, and robustness uses the same suite path for re-runs. That keeps peptide output semantics consistent across the main run, sensitivity checks, and future reporting paths.

### `Mobi_validate.m`

`Mobi_validate.m` provides shared invariant checks for stage-boundary outputs and inputs. It validates parser outputs, TDA inputs and outputs, optimization inputs and outputs, and the standardized optimization-suite output. These checks are intended to fail loudly if a future edit breaks expected structure, dimensionality, label consistency, coverage bounds, or peptide-union semantics.

### `Mobi_utils.m`

`Mobi_utils.m` contains helper routines shared across analytical modules, including HLA normalization, contiguous peptide-overlap checks, and connected-component extraction from adjacency matrices.

### `Mobi_reference.m`

`Mobi_reference.m` provides built-in population-informed HLA reference support. It returns an embedded AFND-derived table containing locus, allele, population, region, frequency, and ancestry-group fields, with persistent caching so the table can be reused efficiently across calls. Its purpose is to give Möbi a self-contained internal source for global, region-based, or ancestry-based HLA filtering when patient-specific typing is unavailable.

### `Mobi_robustness.m`

`Mobi_robustness.m` provides optional sensitivity analysis around the baseline run. It supports none, fast, standard, and full robustness modes. The fast mode performs a small threshold-focused check, the standard mode performs one-factor-at-a-time perturbations, and the full mode performs the exhaustive Cartesian sweep.

It builds baseline TDA and optimization configurations from the main configuration state, defines a sweep over clustering percentile, EL-rank threshold, maximum peptide count, and optimization weight sets, and compares each perturbed run to the baseline in terms of cluster-count stability, assignment stability, peptide-set identity, and peptide overlap. The output includes per-run records and an aggregate summary.

Peptide-set overlap between a perturbed run and the baseline is quantified using the Jaccard index:

$$
J(A,B)=\frac{|A \cap B|}{|A \cup B|}
$$

where $A$ is the baseline peptide set and $B$ is the perturbed peptide set. A value near 1 indicates strong stability of the selected panel under parameter perturbation.

Cluster-assignment stability is additionally quantified by

$$
\text{Agreement}=\frac{1}{n}\sum_{i=1}^{n}\mathbf{1}\{a_i^{(\text{run})}=a_i^{(\text{baseline})}\}
$$

where $a_i$ denotes the assigned group label for sample $i$. This measures the fraction of samples whose group assignments remain unchanged.

Mathematically, this robustness module should be interpreted as a structured local sensitivity analysis over parameter space rather than as a probabilistic uncertainty model. These robustness checks are deterministic re-runs under perturbed parameter settings rather than bootstrap or Monte Carlo uncertainty estimates.

### `Mobi_smoke_test.m`

`Mobi_smoke_test.m` is a compact, non-interactive validation routine. It identifies a small local example, parses two CSV files, assembles the feature matrix, runs TDA without plots, performs a small optimization call, and reports a concise success summary. Its purpose is to validate that the key components of the suite run successfully in the local environment.

### `Mobi_run_all_tests.m` and `tests/`

`Mobi_run_all_tests.m` runs the deterministic MATLAB regression tests stored in the `tests/` folder. These tests check parser feature regression, TDA geometry, deterministic peptide optimization, optimization-suite consistency, and robustness sweep sizing.

These tests complement the smoke test. The smoke test confirms that the pipeline can run on local data, while the regression tests check specific expected behavior in small synthetic cases.

---

## Input Requirements

Möbi expects a folder of NetMHCpan CSV files as its primary input.

### OptiType File Format

When using OptiType-based HLA filtering, the input file must obey the following structure:

- The first column must contain filenames matching the NetMHCpan CSV files exactly, including extensions.
- The next six columns must contain HLA alleles corresponding to that sample.

**Requirements**

- The filename column must match NetMHCpan CSV filenames exactly, including the `.csv` extension.
- Exactly six HLA alleles must be provided per sample (`2×A`, `2×B`, `2×C`).
- If only one allele is available for a given locus, enter the same allele in both columns for that locus.
- The first filename should begin in the first row and first column; do not add column headers.
- If OptiType-based HLA input is used, the selected file must be a `.csv`, `.xlsx`, or `.xls` file.

When OptiType mode is selected, filename matching must be exact between the OptiType table and the NetMHCpan CSV filenames. Möbi uses this filename match to determine which patient-specific HLA set should filter each NetMHCpan file. When reference-panel mode is selected instead, filtering is driven by the internal AFND-derived reference support.

**Example OptiType format**

| sample1.csv | A\*01:01 | A\*02:01 | B\*07:02 | B\*08:01 | C\*07:01 | C\*07:02 |  
|-------------|---------|---------|---------|---------|---------|---------|
| sample2.csv | A\*03:01 | A\*11:01 | B\*15:01 | B\*44:02 | C\*03:04 | C\*05:01 |

### NetMHCpan File Format

Each input CSV should correspond to one sample and must contain the required columns:

- `Identity`
- `Peptide`
- `MHC`
- `%Rank_EL`

Additional NetMHCpan output columns may be present, but Möbi relies on these required fields for parsing and downstream analysis.

If these formats are not satisfied, Möbi will skip or reject the corresponding file.

---

## MATLAB Requirements

Möbi is implemented and tested in MATLAB. The present development environment used MATLAB `R2025b`. The workflow uses MATLAB table and string functionality, along with numerical routines including pairwise distances, PCA, hierarchical linkage, and clustering.

In a standard MATLAB installation, these functions are commonly provided by MATLAB plus the **Statistics and Machine Learning Toolbox**.

---

## Expected Folder Layout

A typical local layout is:

```
Mobi/
  Mobi_frontend.m
  Mobi_config.m
  Mobi_parsing.m
  Mobi_tda.m
  Mobi_optimization.m
  Mobi_run_optimization_suite.m
  Mobi_validate.m
  Mobi_utils.m
  Mobi_reference.m
  Mobi_robustness.m
  Mobi_smoke_test.m
  Mobi_run_all_tests.m
  tests/

Final_CSVs/
  sample1.csv
  sample2.csv
  ...

Optional OptiType file:
  OptiType.xlsx
```
## Outputs

By default, each run writes a detailed text log named `Mobi_DetailedView_<timestamp>.txt` into the selected NetMHCpan folder. This log records the detailed run transcript and functions as the primary reproducibility artifact. Structured CSV and MAT exports are disabled by default in the current configuration, though the architecture leaves room for enabling them in future releases. The TDA stage also returns internal objects such as distances, edge lists, thresholds, raw cluster labels, practical optimization groups, and practical-grouping audit records, while the optimization stage returns selected peptides, peptide feature tables, and selection histories.

The detailed text log includes menu selections, merge events, cluster assignments, practical grouping audit information, normalized optimization weights, practical per-cluster selected peptides, and the robustness summary when robustness is run.

Furthermore, four plots are produced per run:
- a PCA projection of the parser-derived sample feature space. Note that This is a projection of 13-dimensional space onto 2-dimensional space, so any apparent groupings may not be entirely accurate
- an $H_0$ barcode plot summarizing connected-component persistence across filtration scale. *See also* the detailed text log.
- a finite $H_0$ death-value plot showing the distribution of merge thresholds and the selected clustering cutoff  
- a pairwise distance heatmap of the weighted Euclidean feature space 

## Robustness Modes

The optional robustness check supports four modes:

- `0`: None. Skips robustness checks entirely and reports only the main TDA and optimization run.
- `1`: Fast. Runs a minimal threshold-focused sensitivity check around the baseline clustering setting. This mode keeps optimization settings fixed and perturbs the TDA percentile around the baseline value. It is intended as a quick check of whether the sample grouping is stable to small changes in the clustering threshold, and is the most practical option when the user wants a rapid reassurance step.
- `2`: Standard. Runs a one-factor-at-a-time sensitivity analysis around the baseline configuration. In this mode, Möbi perturbs one parameter family at a time while holding the others fixed. The current implementation tests nearby TDA percentiles, alternative EL-rank thresholds, nearby maximum peptide limits, and several alternate optimization weight sets. This mode is intended to answer whether the final grouping and peptide outputs are reasonably stable under plausible local parameter changes without paying the cost of the full exhaustive sweep.
- `3`: Full. Runs the exhaustive cartesian robustness sweep. In this mode, Möbi evaluates all combinations across the selected percentile settings, EL-rank thresholds, peptide-count limits, and alternate optimization weight sets. This gives the most complete sensitivity picture, but it is also the slowest option and is best reserved for final methodological checking rather than quick exploratory use.

For each robustness run, Möbi compares the perturbed result to the baseline in terms of cluster-count stability, exact cluster-assignment agreement, exact peptide-set agreement, and peptide-set overlap. The robustness summary therefore helps distinguish whether a parameter change mainly affects sample grouping, peptide selection, or both.

## Peptide Output Semantics

Möbi distinguishes peptide outputs by optimization mode:

- `globalSelectedPeptides`: peptides selected from the full cohort-wide optimization.
- `honestSelectedPeptideUnion`: union of peptides selected across raw TDA-group optimizations.
- `practicalSelectedPeptideUnion`: union of peptides selected across practical-group optimizations.
- `selectedPeptideUnionAcrossModes`: union across all executed optimization modes.

## Quick Setup

1. Open MATLAB in the folder containing the full Möbi file suite.
2. Ensure that the required files are present in the same working directory:
   ```text
    Mobi_frontend.m
    Mobi_config.m
    Mobi_parsing.m
    Mobi_tda.m
    Mobi_optimization.m
    Mobi_run_optimization_suite.m
    Mobi_validate.m
    Mobi_utils.m
    Mobi_reference.m
    Mobi_robustness.m
    Mobi_smoke_test.m
    Mobi_run_all_tests.m
    tests/
   ```

3. Prepare a folder containing the NetMHCpan CSV files to be analyzed.
4. Confirm that each NetMHCpan CSV contains the required columns:
   ```
    Identity
    Peptide
    MHC
    %Rank_EL
   ```

5. If you plan to use patient-specific HLA typing, prepare an OptiType `.csv`, `.xlsx`, or `.xls` file with the aforementioned required formatting. 
6. In MATLAB, run:

   ```
   Mobi_frontend
   ```

7. Choose the HLA source mode when prompted:
   ```
    Global/common
    Region-based
    Ancestry-group
    OptiType file
   ```

8. Select the folder containing the NetMHCpan CSV files.
9. Follow the prompts for TDA grouping, peptide optimization, and optional robustness analysis.
10. After the run completes, review the generated .txt log in the selected NetMHCpan folder for a detailed file output.

## Validation

For a quick runtime check, run:

```
Mobi_smoke_test
```

For deterministic regression tests, run:

```
Mobi_run_all_tests
```

The smoke test checks that the core parser, TDA, and optimization pipeline can run on local example files. The regression tests check deterministic behavior for parser features, TDA geometry, optimization selection, optimization-suite consistency, and robustness sweep sizing.


## References

Gonzalez-Galarza FF, McCabe A, Santos EJM, Jones J, Takeshita L, Ortega-Rivera ND, Del Cid-Pavon GM, Ramsbottom K, Ghattaoraya G, Alfirevic A, Middleton D, Jones AR. Allele frequency net database (AFND) 2020 update: gold-standard data classification, open access genotype data and new query tools. Nucleic Acids Research. 2020;48(D1):D783-D788.

Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. Nucleic Acids Research. 2020;48(W1):W449-W454.

Jurtz V, Paul S, Andreatta M, Marcatili P, Peters B, Nielsen M. NetMHCpan-4.0: improved peptide-MHC class I interaction predictions integrating eluted ligand and peptide binding affinity data. Journal of Immunology. 2017;199(9):3360-3368.

Otter N, Porter MA, Tillmann U, Grindrod P, Harrington HA. A roadmap for the computation of persistent homology. EPJ Data Science. 2017;6:17.

## How to cite Möbi

If you use Möbi in research or downstream analyses, please cite it as software:

McLemore J. *Möbi: a Topological data analysis assisted optimization framework for designing broad-spectrum multivalent neopeptide pools* [software]. Version 1.0.0; 2026. Available from: GitHub repository.
