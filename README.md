## Context Predictions — Variant Margin Analysis

This repo generates codon-variant sequences around a mutation, runs Evo inference to produce logits, and creates figures analyzing entropy and log-likelihood across different left-margin sizes.

### Directory structure

Key folders and files:

```
context_predictions/
  input/
    codon_table                 # TSV: codon\taa
    human_contigs_src.fasta     # FASTA: contigs used as sources
    updated_data.tsv            # Master table of variants to run

  scripts/
    analyze_codon_pos_marg.sh   # Orchestrates per-position run (generate → evo_gcp → plots)
    gen_for_cod_var_marg.py     # Generates FASTA with all codon variants at a pos, info + query tables
    plot_margins.R              # Builds figures from Evo logits for selected variants and margins
    utils.r                     # R helpers for loading logits, computing metrics, and plotting

  run_mutation_analysis.sh      # Batch runner over rows in input/updated_data.tsv

  jobs/                         # Downloaded Evo outputs per job/job_version
  output/                       # Generated FASTA + tables per AID (created by scripts)
  figures/                      # All generated plots organized by AID and margin
  requirements.txt              # Python deps (for gen_for_cod_var_marg.py)
```

### What the scripts do

- **`scripts/gen_for_cod_var_marg.py`**: Given a gene interval and an amino-acid coordinate:
  - Reads `input/codon_table` and the target contig from `input/human_contigs_src.fasta`.
  - Generates a FASTA with all 64 codon substitutions at the target position across a set of left margins `[0, 10, 100, 1000, 10000, 100000]` (bounded by gene start).
  - Writes a sequence info table (`.tab`) with the original codon/AA, target codon/AA, synonymous codon, gene bounds, AID, and left margin.
  - Writes a query table (`_query.tab`) mapping each generated sequence ID to `start`/`end` coordinates for Evo.

- **`scripts/analyze_codon_pos_marg.sh`**:
  - Calls `gen_for_cod_var_marg.py` to produce: `output/<aid>/margins_<seq_id>_<pos>.fasta`, `.tab`, and `_query.tab`.
  - Submits an Evo job via `evo_gcp submit` with `--output_type logits` and `--wait`.
  - Downloads results into `jobs/<aid>-<job_version>/output` via `evo_gcp download`.
  - Calls `scripts/plot_margins.R` to generate figures for the variants and margins.

- **`scripts/plot_margins.R`**:
  - Loads logits from Evo output files named `input_<variant_name>_logits.npy`, where `variant_name` is `<aaCoord>_<aa>_<codon>_<margin>`.
  - Computes per-position entropy and log-likelihood against the no-margin sequence (`_0`).
  - Produces per-margin “full” plots and stacked/comparison summaries under `figures/<aid>/`:
    - `margin_<size>/full/*.pdf` (entropy, log-likelihood per variant)
    - `stacked/*.pdf` (stacked by margin)
    - `compare/*.pdf` (line comparisons across margins)
    - `total_loglik/total_loglik_all_variants.pdf`

- **`scripts/utils.r`**: R helpers for loading numpy logits, computing probabilities/entropy, building plot data, and composing figures.

- **`run_mutation_analysis.sh`**:
  - Batch runner over `input/updated_data.tsv`.
  - Accepts optional filters: `AID` and `POSITION`.
  - For each row (respecting filters), invokes `scripts/analyze_codon_pos_marg.sh` with the appropriate arguments.

### End-to-end workflow to create figures

1) Install dependencies

- Python packages (used by `gen_for_cod_var_marg.py`):

```bash
pip install -r requirements.txt
```

- R packages (used by `plot_margins.R`):

```bash
Rscript -e 'install.packages(c("optparse","Biostrings","ggplot2","patchwork"), repos="https://cloud.r-project.org")'
```

- Evo CLI: Ensure `evo_gcp` is installed/authenticated and that you have access to submit and download jobs.

2) Prepare inputs

- `input/human_contigs_src.fasta`: FASTA containing the contig(s) by `seq_id` referenced in `updated_data.tsv`.
- `input/codon_table`: Tab-separated with headers `codon` and `aa`.
- `input/updated_data.tsv`: Tab-separated table with columns including: `aid`, `gene`, `contig`, `start`, `end`, `strand`, `flipped`, `src_codon`, `tgt_codon`, `mut_pos`.

3) Run per-position (single example)

```bash
scripts/analyze_codon_pos_marg.sh \
  <AA_POS> <GENE_START> <GENE_END> <SEQ_ID> <AID> \
  input/human_contigs_src.fasta <TARGET_CODON> [left_margin=2000] [right_margin=1000]
```

This will:
- Generate `output/<aid>/margins_<seq>_<pos>.fasta`, `.tab`, and `_query.tab`.
- Submit to Evo with the FASTA and query table, wait for completion, download results to `jobs/`.
- Create figures under `figures/<aid>/...`.

4) Run batch over `updated_data.tsv`

```bash
./run_mutation_analysis.sh              # all AIDs and positions
./run_mutation_analysis.sh BAA          # only for AID BAA
./run_mutation_analysis.sh BAA 154      # only position 154 for AID BAA
```

### Evo outputs and file naming

- Evo outputs used by the plotting step are numpy arrays named like:
  - `jobs/<aid>-<job_version>/output/input_<aaCoord>_<aa>_<codon>_<margin>_logits.npy`
- The R code expects exactly this naming and uses the no-margin (`_0`) sequence as the baseline for log-likelihood.

### Notes on margins and query table

- Margins generated are `[0, 10, 100, 1000, 10000, 100000]` filtered to be less than `gene_start`.
- The query table (`*_query.tab`) maps each generated variant ID to `start`/`end` coordinates; this is what enables Evo to score the gene region with and without additional upstream context (margin).

### Troubleshooting

- If downloaded outputs contain `.gstmp` files in `jobs/.../output`, the download may be incomplete; re-run the download step. (.gstmp means partial download)
- Ensure `seq_id` in the FASTA matches the target contig ID and that `start`/`end` bounds are correct.
- Confirm `input/codon_table` has header `codon\taa` and includes all 64 codons including stops (`*`).


