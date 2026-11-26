# GATK4 germline variant-calling pipeline + BAM→VCF web app

A reproducible Snakemake pipeline implementing the [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)
germline short-variant workflow — multi-lane FASTQ in, a filtered VCF and
a QC report out — plus a small Flask web app that runs the same variant
calling core on a single BAM someone hands it, for when the full
pipeline is more than they need.

Demonstrated end-to-end on a small public multi-lane human dataset, so
every number in this README is a real result from a real run, not a
description of what the pipeline is supposed to do.

## Why this project

Built to cover a gap that shows up constantly in bioinformatics job
specs: NGS variant calling with the GATK toolchain. It sits alongside
two other projects — a [ChIP-seq peak-calling pipeline](https://github.com/amaaniu/chipseq-peak-calling-pipeline)
and my [MSc dissertation](https://github.com/amaaniu/BIO702P-Research-Project)
on long-read scRNA-seq — as the third leg: NGS/short-read variant
calling, ChIP-seq/epigenomics, and single-cell long-read, covering the
main data types that come up in cancer genomics and computational
biology roles.

## Pipeline overview

```
Per-lane FASTQ (5 lanes, one sample)
    │
    ├─▶ fastp                               (adapter/quality trimming + QC)
    │
    ▼
BWA-MEM alignment (per-lane read groups) ──▶ sort/index
    │
    ▼
Merge lanes ──▶ MarkDuplicates
    │
    ▼
BaseRecalibrator + ApplyBQSR   (known sites: dbSNP, Mills indels)
    │
    ├─▶ samtools flagstat                   (alignment QC)
    │
    ▼
HaplotypeCaller  ──▶  raw VCF
    │
    ▼
Split SNPs/indels ──▶ hard-filter each (GATK best-practices thresholds) ──▶ merge
    │
    ▼
bcftools stats + pipeline summary  ──▶  MultiQC report
```

Every step is a Snakemake rule in [`workflow/Snakefile`](workflow/Snakefile).
Reads are aligned and recalibrated per-lane before merging, because BQSR
and duplicate marking need correct per-lane read groups to properly
account for lane-specific batch effects — a detail that's easy to skip
with single-lane demo data but matters on anything with real sequencing
runs behind it.

## Results (demo dataset)

A small multi-lane human dataset (GRCh37 chr1/2/3/8/11/X subset, 5
sequencing lanes, ~4,600 read pairs total), from
[`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/sarek)
(the same small reference + GATK known-sites bundle — dbSNP, Mills
indels — used in nf-core/sarek's own CI tests).

| Metric | Value |
|---|---|
| Total read pairs (5 lanes merged) | 4,602 |
| Mapping rate | 99.98% |
| SNPs called (PASS) | 31 |
| SNPs filtered out (QD/FS/MQ) | 6 |
| Indels called | 0 |

(No indels is expected: the demo region is tiny, and indel calling
needs read support that a ~200kb slice with ~4,600 read pairs mostly
doesn't have — this is a dataset-size artefact, not a pipeline bug. The
hard-filter thresholds and both a SNP and an indel branch are in place
and exercised by the code regardless.)

Every called variant is a real heterozygous SNP with real GATK
annotations — allele depth, genotype quality, strand-odds ratio, the
lot:

```
1  133129  .  G  A  795.64  PASS  AC=1;AF=0.500;...;QD=13.72;...  GT:AD:DP:GQ:PL  0/1:27,31:58:99:803,0,678
```

Full output: [`results/variants/NA_tiny.filtered.vcf.gz`](results/variants/),
[`results/summary/NA_tiny.variant_summary.txt`](results/summary/) (bcftools
stats), and the MultiQC report at `results/qc/multiqc_report.html`
(open locally — GitHub doesn't render HTML inline).

## The web app

`app/` is a small Flask front end over the same variant-calling core:
upload an aligned, coordinate-sorted BAM, it runs HaplotypeCaller and
hard-filtering against the bundled reference, and returns a summary
table plus a downloadable filtered VCF.

```bash
cd app
pip install -r requirements.txt
export FLASK_APP=app.app
flask run --debug
```

It's deliberately a thin layer, not a rebuild of the pipeline logic —
the point is one tested variant-calling core (aligner choices,
filter thresholds) with two front ends for two different jobs: a batch
pipeline for many samples, and a quick single-BAM tool for one. It skips
BQSR (a BAM handed to a one-off web tool has usually already been
through whatever alignment/recalibration its owner wanted), which is
the main functional difference from the full pipeline — noted in the
app's own docstring, not hidden.

Tested end-to-end against the pipeline's own `NA_tiny.recal.bam` output
before being included here — see the "Running it" section for how to
reproduce that.

## Repository structure

```
├── config/
│   ├── config.yaml           # reference paths, filter thresholds
│   └── samples.tsv           # sample/lane sheet
├── data/
│   ├── raw/                  # input FASTQs (not committed — see below)
│   └── reference/            # small GRCh37 reference + GATK known sites
├── workflow/
│   ├── Snakefile
│   ├── envs/environment.yml
│   └── scripts/summarize_variants.py
├── app/                      # Flask BAM→VCF web app
│   ├── app.py
│   └── templates/
├── results/                  # pipeline output (small files committed, BAMs excluded)
└── logs/
```

## Running it

```bash
git clone https://github.com/amaaniu/gatk-variant-calling-pipeline.git
cd gatk-variant-calling-pipeline

conda env create -f workflow/envs/environment.yml
conda activate variant-calling-pipeline

# GATK4 isn't on conda's default channel resolution path for every
# platform — if the conda install doesn't pull it, grab the release
# directly: https://github.com/broadinstitute/gatk/releases

# fetch demo data + reference (see below), then:
snakemake -s workflow/Snakefile -n -p        # dry run
snakemake -s workflow/Snakefile --cores 2 -p # full run, ~1 minute
```

**Memory note:** `fastp` allocates large buffers regardless of input
size. On a memory-constrained machine, cap concurrency
(`--cores 2` rather than a higher number) — this is exactly the kind of
thing that only shows up when you actually run a pipeline rather than
just writing it, and is worth knowing about before it causes a
confusing OOM kill on a shared cluster node.

### Getting the data

```bash
git clone --depth 1 --branch sarek https://github.com/nf-core/test-datasets.git tmp-testdata

mkdir -p data/raw data/reference
cp tmp-testdata/testdata/tiny/normal/*.fastq.gz data/raw/
cp tmp-testdata/reference/human_g1k_v37_decoy.small.fasta* data/reference/
cp tmp-testdata/reference/dbsnp_138.b37.small.vcf.gz* data/reference/
cp tmp-testdata/reference/Mills_and_1000G_gold_standard.indels.b37.small.vcf.gz* data/reference/
rm -rf tmp-testdata

# rename to match config/samples.tsv (NA_tiny_L001_R1.fastq.gz, etc.)
# and reference.fasta / reference.dict / known_sites_*.vcf.gz — see
# config/config.yaml for the exact names expected.

# regenerate the sequence dictionary (the bundled one is occasionally
# a stripped placeholder rather than a full dict):
gatk CreateSequenceDictionary -R data/reference/reference.fasta \
    -O data/reference/reference.dict
```

## Swapping in real data

Nothing in the Snakefile is specific to this demo region. To run this
on a real sample:

1. Update `config/samples.tsv` with real sample/lane names.
2. Point `config/config.yaml` at a full reference (e.g. GRCh38) and
   full-size known-sites files (GATK's [resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)).
3. Run with more cores and more memory — a real WGS/WES sample needs
   considerably more RAM for BWA indexing/alignment and GATK than this
   200kb demo region does.
4. For a real cohort (not a single sample), the natural next step is
   joint genotyping with `GenomicsDBImport` + `GenotypeGVCFs` across
   samples, and VQSR instead of hard filters once there are enough
   variants to train the model on.

## Known limitations

- **Demo scale.** ~4,600 read pairs across a ~200kb region is enough to
  exercise every step of the pipeline correctly, but far below real
  sequencing depth — variant counts here are illustrative.
- **Hard filters, not VQSR.** GATK's recommended approach for a real
  cohort is Variant Quality Score Recalibration, a model trained on
  known-good/known-bad variant sets — it needs many more variants than
  this demo has to train on, so hard filtering (also a standard,
  GATK-documented approach for smaller callsets) is used instead.
- **Single sample, no joint genotyping.** The pipeline calls one sample
  at a time; a real cohort study would use GenomicsDBImport +
  GenotypeGVCFs to jointly call across samples.
- **Flask app has no auth or job persistence** — it's a demo of the
  variant-calling core, not a deployable service.

## Tools

Snakemake · BWA-MEM · fastp · GATK4 · samtools · bcftools · MultiQC ·
Flask · pandas

## License

MIT — see [LICENSE](LICENSE).
