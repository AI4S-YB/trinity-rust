# trinity-rust

English | **[中文](README.md)**

A Rust port of the main assembly pipeline of the Trinity RNA-Seq de novo assembler (upstream [Trinity v2.15.2](https://github.com/trinityrnaseq/trinityrnaseq)) — from in silico read normalization to Butterfly transcript output, end-to-end in a single binary with no external Perl/Java/C++ dependencies (jellyfish / seqtk / ParaFly / Butterfly.jar are no longer needed).

## Usage

```bash
cargo build --release
target/release/trinity-cli \
  --seqType fq --left reads.left.fq.gz --right reads.right.fq.gz \
  --CPU 8 --max_memory 2G --output out
# outputs: out.Trinity.fasta + out.Trinity.fasta.gene_trans_map
```

Common options (same names as the upstream `Trinity` script):

| Option | Description |
|---|---|
| `--seqType fq\|fa` | Input type (`.gz` decompressed automatically) |
| `--left/--right` | PE reads (comma-separated for multiple files); `--single` for SE |
| `--SS_lib_type F\|R\|RF\|FR` | Strand-specific library (RF/FR = PE) |
| `--CPU/--max_memory` | Thread count / memory guardrail |
| `--KMER_SIZE` (alias `__KMER_SIZE`, default 25) | inchworm k-mer size |
| `--min_kmer_cov` (alias `--min_kmer_count`) | k-mer count cutoff |
| `--normalize_max_read_cov` / `--no_normalize_reads` | Normalization cap / skip normalization |
| `--bfly_stack_mb` | Butterfly thread stack (extension of this port; defaults to match upstream JVM stack behavior) |
| `--no_cleanup` | Keep intermediate files (both.fa / left.fa / right.fa etc.) |

Resume: each stage writes an `.ok` checkpoint file; rerunning with the same `--output` skips completed stages automatically.

A `--recursive` mode (partitioned re-assembly) is also available.

## Architecture

Five crates plus an orchestrating CLI; the data flow mirrors the upstream pipeline top to bottom:

```
                 trinity-cli (orchestration: CLI surface / checkpoints /
                              per-stage invocation / harvesting)
                                    |
   normalization (diginorm, K25 maxC200) → both.fa
                                    |
   trinity-kmer        jellyfish count/dump equivalent (counting + -L filter
                       + multi-file merge)
                                    |
   trinity-inchworm    linear contig assembly (k-mer seed extension,
                       PARALLEL mode)
                                    |
   trinity-chrysalis   GraphFromFasta → BubbleUpClustering → ReadsToTranscripts
                       → FastaToDeBruijn → QuantifyGraph → component partition
                                    |
   trinity-butterfly   per-component graph search + EM + path output
                       (allProbPaths)
                                    |
   harvest             <out>.Trinity.fasta + gene_trans_map
                                    |
   trinity-common      low-level primitives: 2-bit k-mer encoding /
                       FASTA·FASTQ reading / sdbm seq_hash / drand48 /
                       seqtk read renaming, etc.
```

## Differences from upstream (established list)

**Not ported (explicitly unsupported)**: salmon abundance estimation, bowtie (`--no_bowtie` semantics always hold), jaccard-clip read clustering, long reads (`--long_read`), DNA mode (DNA assembly outside `--genome_guided`). These entry points either error out directly or are ignored with a warning.

**Option aliases**: `__KMER_SIZE` → `--KMER_SIZE`, `--min_kmer_cov` → `--min_kmer_count` (both spellings accepted, compatible with the upstream main program).

**Extension**: `--bfly_stack_mb` explicitly controls the Butterfly worker thread stack limit (upstream derives this implicitly from the JVM `-Xss` setting).

**Known tie-break differences (accepted divergence band)**:

- The inchworm parallel seed tie-breaking order is non-deterministic, exactly like upstream `--PARALLEL_IWORM` — even two runs of the same implementation drift in their transcript multiset (observed coverage band 78–94%); cross-implementation end-to-end bidirectional coverage measured at 57–71% (see the threshold calibration notes in `docs/xcheck-trinity-report.md`). The differences are dominated by end-truncated variants of identical sequences (100%-identity containment), not porting errors.
- Byte-for-byte comparison does not apply to full-pipeline outputs; equivalence is judged by "sequence multiset (with revcomp normalization) + ≥99% clustering" (`cargo xtask eval-trinity`).
- Output traversal order (hash order) is not guaranteed to match upstream; equivalence is always judged by multiset.

## Performance summary

See `docs/benchmarks.md` (per-stage cross-checks) and `docs/opt-log.md` (end-to-end optimization log) for details.

**Compared to upstream Trinity v2.15.2 (C++/Perl/Java toolchain):**

- k-mer count+dump: ~3.2× faster than jellyfish count+dump at ~1/7 peak memory, with byte-equivalent output multisets (sample_data PE 100k reads, 4 threads).
- inchworm: ~2× faster, ~14% lower peak RSS (single-threaded).
- chrysalis full chain: ~1/3.2 peak RSS (147 MB vs 466 MB).
- butterfly: nearly an order of magnitude faster on small components with ~27× lower peak RSS (no JVM startup/baseline overhead).
- End-to-end full sample_data run (30575 PE reads, 8 threads): this port ~9 s vs upstream (incl. Perl/JVM) ~26 s.

**Compared to this project's v0.0.1 straight port (end-to-end speed work in the current version):**

- Arabidopsis real data (PE, `--CPU 8`) end-to-end wall time **54 min → 17 min (-69%)**, with assembly metrics identical to baseline.
- Four optimizations: k-mer dump elimination (direct binary bridge, saves ~8 min), ReadsToTranscripts query structure 16.3×, normalization 3.6×, component partition writing 5.5× — all verified oracle/byte-level equivalent.

## Validation

- `cargo test --workspace`: 534 unit tests / component cross-checks, all green.
- Three-layer cross-validation (`cargo xtask`):
  - `xcheck-kmer / xcheck-inchworm / xcheck-chrysalis / xcheck-butterfly` — layers 1 and 2 (unit/pipeline level, checked against upstream binaries or goldens);
  - `xcheck-trinity [--full]` — layer 3 end-to-end: truncated (default 50000 PE reads) / full-scale full-pipeline cross-checks in both directions + eval statistics + both.fa cross-feed spot checks + SS(RF) synthetic small set; decision thresholds and calibration rationale in `docs/xcheck-trinity-report.md`;
  - `eval-trinity <ours> <orig>` — standalone eval report (exact match / 99% clustering / bidirectional coverage / gene counts).

Upstream toolchain environment (for cross-checking): `TRINITY_SRC` (upstream source tree), jellyfish/java in conda env `trinity`; see `docs/setup.md`.

## Development docs

- `docs/porting-map.md` — function-by-function Rust ↔ upstream source mapping table (with line numbers and proven differences)
- `docs/benchmarks.md` — benchmark records
- `docs/backlog.md` — backlog and decision log
- `docs/setup.md` — environment and upstream toolchain preparation
- `docs/plans/` / `docs/specs/` — stage implementation plans and design documents (P0–P5 porting, M0/M1 benchmarking & recursive partitioning, M7 speed optimization)
