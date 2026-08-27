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
| `--no_cleanup` | Keep intermediate files (both.fa / left.fa / right.fa, normalized norm.fq, etc.; located in the scratch directory) |
| `--scratch <dir>` | Write intermediates to `<scratch>/<output dir name>/` (no space check when given explicitly) |

Resume: each stage writes an `.ok` checkpoint file; rerunning with the same `--output` skips completed stages automatically.

**Default location of intermediates (this port's default, differs from upstream)**: without `--scratch`, all intermediates (normalization dir / both.fa / mer_counts / inchworm / chrysalis / butterfly) go to a local temp directory `${TMPDIR:-/tmp}/trinity-rust-<hash of the output path>/`, and the `--output` directory only holds the final products (`<output>.Trinity.fasta` + gene_trans_map are written next to it as usual). The path is derived deterministically from the output directory, so rerunning with the same `--output` hits the same scratch directory and resume semantics are unchanged. At startup a statvfs space check runs (rough estimate: 6× the total compressed input size, floor 4 GB); if TMPDIR is unwritable or too small, a warning is printed and the run falls back to the old behavior (intermediates inside `--output`). Scratch is treated as volatile: deleting it means recomputing everything; it is not auto-cleaned after a successful run and can be removed manually.

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

**Intermediate layout (performance, accepted difference)**: `chrysalis/Component_bins/Cbin*/` no longer holds per-component `cN.graph.out` / `cN.graph.reads` / `cN.graph.allProbPaths.fasta` files (28k components ≈ 85k small files — prohibitive metadata round-trips on NFS). Instead each Cbin has two packed stores: `graphs.pack`/`graphs.idx` (QuantifyGraph output) and `butterfly.pack`/`butterfly.idx` (Butterfly output); the `.idx` is a fixed-size 32-byte record stream (component id / kind / offset / len, little-endian). The normalization `left.norm.*` convenience names are hard links to the long-named files (upstream uses symlinks) instead of physical duplicate writes.

**Normalization output format (fq input)**: fq input now yields normalized **fa** directly (`left.norm.fa` / `{name}.normalized_K….fa`), skipping the full materialization of `norm.fq` and the downstream fq→fa pass — normalization already converts fq→fa internally, so selected records are extracted from the in-memory fa, byte-identical to "write fq then convert" (locked by tests). Upstream leaves the normalized fq for the user; this port **no longer keeps norm.fq by default** — use `--no_cleanup` if you need it (both fq and fa are then written).

**Known tie-break differences (accepted divergence band)**:

- The inchworm parallel seed tie-breaking order is non-deterministic, exactly like upstream `--PARALLEL_IWORM` — even two runs of the same implementation drift in their transcript multiset (observed coverage band 78–94%); cross-implementation end-to-end bidirectional coverage measured at 57–71% (see the threshold calibration notes in `docs/xcheck-trinity-report.md`). The differences are dominated by end-truncated variants of identical sequences (100%-identity containment), not porting errors.
- Byte-for-byte comparison does not apply to full-pipeline outputs; equivalence is judged by "sequence multiset (with revcomp normalization) + ≥99% clustering" (`cargo xtask eval-trinity`).
- Output traversal order (hash order) is not guaranteed to match upstream; equivalence is always judged by multiset.

## Performance summary

See `docs/benchmarks.md` (per-stage cross-checks and thread scaling) and `docs/opt-log.md` (end-to-end optimization log) for full records.

**End-to-end (real data: Arabidopsis DRR095137 prefix subset, 10M PE reads, ~1.1 GB gz):**

| Version | 8 threads | 64 threads | Notes |
|---|---|---|---|
| v0.0.1 (straight port) | ~54 min | not measured | porting baseline, metrics anchored to upstream |
| v0.3.0 | ~17 min | 12:29 | algorithmic optimizations (see below) |
| **v0.4.0** | **11:20 (-79%)** | **7:37** | IO/system-overhead optimizations |

Small-set end-to-end (upstream sample_data, full 30575 PE reads, 8 threads): this port ~9 s vs upstream (incl. Perl/JVM) ~26 s.

**Per-stage vs upstream Trinity v2.15.2 (equal thread counts, best-of-3):**

| Stage | Upstream wall / peak RSS | This port wall / peak RSS | Speedup | Memory |
|---|---|---|---|---|
| k-mer count+dump (100k reads, 4 threads) | jellyfish 0.95 s / 626 MB | 0.30 s / 86 MB | ~3.2× | ~1/7 |
| inchworm (288k kmers, single-threaded) | 0.72 s / 26.1 MB | 0.37 s / 22.5 MB | ~2.0× | -14% |
| chrysalis full chain (55 components, single-threaded) | 11.6 s / 466 MB | 9.9 s / 147 MB | ~1.2× | ~1/3.2 |
| butterfly small component (c0, N=4342) | 1.15 s / 244 MB | 0.03 s / 9.2 MB | ~40× | ~1/27 |
| butterfly mid components (c1/c2, N≈10k) | 1.7 s / ~420 MB | 0.6 s / ~380 MB | ~2.5× | on par |

Note: the chrysalis ReadsToTranscripts stage initially did not scale with threads (slower than upstream at `-t 4`); its query structure was rebuilt in v0.3.0 (16.3×). The dramatic small-component butterfly gap is mostly upstream JVM startup and baseline overhead.

**Thread scaling (v0.4.0, real data, `--max_memory 64G`):**

| Threads | Wall | Peak RSS | Speedup |
|---|---|---|---|
| 8 | 11:20 | 22.7 GB | 1.00× |
| 16 | 9:12 | 28.5 GB | 1.23× |
| 32 | 8:00 | 33.7 GB | 1.42× |
| 64 | 7:37 | 37.2 GB | 1.49× |

Scaling saturates beyond 32 threads (32→64 yields only 1.05×); peak RSS grows with thread count (k-mer counter concurrent-accumulator effect). At this data scale, 16-32 threads is the sweet spot.

**Optimization breakdown:**

- v0.3.0 (algorithmic): k-mer dump elimination (direct binary bridge, saves ~8 min), ReadsToTranscripts packed hash index + sliding-window entropy (16.3×), parallel coverage statistics + rolling revcomp normalization (3.6×), batched component partition writes (5.5× — the original issued ~45M raw per-line syscalls).
- v0.4.0 (IO/system overhead): normalization double-write elimination via hard links (-2.6 GB written), intermediates on a local scratch dir by default (`--scratch` is now the default; output dir stays clean), full component-data packing (~140k per-component files eliminated, butterfly stage 136 s → 31 s), direct fa emission from normalization, mimalloc as the global allocator (sys time -85%, page faults -96%).

All optimizations carry oracle/byte-level equivalence verification; assembly metrics are identical to the porting baseline (see the tie-break divergence band above).

## Validation

- `cargo test --workspace`: 543 unit tests / component cross-checks, all green.
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
