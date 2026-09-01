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
| `--KMER_SIZE` (alias `__KMER_SIZE`, default 25) | inchworm k-mer size (**ignored when multi-k is on**) |
| `--multi-k auto\|k1,k2,...` (default auto) | **Multi-k iterative assembly**: trusted contigs from lower-k rounds are injected into higher-k rounds (k chosen by the rnaSPAdes formula, supported up to k≤63); `--no-multi-k` falls back to single-k |
| `--inchworm_cpu` (default min(6,CPU)) | inchworm parallelism (upstream-faithful default; raise it on many-core machines) |
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

## Multi-k iterative assembly (default since v0.5.0)

Upstream Trinity builds the graph with a single k-mer (k=25): k-mers of low-coverage transcripts fall below pruning thresholds and are lost, while high-coverage regions suffer from excessive branching. This port adopts the rnaSPAdes multi-k strategy: `--multi-k auto` picks the k sequence from the read length (`lower=int(RL/3)`, `upper=int(RL/2)-1`, rounded down to odd, clamped to the supported domain k≤63), runs the lower-k rounds first (graph building only, no transcript output), and injects all k-mers of their inchworm contigs into the next round's count table in presence-only form (count=1, leaving depth-pruning semantics untouched). `--no-multi-k` restores the upstream single-k behavior; an explicit `--multi-k 29,49,63` overrides the sequence (2–3 values recommended).

Companion change: k-mer keys widened from u64 to u128 (`KmerId`), supporting k≤63 (was 32); the counts_bin cache format is versioned as TRM2 (old TRMC files are rejected and recomputed — do not reuse workdirs across versions).

Measured recall gains (R_full = fraction of reference transcripts assembled at ≥95% coverage and ≥95% identity; see `docs/opt-log.md`):

| Data | single k=25 | two-k | three-k |
|---|---|---|---|
| Arabidopsis simulated (76bp, cross-validated on multiple datasets) | baseline | +1.3~2.0pp, P −0.1~−2.2pp | — |
| Real RNA-seq (150bp, 21.35M PE pairs, ~436Mb-genome species) | 31.5% | **43.9%** ([49,63], +12.4pp) | **45.3%** ([29,49,63], another +1.4pp) |

Longer reads and more complex transcriptomes benefit more; on 150bp real data the recall gain is +44% relative.

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
| v0.4.0 | 11:20 (-79%) | 7:37 | IO/system-overhead optimizations |
| **v0.5.0** | see note below | — | multi-k by default + hotspot optimizations |

**v0.5.0 end-to-end (larger real data: 21.35M pairs of 150bp PE reads, ~436Mb-genome species, 32 threads):**

| Version | wall | key changes |
|---|---|---|
| v0.4.0 basis | ~137 min | — |
| **v0.5.0** | **~31 min (4.2×)** | multi-k [49,63] iteration + hotspot optimizations below |

v0.5.0 hotspot optimizations (all verified byte-identical or golden-file equivalent; assembly metrics unchanged):

- **kmer_load elimination**: the inchworm count table is passed in-memory and move-constructed (skipping the ~4 GB counts_bin read-back and a second hash-table build; ~10 min → ~0 per round).
- **GraphFromFasta**: zero-clone two-pointer candidate intersection; the deliberately single-threaded Phase 2 split into "parallel collection of pure-function events + serial in-order application" (order-sensitive toasted/pools semantics preserved); weldmer counting via a 2-bit hash index (2.7 G binary searches → hash lookups, 314 s → 74 s); flattened template collection. Total 20 min → 2.9 min.
- **FastaToDeBruijn**: per-bundle independent graph construction parallelized (95 s → 5 s).
- **Normalization**: two-sided parallel decompress/prep/emit + copy-free counting (both.fa copy eliminated, smaller-into-larger merges), 52 min → 12.5 min.
- **inchworm threads**: `--inchworm_cpu` (default min(6,CPU), upstream-faithful) should be raised on many-core machines — two extension rounds 26.5 → 6 min.

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

**IO overhead work (v0.4.0, before/after on real data at 64 threads):**

| Metric | v0.3.0 | v0.4.0 | Change |
|---|---|---|---|
| Intermediate write volume | ~13.4 GB | ~9.6 GB | -29% |
| Per-component small files | ~140k files | 4-6 pack files per Cbin | ≈ eliminated |
| sys time share | 2056 s (31%) | 561 s (13%) | -73% |
| Minor page faults | 197M | 35.4M | -82% |
| butterfly stage wall | 136 s | 31 s | -77% |
| wall (scratch is default) | 12:29 | 7:23 | -41% |

Background: assembly output directories often live on network storage (NFS), where the upstream-style "one set of small files per component + full-dataset file handoff between stages" turns every create/write/close into a network metadata round-trip. v0.4.0 attacks this three ways: moving all intermediates off the network disk (`--scratch` by default, see "Default location of intermediates" above), packing small files into sequentially-written pack files with fixed-size offset indexes, and hard-linking duplicate content. Peak RSS is unaffected (flat before/after mimalloc).

All optimizations carry oracle/byte-level equivalence verification; assembly metrics are identical to the porting baseline (see the tie-break divergence band above).

## Validation

- `cargo test --workspace`: 566 unit tests / component cross-checks, all green.
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
