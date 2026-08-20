# p4 align golden (jaligner NW / NW-banded f32)

- `align_golden_input.tsv`: lines of `mode\ts1\ts2[\tbandwidth]` (mode B=banded, N=plain NW),
  scoring MatrixGenerator.generate(4, -5), gap open 10, extend 1.
- `align_golden.tsv`: per line `aligned1\taligned2\tscore\tstart1\tstart2` from the Java jar.

Regenerate:
```
TRINITY_SRC=/storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2
BJAR=$TRINITY_SRC/Butterfly/Butterfly/Butterfly.jar
javac -cp $BJAR xtask/fixtures-src/JAlignGolden.java   # run from a scratch dir
java  -cp $BJAR:. JAlignGolden < fixtures/p4/align/align_golden_input.tsv > fixtures/p4/align/align_golden.tsv
```
(javac/java = /public/home/senior007/miniconda3/envs/trinity/bin/)
