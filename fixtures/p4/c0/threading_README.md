# p4 c0 read 穿线黄金（getReadStarts）

`threading_golden.tsv`：每行 `read名\t路径(逗号分隔节点id)`，顺序 = reads 文件中
"Threaded Read as" 出现序（-V 17 --stderr 的 err 输出逐行解析）。
`threading_stats.txt`：jar stderr 的统计行（"number of reads threaded = ..."）。

再生成（scratch 目录内）：
```
TRINITY_SRC=/storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2
FIX=fixtures/p3/quantify/c0
ln -s $FIX/orig.graph.out c0.graph.out
ln -s $FIX/orig.reads.out c0.graph.reads
java -Xmx4G -jar $TRINITY_SRC/Butterfly/Butterfly/Butterfly.jar \
  -N 4342 -L 200 -F 10000 -R 2 -C c0.graph -V 17 --stderr > out.txt 2> err.txt
tr '\r' '\n' < err.txt | grep "Threaded Read as" \
  | sed -E 's/^Threaded Read as: (\S+) : \[(.*)\]$/\1\t\2/' | sed 's/, /,/g' > threading_golden.tsv
tr '\r' '\n' < err.txt | grep "number of reads threaded" > threading_stats.txt
```
（javac/java = /public/home/senior007/miniconda3/envs/trinity/bin/）
