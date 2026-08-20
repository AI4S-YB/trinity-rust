# 移植映射表（Rust ↔ 原版 Trinity v2.15.2）

持续维护：每个新移植模块追加一行。行号引用以仓库内代码注释实际标注为准
（原版源码树：`/storage/home/senior007/test/trinity_rust/trinityrnaseq-v2.15.2`）。

## trinity-common

| Rust | 原版 | 说明 |
|---|---|---|
| kmer.rs::kmer_to_intval / base_to_int | Inchworm/src/sequenceUtil.cpp:258（表 :10-24） | 2-bit 编码 G=0,A=1,T=2,C=3；小写同表；非 gatc 报错 |
| kmer.rs::decode_kmer_from_intval | Inchworm/src/sequenceUtil.cpp:298 | 低位端逐 2-bit 解出，逆序写入 |
| kmer.rs::revcomp_val | Inchworm/src/sequenceUtil.cpp:181 | `~kmer` 互补 + 循环移位反转 2-bit 组 |
| kmer.rs::get_ds_kmer_val | Inchworm/src/sequenceUtil.cpp:376 | canonical = max(kmer, revcomp(kmer)) |
| kmer.rs::compute_entropy | Inchworm/src/sequenceUtil.cpp:316 | log2 香农熵，f32 逐项累加（浮点序不可变） |
| kmer.rs::KmerId / MAX_KMER_LENGTH | Inchworm/src/sequenceUtil.hpp:20 | `kmer_int_type_t = unsigned long long`，上限 32 |
| fasta.rs::FastaReader/FastaRecord | Inchworm/src/Fasta_reader.cpp:121-122 + Fasta_entry.cpp:22-30 | header/acc 语义 + 序列大写化去空白；4 处已证差异见文件头注释（有意为之） |
| fasta.rs::add_fasta_seq_line_breaks | Inchworm/src/IRKE.cpp:386 | 每 interval 字符插换行，末组不换 |
| fastq.rs::FastqReader | （原版由 Perl 脚本逐行转换，无独立 C++ 读者） | 严格 4 行记录；校验比原版严格 |
| io_util.rs::open_maybe_gz | （无原版对应，纯 Rust 基础设施） | gzip 魔数 1f 8b 嗅探，flate2 流式解压 |
| seq_hash.rs::generate_hash | Inchworm/src/sequenceUtil.cpp:422-445 | sdbm 型 `hash = 65599*hash + c`（u32 回绕）+ base_val 求和（G=1,A=2,T=3,C=4,非 gatc=0）占高 32 位；fold `hash ^= hash>>16`；inchworm 去重 key 接收方再 `as u32` 截断复现 unsigned int 隐式收窄；黄金 fixtures/p2/hash_golden.tsv（C harness 直链 sequenceUtil.cpp） |
| error.rs | sequenceUtil.cpp:260-261 / 276-280 | "error, ..." 消息格式镜像 |

## trinity-kmer

| Rust | 原版 | 说明 |
|---|---|---|
| counter.rs::for_each_kmer / KmerCountTable | jellyfish count（2.3.x）语义 | 多重集等价（xcheck-kmer [1/3]）；小写同大写、非 gatc 断窗续滑、序列段内滤 `\n`/`\r` 跨行拼接、DS 键 = max(kmer, revcomp)；rayon 分片 → 线程局部表 → reduce |
| dump.rs::write_dump | jellyfish dump -L | ">count\nKMER" FASTA；DS 代表取词典序较小串（jellyfish 编码语义）以逐字节对拍；遍历序哈希序不定，等价性靠多重集 |
| fmt.rs::format_g6 | C++ `ostream <<` 默认浮点（≈ printf %g，6 位有效） | "nan"/"-nan"/"-0" 特例；12k 值差分 0 失配 |
| coverage_stats.rs::populate_kmer_counter | Inchworm/src/fastaToKmerCoverageStats.cpp:181-228 | 解析 jellyfish dump 建计数表 |
| coverage_stats.rs（窗口循环） | 同上 cpp:300-335 | 每个滑窗都产覆盖值（非 gatc 窗按 1），与计数器跳过语义不同，勿合并 |
| coverage_stats.rs::median_cov / mean_f32 / stdev_f32 | 同上 cpp:337-347 / 371-387 / 389-402 | 无符号截断中位、long 求和转 f32、样本标准差（n=0 → -0，n=1 → "-nan"） |
| coverage_stats.rs::CoverageStatsRow | 同上 cpp:140-148（stats_text） | 数值逐字段相等（xcheck-kmer [2/3]，58 行） |
| drand48.rs | perl ≥5.20 rand() = POSIX drand48（perl 5.38.2 实测） | 48-bit LCG a=0x5DEECE66D c=0xB；srand(12345) 序列前 1000 值位级一致 |
| nbkc.rs::select | util/support_scripts/nbkc_normalize.pl:67-127（过滤 L81-115） | 随机数消耗序是核心不变量（过滤行不消耗 rand）；NaN 比较语义自然镜像；名单逐 acc 相等 |
| nbkc.rs::merge_pairs | util/support_scripts/nbkc_merge_left_right_stats.pl:89-148（--sorted 路径） | 双指针按 core acc 合并；%.1f 合成列 |
| read_names.rs::update_read_name_for_trinity | trinity-plugins/seqtk-trinity/seqtk.c:201-290（子段 :209-219/:222-235/:243-267/:276-284） | 400 例差分 0 分歧；`_forward` 截断语义 |
| read_names.rs::comp_tab / revcomp | 同上 seqtk.c:163-172 / 538-551 | A<->T、C<->G、U->A/a，大小写保留 |
| read_names.rs::fq_records_to_fa | `seqtk-trinity seq -A -R <1|2> [-r]`（insilico_read_normalization.pl L688-695 调用形态） | 黄金对拍 fixtures/p1/seqtk_names.fa 逐字节相等；空文件 exit(5) 语义（seqtk.c:564-567） |
| read_names.rs::core_read_name | PerlLib/Fastq_reader.pm:113-116 + insilico_read_normalization.pl:534-536 | 抽取侧子串**删除**语义（perl 实测与 seqtk 截断分叉） |
| diginorm.rs::run（together 模式编排） | util/insilico_read_normalization.pl：prep L318-368/L660-719、both.fa L350、count+dump L588-654、stats L816-898、merge L950-985、抽取 L505-584 | 端到端逐字节一致 ×8 配置（xcheck-kmer [3/3]） |
| diginorm.rs（SS 链路接线） | fastaToKmerCoverageStats.cpp:75 + insilico_read_normalization.pl:838 | 计数阶段按 SS 关 --canonical；stats 阶段恒 DS（原版从不传 --SS） |
| diginorm.rs（dump -L 语义） | insilico_read_normalization.pl:46（固定 -L 2） | 计数 < 2 的 kmer 不进查表，查表侧缺失按 1 |
| diginorm.rs（stats 装载） | Inchworm/src/KmerCounter.cpp:476-487（add_kmer） | DS 装载键 canonical 化 + 求和合并 |
| diginorm.rs::extract_fq / extract_fa | PerlLib/Fastq_reader.pm + Fasta_reader.pm（$/="\n>"） | 字节级 4 行扫描 / 块内控制字符删除，命中原样写出 |

## trinity-inchworm

原版行号缩写：KC = Inchworm/src/KmerCounter.cpp。键型 u64（kmer_int_type_t）× u32 计数、
**删除=置 0 不 erase**（size 含 0 值键）、DS 模式 add/find/clear 全先 canonical = max(k, rc(k))。

| Rust | 原版 | 说明 |
|---|---|---|
| kmer_counter.rs::KmerCounter::new | KC:13-17 | kmer_length > 32 断言（消息原样） |
| kmer_counter.rs::canonical | KC:367-375（find_kmer 前置折叠） | DS 折叠 max(kmer, revcomp)，add/get/clear 共用 |
| kmer_counter.rs::add_kmer | KC:476-488 | canonical 后 `counter[key] += count`（unsigned 溢出回绕 → wrapping_add） |
| kmer_counter.rs::get_kmer_count | KC:448-457 | canonical 后查，缺省 0 |
| kmer_counter.rs::clear_kmer | KC:420-436 | 惰性删除：键存在才置 0（:428 find==end 什么都不做，不插入）、不 erase、size 不变 |
| kmer_counter.rs::size | KC:51-53 | 键总数（含 count=0 键） |
| kmer_counter.rs::get_kmer_string | KC:460-465 | 按存储值解码（不做 canonical） |
| kmer_counter.rs::get_kmer_intval | KC:467-473 | 仅编码不折叠（canonical 在 add_kmer :479）；非 gatc → Err 镜像 throw |
| kmer_counter.rs::KmerCatalog::get_forward_kmer_candidates | KC:507-533 | 纯位运算 `(kmer << (33-K)*2) >> (32-K)*2 \| i`（i=0..3 = G,A,T,C 序），count==0 过滤（:528）；K=32/K=1 移位边界按 x86 回绕（wrapping_shl） |
| kmer_counter.rs::KmerCatalog::get_reverse_kmer_candidates | KC:549-575 | `(i << (2K-2)) \| (kmer >> 2)`，过滤与排序同上 |
| kmer_counter.rs（候选排序） | KC:823-831（发布版比较器） | 原版仅 count 降序 + std::sort 不稳定（平局 introsort 实现定义序）→ 显式化为稳定排序，平局保持 G,A,T,C 收集序 |
| kmer_counter.rs::iter_nonzero / iter_all | KC:719-728（收集阶段）/ KC:143（prune 遍历域） | count>0 快照迭代 / 全键含 0 值（C++ map 迭代域；序为 FxHashMap 序 ≠ 原版哈希序） |
| kmer_counter.rs::KmerCatalog（trait） | （原版单一 KmerCounter& 代码，单线程/omp 共享同对象） | 只读视图泛型——单线程 KmerCounter 与 PARALLEL SyncKmerCounter 复用 irke.rs 同一份直译贪心核心；候选生成作默认方法两实现共用 |
| visitor.rs::KmerVisitor | KC:756-799（new :760-764 / add :766-772 / exists :774-780 / erase :782-789 / clear :791-794 / size :796-799） | std::set → HashSet（4^32 键空间无行为级差异）；DS 下 add/exists/erase 先 canonical（:766-789）；erase-重入支持递归回溯 |
| glibc_rand.rs::GlibcRand | glibc stdlib/random_r.c __srandom_r + __random_r（TYPE_3 additive feedback，degree 31/separation 3） | 原版从不 srandom → rand() 即 srand(1) 序列逐位复现；Schrage 初始化（i64 直乘等价）、r[31..33]=r[0..2] 别名展开、34 元环递推 r[i]=r[i-31]+r[i-3]（i32 回绕）、丢弃前 310 输出、返回 `(u32)>>1`；黄金 fixtures/p2/glibcrand_{seed1,mod2}.txt（C harness 位级一致） |
| irke.rs::c_atoi | IRKE.cpp:129（`atoi(header.c_str())`） | glibc `(int)strtol` 语义：跳 C 空白、±、首非数字截断、无数字=0、long 饱和 + (int) 回绕按位解释 u32（黄金值 C harness 锁定） |
| irke.rs::populate_from_kmers | IRKE.cpp:81-154 | jellyfish dump 装载：**空序列记录即终止解析**（:107-108，原版 EOF 哨兵同款判定且不计 parsed）；len!=K 严格不等跳过但计 parsed（:110→:124）；atoi 计数（:129）；非 gatc → Err 镜像 throw（:128） |
| irke.rs::populate_from_reads | IRKE.cpp:157-280（非 reassemble 路径） | hasNext 循环**无空记录终止**（:183-184，与 kmers 模式不同）；len<K 严格 <（:207-209）；每窗口 add_kmer(kmer,1)；reassembleIworm（:211-249）与 PRUNE_SINGLETON_READ_INTERVAL（:255-264，默认 0 无 CLI 入口）不在主线未移植 |
| irke.rs::add_sequence | KC:34-44 + KC:493-505 | 逐窗 substr(i,K)；含非 gatc 窗口被 add_kmer(string) 的 contains_non_gatc 前置 return **静默跳过**（跳过点不在 add_sequence），窗口继续前滑——与 P1 counter 同一窗口集合 |
| irke.rs::write_kmer_count_report | IRKE.cpp:144-147 | 单行 size()（读入后剪枝前；原版固定写 CWD，路径参数化） |
| irke.rs::prune_some_kmers | KC:135-280 | 三段：count<min 即时置 0（:156）；熵<min_entropy 即时置 0（:164）；错误 kmer（:172-258）前向+反向对称，dominant=candidates[0]，`c/dom < r && c/self < r` 双严格 < 入 deletion_list**延迟统一置 0**（:267-273，可重复入列幂等）；遍历域全键含 0 值、序为 FxHashMap 序（≠原版哈希序，即时置 0 影响本趟后续候选查询——by-design） |
| irke.rs::sorted_seed_list / collect_nonzero_seeds / sort_seeds_desc | KC:701-752（收集 :719-728、排序 :739-743、__DEVEL_no_kmer_sort :732-737）+ IRKE.cpp:1332-1338 | count>0 收集；单线程排序 count 降序 + 平局 kmer 值降序（KC:823-831 _DEBUG 比较器的显式化——发布版平局为不稳定排序实现定义序，此处选定确定性全序）；PARALLEL 不排序保持容器迭代序 |
| irke.rs::is_good_seed_kmer | IRKE.cpp:719-760 | 四段闸门按原序严格比较：count==0（:725-727）→ 回文 kmer==revcomp 精确等值（:730-739）→ count<MIN_SEED_COVERAGE 严格 <（:742-748）→ 熵<MIN_SEED_ENTROPY 严格 <（:750-758）；未用的 min_connectivity 形参不移植 |
| irke.rs::build_inchworm_contig_from_seed | IRKE.cpp:764-818 | 实时查种子 count（:769）→ 前向 inchworm('F') → visitor 清空后回放前向路径再 erase 种子（:775-795）→ 反向 inchworm('R') → total = fwd+rev+种子（int+int+unsigned 回绕加 :801）→ join（:1150）；种子不在两条 path 中 |
| irke.rs::inchworm | IRKE.cpp:821-918 | 外层循环：每轮 eliminator 清空（PACMAN=false 恒空）、轮数>num_total_kmers 抛错（:851-853）、visitor.erase 当前 kmer（:864）、step(depth=0, MAX_RECURSION)（:870）；`best.count > 0` 严格 > 才延伸（:878）；step 返回深→浅压栈 path **倒序**消费（:884-897）；当前 kmer←entire 末尾（:898）；count 回绕累加（:901） |
| irke.rs::inchworm_step | IRKE.cpp:933-1125 | 递归贪心体：visited/eliminated → 空（:953-968）；visitor.add（:970）；`depth < max_recurse` 取候选（count>0、!visitor、connectivity 恒过）逐个递归（传 depth+1 与**本层 recurse_cap** :1016）后 visitor.erase 回溯（:1022）；paths(size≥1) 排序（:1040，比较器 :921-929 发布版仅 count 降序 → 稳定排序显式化，平局保持候选 G,A,T,C 收集序）；真 tie（同分且**最远端** path[0] 不同 :1049）：recurse_cap≥50 → **rand()%2 二选一**（:1057-1067，唯一 rand 调用点，glibc 复刻）；elif paths[0] 长>best_path_length → cap++ 继续（:1069-1074，唯一保持 tie 分支）；else 取 paths[0]（:1075-1079）；同分同端点任取（:1081-1093）；无 tie 取 paths[0]（:1095-1098）；单路径取之（:1100-1103）；零路径空（:1104-1106）；尾部 depth>0 才 push 自身 + count+=进入时快照（:1112-1117，种子 depth=0 不入不计） |
| irke.rs::join_forward_n_reverse_paths | IRKE.cpp:1128-1150 | reverse 逆序 + 种子 + forward 原序 |
| irke.rs::reconstruct_path_sequence | IRKE.cpp:1159-1176 | 首 kmer 全串 + 后续各末碱基（substr(len-1,1)）；cov_counter 逐 kmer 快照出参 |
| irke.rs::exceeds_min_connectivity | IRKE.cpp:1179-1213 | 首行 `min_connectivity < 1e5 → true` 短路（:1185-1188，默认 0 恒过）；后续死代码逐行镜像保留 |
| irke.rs::MAX_RECURSION / MAX_RECURSION_HARD_STOP | IRKE_run.cpp:79 / IRKE.cpp:29 | 1 / 50（硬停后 rand 打破） |
| irke.rs::IrkeParams | IRKE.hpp:102-105 + IRKE_run.cpp:101-102 | min_connectivity 0.0（恒过）/ MIN_SEED_ENTROPY 1.5 / MIN_SEED_COVERAGE 2 |
| irke.rs::AssemblyParams | IRKE_run.cpp:90-92 | MIN_ASSEMBLY_LENGTH=25 在 -K 解析（:190）**之前**定值——默认恒 25 不随 -K 变（镜像该 quirk）；MIN_ASSEMBLY_COVERAGE=2 |
| irke.rs::compute_sequence_assemblies（单线程） | IRKE.cpp:426-716 | 主循环：哈希只减不增护栏（:441-449）；omp_set_num_threads(1)（:494-501）；种子实时查 count（:546，排序快照已过期）→ 闸门（:548-551）→ build（:553-554，无 TWO_PHASE——仅 PARALLEL 分支 :549）；avg_cov = (f32)total/(len-K+1) + 0.5 截断（:563，f32 除 + double 加 0.5）；记录条件 len≥L && avg≥cov（:569-573）；**无论是否记录** joined_path 全部 clear（:564-574）；tmp.iworm.fa 中转（:505-517）→ 回读去重 key=generateHash 低 32 位（:624-686）；header `>a{i};{avg} total_counts: {tc} Seed: {seed} K: {K} length: {len}`（:662-670）+ 60 列折行（:671）；rand 进程级 srand(1) 全组装共用一实例 |
| irke.rs::compute_sequence_assemblies_parallel | IRKE.cpp:426-686 的 PARALLEL 分支（:460-620） | 种子**不排序**（KC:732-737，容器迭代序）；`#pragma omp parallel for schedule(dynamic,1000)`（:504）→ rayon par_chunks(1000) 动态分发；per-thread tmp 文件（:474-489）→ 内存缓冲按全局种子序收集去重；目录弱一致读写（无锁 hash_map 竞态）→ SyncKmerCounter；TWO_PHASE（:566-577）extract_best_seed 实时查，zapped（:552-557）放弃不清零；omp_set_num_threads（:464-465）→ 自建 rayon 池；rand 原版多线程交错 nondeterministic → **每 chunk 独立 srand(1)**（同输入同 chunk 划分 → tie 序列确定，语义等价的任意平局打破） |
| irke.rs::parallel_assemble_seed | IRKE.cpp:504-620（循环体） | draft build → TWO_PHASE 重取种子重建（Seed 字段仍主循环快照的原种子 :546 捕获点）→ 记录条件（:585-595）→ 路径清零（:597-618，dashmap 原子置 0 弱一致） |
| irke.rs::extract_best_seed | IRKE.cpp:1348-1371 | `count > best && is_good_seed`（&& 短路严格 >）平局取路径更早者；找不到返回 0 |
| counter_sync.rs::SyncKmerCounter | IRKE.cpp:504-620（多线程共享同一 KmerCounter，C++ 侧无同步——竞态 UB 实践可用） | dashmap 分片 RwLock 承载组装期目录：**单键操作原子、跨键无一致性**（另一线程可读到置 0 前旧值——正是原版无锁竞态语义，仅消除 UB）；clear_kmer 为 KC:420-436 的并发版（键存在才原子置 0、不插不删）；装载+剪枝在单线程完成后 from_counter 整表转入 |
| bin/inchworm.rs | Inchworm/src/IRKE_run.cpp（run_IRKE） | CLI 同名同参：必填校验（:172-189）；参数处理序镜像（:191-406）：-K（:89/:193）、--minKmerCount（:90/:199）、-L（:204）、--min_assembly_coverage（:208）、--monitor（:231）、--keep_tmp_files（:237，no-op——不落 tmp）、--DS（:244）、--min_seed_entropy（:248）、--min_seed_coverage（:253）、--min_any_entropy（:258）、--num_threads（:283，非 PARALLEL 下组装恒单线程 IRKE.cpp:467）、prune_error_kmers 默认 true（:354-359）、--min_ratio_non_error（:362）、--PARALLEL_IWORM/--SINGLE_PHASE（:371-380，后者仅前者下生效 :436-438）；IRKE 构造（:448）→ catalog（:454-472）→ 剪枝（:493-506）→ 组装（:520-536）→ TIMING 行（:541-545，整秒）。stderr 文案镜像：`-reading Kmer occurrences...`（IRKE.cpp:96）、`done parsing N Kmers`（IRKE.cpp:140）、TIMING KMER_DB_BUILDING（:142）、`Pruned N kmers`（KC:266）、reads 模式 `done parsing N sequences`（IRKE.cpp:276）。已知偏差：usage 退出码 2 vs 原版 return 1（仓库 CLI 约定）；gzip 嗅探读入是原版超集 |

## trinity-chrysalis

原版行号缩写：CA = Chrysalis/analysis。行号引用以各 Rust 文件头部/条目
doc 注释实际标注为准。核心怪癖（weld 输出序非契约、BubbleUp +2 合并、
EOF 伪边、PrintSeq 空行、ReadsToTranscripts 最长 run、QuantifyGraph 越界
strncpy 熵 / SortPrint 双 tab）均逐位复刻并有测试锁定。

| Rust | 原版 | 说明 |
|---|---|---|
| dna_vector.rs::nuc_index / plain_table | CA/DNAVector.h:67-84 plain_table、:115-118 NucIndex | A=0,C=1,G=2,T=3,N=4（与 Inchworm 2-bit 表不同，勿混用） |
| dna_vector.rs::bases_to_number / translate | CA/aligns/KmerAlignCore.h:46-64 TranslateBasesToNumberExact::BasesToNumber | m_range=12 双表 → 24-mer 字典序键；N/越界打断 |
| dna_vector.rs::read_fasta | CA/DNAVector.cc:856-971 vecDNAVector::Read | shortName/allUpper 语义；序列行只取行内第一个空白 token（DNAVector.cc:938 AsString(0) quirk）；FlatFileParser token = 空白(' ' 与 '\t')分隔非空段（mutil.cc Tokenize） |
| dna_vector.rs::read_fasta_short_names | 同上（QuantifyGraph.cc:346 调用形态） | `Read(file, false, shortName=true, allUpper=true, ...)` |
| dna_vector.rs::DnaStringStream 等价批量读 | CA/DNAVector.cc:1456-1501 DNAStringStreamFast | 流式读的批量等价物（内存整读，行为级一致） |
| dna_vector.rs::revcomp | CA/sequenceUtil.cc:124-177 | A<->T、C<->G，大小写各自互补 |
| dna_vector.rs::compute_entropy（string 版） | CA/sequenceUtil.cc:326-355 / GraphFromFasta.cc:166-195 | Chrysalis 侧字符串熵（与 Inchworm 位级版不同实现） |
| dna_vector.rs::is_simple / SIMPLE_REPEAT 阈值 | CA/GraphFromFasta.cc:197-214 IsSimple、:21 MIN_KMER_ENTROPY=1.3、:24 MAX_RATIO_INTERNALLY_REPETITIVE=0.85 | `compute_entropy(d) < 1.3`（严格 <） |
| dna_vector.rs::is_simple_repeat | CA/GraphFromFasta.cc:70-160 | mid=len/2、i∈[0,mid) 对称判 repeat |
| dna_vector.rs::simple_halves_with | CA/GraphFromFasta.cc:239-276 SimpleHalves | DISABLE_REPEAT_CHECK=false 默认语义 |
| kmer_align.rs::KmerAlignCore | CA/aligns/{KmerAlignCore.h,KmerAlignCore.cc} | 24-mer 谱对齐索引（m_range 12 × 2 表 = GetWordSize 24，KmerAlignCore.h:16/:293；GetBoundValue=4^12 桶，:26-34） |
| kmer_align.rs::add_data | KmerAlignCore.cc:18-24 → 三参重载 | 逐 contig 入索引 |
| kmer_align.rs::KmerAlignCoreRecord | KmerAlignCore.h:156-204 | 比较只看 (contig, pos) |
| kmer_align.rs::merge_sort_filter | KmerAlignCore.cc:287-340 | 两指针线性交，保 one 的顺序 |
| nonred_table.rs::NonRedKmerTable::set_up_templates | CA/NonRedKmerTable.cc:12-96 SetUp(templ, noNs) | 跨 'X'/含非大写 ACGT 窗口不入表 |
| nonred_table.rs::add_counts_from_reads | NonRedKmerTable.cc:161-200（omp 并行版） | 逐 bundle 逐位置 set_count——count 字段存 bundle 行下标，共享 k-mer 后写者赢（-t 1 下 = i 升序） |
| nonred_table.rs::get_count | NonRedKmerTable.h:40-45 GetCount | miss 返回 0 |
| graph_from_fasta.rs（整体） | CA/GraphFromFasta.cc | Phase 1 :1304-1403、weldmer 计数 :1405-1427、Phase 2 :1435-1754、report :672-747 |
| graph_from_fasta.rs::TOO_SIMILAR | GraphFromFasta.cc:29 | 97 |
| graph_from_fasta.rs::coverage / is_good_coverage | GraphFromFasta.cc:396-408 / :410-425 | name `...total_counts:N...` 抽取；min/max > min_iso_ratio（严格 >） |
| graph_from_fasta.rs::encapsulates | GraphFromFasta.cc:328-347 | 包含判定 |
| graph_from_fasta.rs::align_get_per_id | GraphFromFasta.cc:351-386 | 锚定对齐（无 gap 对角线） |
| graph_from_fasta.rs::WeldableKmer / Weldable | GraphFromFasta.cc:481-518 / :521-577 | flank=(kk-k)/2；weldmer 越界 → false |
| graph_from_fasta.rs::sort_weld_graph | shell `sort -k9,9gr`（Trinity 管线） | 无 -s 时 last-resort 整行比较 tie-break（簇序依赖，s_3/s_4 会互换若丢掉） |
| graph_from_fasta.rs（report 输出序） | GraphFromFasta.cc:672-747（非稳定 std::sort + OMP map 插入） | 本版取确定性 (pool_size 升序, pool_id 升序, 成员插入序)；行序非下游契约（BubbleUp 消费前过 sort_weld_graph） |
| graph_from_fasta.rs（Phase 2 单线程决策） | GraphFromFasta.cc:1435-1754（OMP toasted 无锁竞态） | 按 i↑ j↑ FW 后 RC 命中序串行实现 = 原版 -t 1 确定值 |
| graph_from_fasta.rs（add_scaffolds_to_clusters） | GraphFromFasta.cc:903-995 | 不移植（P3 无 PE scaffolding 输入，scaff_pairs 恒 0，占位注释） |
| bubble_up.rs（整体） | CA/BubbleUpClustering.cc | grow_prioritized_clusters :117-214 + 未聚类补齐 + 长度和过滤 + COMPONENT 块输出；参数默认 :26-27（管线实传 -min_contig_length 200 -max_cluster_size 25） |
| bubble_up.rs（+2 合并判据 quirk） | BubbleUpClustering.cc 合并分支 | `sizeA+sizeB+2 <= MAX`（非 `<= MAX`），与单端 `size < MAX` 不一致，按原样保留 |
| bubble_up.rs（EOF 伪边 quirk） | BubbleUpClustering.cc `while(!in.eof()) getline` | 尾换行多跑一轮空行 → `0 -> 0` 伪边 → iworm 0 可能输出两次（Pool::add 不去重）；split('\n') 尾空串天然复刻 |
| bubble_up.rs（PrintSeq 空行 quirk） | BubbleUpClustering.cc PrintSeq | 80 整倍长序列折行后多一个 "\n" 空行 |
| bundle.rs（整体） | CA/CreateIwormFastaBundle.cc | COMPONENT 块 → `>s_<no> <cov>...` + X 连接多序列一行；过滤后组件号留空洞不重编号 |
| bundle.rs::get_iworm_coverage | CreateIwormFastaBundle.cc:112-137 | `[iworm>a...;cov_...]` 抽 cov；start>end → Err 镜像 exit(5) |
| reads_to_transcripts.rs（整体） | CA/ReadsToTranscripts.cc（413 行） | k=25 固定（cc:130）；min_kmer_entropy 默认 1.5（cc:76） |
| reads_to_transcripts.rs（bundle 索引） | cc + NonRedKmerTable | set_all_counts(-1) 后逐 bundle set_count——count 存 bundle 下标、后写者赢 |
| reads_to_transcripts.rs（read 枚举） | cc | toupper 后正向枚举：熵 < min 跳过、miss(-1) 跳过；!strand 时整条 revcomp 再枚举入同一 comp |
| reads_to_transcripts.rs::format_read_name | CA/DNAVector.cc:1504-1514 formatReadNameString | 保留 '>'、内部空格→'_' |
| reads_to_transcripts.rs::best_component | cc:254-268 最长 run 怪癖 | `comp[j]!=comp[j-1] \|\| j+1==len` 才结算——非末组 run=m-1、末组 run=m-2、大小 1 的组永不成 best（勿"修复"） |
| reads_to_transcripts.rs（pct） | cc:271 | `num_kmer_pos = len-25+1`（仅正向）；`pct = (int)(max/pos*100 + 0.5)`（四舍五入） |
| reads_to_transcripts.rs（行序） | cc multimap<int,int>(best, read_idx) | best 升序、同 best 内 read 下标升序 |
| debruijn.rs（trinity-inchworm crate） | Inchworm/src/{FastaToDeBruijn.cpp, DeBruijnGraph.cpp, DeBruijnGraph.hpp} | F2DB 主流程 :152 createGraphPerRecord（逐捆绑记录建图输出）；4-bit 邻接掩码 DeBruijnGraph.cpp:11-14 |
| debruijn.rs::DeBruijnKmer | DeBruijnGraph.hpp | 节点 = kmer 值 + 1-based id + 覆盖计数 |
| debruijn.rs::add_prev/next_kmer | DeBruijnGraph.cpp:176 / :214 | 取 k 的首/末碱基置位 |
| debruijn.rs::get_prev/next_kmers | DeBruijnGraph.cpp:105 / :145 | prev = 首碱基拼 k 后缀；next = k 左移一位碱基 |
| debruijn.rs::add_sequence / get_kmer_node | DeBruijnGraph.cpp:271 / :341 | 逐 kmer 入图；miss 新建 id=++counter；窗口序列转换 sequenceUtil.cpp:387 |
| debruijn.rs::get_root_kmers | DeBruijnGraph.cpp:639 | 起点 = 无 prev 者 |
| debruijn.rs::to_chrysalis_format | DeBruijnGraph.cpp:556 toChrysalisFormat | Chrysalis 图格式输出；header tokenize（string_util.cpp tokenize：跳前导分隔符、连续分隔符不产空 token） |
| debruijn.rs（kmer 计数） | DeBruijnGraph.cpp:545 kmer_count_comparer | max-heap 优先队列元素 |
| debruijn.rs（图容器） | DeBruijnGraph.hpp std::map | → BTreeMap（kmer 升序遍历） |
| debruijn.rs::contains_non_gatc / replace_non_gatc | sequenceUtil.cpp | 同名语义 |
| quantify.rs（整体） | CA/QuantifyGraph.cc | 读索引 :225-235、第一遍 :367-380、第二遍 :386-473 |
| quantify.rs（读索引 add_kmers） | QuantifyGraph.cc:225-235 | 每 read 每 25-mer 位置入 KmerEntry{read,pos}，含 N 照入、无熵过滤；排序键 = KmerEntryCompare（:186-223）对 read 原始字节从 pos 起比较 25 字符 |
| quantify.rs::KmerEntry / IDS | CA/KmerTable.h:31-59 / :63-135 | ori=-1 表示 revcomp 命中（坐标已翻转） |
| quantify.rs（熵怪癖） | QuantifyGraph.cc:28 + strncpy(&d[1], kmer_length) 越界 | 25 长 sub 从下标 1 拷 25 字节——越界一字节，strncpy NUL 补位等价于 `sub[1..25]` 24 字符 + NUL 进 compute_entropy：NUL 不进分子但计分母（24 计数字符 / 分母 25）；阈值 `< 1.0`（本地常量，非 GFF 的 1.3）低复杂度边跳过。原版该处真 UB（读堆尾巴陈旧字符），见刀口边契约 |
| quantify.rs（二分收集） | QuantifyGraph.cc:237-287 BasesToNumberCountPlus | 拼 25-mer = first[prevNode]+24mer 二分收集；非 strand 下 revcomp 再查 + ori=-1 坐标翻转；第 3 列改写 n1+n2 |
| quantify.rs::ReadsExt | QuantifyGraph.cc:163-184 | 只在文件名**最后 6 个字符**内找 '.'（第 7 位及以外不算） |
| quantify.rs::SortPrint | QuantifyGraph.cc:45-146 | ids 按 (ori,id,start) 升序（KmerTable.h:94-109，ori=-1 组在前）；按 (id,ori) 分组，`lastStart > 组首 start`（严格 >，单一 kmer 位置 read 丢弃）；node2 与 seq 间**双 tab**（fprintf("%s\t") 再补一个）——已与原版产物逐字节对拍 |
| partition.rs（整体） | util/support_scripts/partition_chrysalis_graphs_n_reads.pl（249 行） | 图分区按 Component 行切块；`num_kmers+24 < L` 跳过（25-mer 假定）；Cbin<通过数/1000>；块体为空组件终止读取循环（perl `if(@lines) else undef`） |
| partition.rs::sort_reads_to_components | shell `LC_ALL=C sort -k1,1n -k3,3nr -k2,2` | **必须 C locale**（en_US.UTF-8 collation 会重排 key2） |
| partition.rs（reads 分区） | 同 perl | comp 变化开新文件；登记表没有的 comp 静默丢；每条 `>acc pct\nread\n` |
| partition.rs::component_base_listing | 同 perl :124 | comp id 数值升序，graph.tmp 与 reads.tmp 都存在且非空才输出 `id\tbase` |
| bin/trinity-chrysalis.rs（六子命令） | Chrysalis/bin/{GraphFromFasta,BubbleUpClustering,CreateIwormFastaBundle,ReadsToTranscripts,QuantifyGraph} + Inchworm/bin/FastaToDeBruijn | 参数名对齐原版二进制；错误处理沿用本仓库约定（参数错 exit 2 / 运行错 exit 1） |
| bin/trinity-chrysalis.rs::sort_welds / sort_rtc | shell `sort -k9,9gr` / `sort -k1,1n -k3,3nr -k2,2` | 管线便利命令，stdin → stdout |
| bin/trinity-chrysalis.rs::chrysalis_all | Trinity 管线 shell/perl 串联（Chrysalis 段全链） | 输出布局镜像原版 Chrysalis 目录（bundled_iworm_contigs.fasta / Component_bins/Cbin*/c*.{graph,reads}{,.out}）；QuantifyGraph 成功后删输入（QuantifyGraph.cc:489-493 unlink）照搬 |

## trinity-butterfly

原版行号缩写：TA = Butterfly/Butterfly/src/src/TransAssembly_allProbPaths.java（15,920 行）。
行号引用以 crates/trinity-butterfly/src 各文件头注释实际标注为准。
**Oracle 裁定**：发布版 jar `$TRINITY_SRC/Butterfly/Butterfly.jar`（md5 312f…），
源码树与内层 jar 在 combinePaths / DFS_add_path 上有偏差——详见 docs/setup.md
「Butterfly.jar 裁定」节。

| Rust | 原版 | 说明 |
|---|---|---|
| context.rs::BflyContext | TA L49-170 + SeqVertex.java 静态 tracker | 全部 Java 静态字段收编为实例上下文：LAST_ID/LAST_REAL_ID、EDGE_THR=0.02/FLOW_THR=0.02、MAX_MM_ALLOWED（每 read 重写→穿线上下文参数）、MAX_READ_SEQ_DIVERGENCE=0.05 等 |
| graph.rs::DiGraph/SeqVertex/SimpleEdge | JUNG DirectedSparseGraph + SeqVertex.java + SimpleEdge.java | JUNG 语义：禁平行边、O(1) findEdge；圈标记/loop 计数/repeat unroll 权重；节点相等按 name 内容 |
| graph_io.rs | TA preProcessGraphFile L12694 / buildNewGraphUseKmers L12742 / getReadStarts+readAndMapSingleRead 解析 L11770 / writeDotFile L12837 | graph.out 逐行解析（header token）；graph.reads 的 `endInRead = fields[3] + KMER_SIZE` off-by-one FIXME 保留；DOT 检查点写出（java_math_round 复刻 `(int)Math.floor(x+0.5)`） |
| prune.rs（剪枝链 + main 调用序） | TA L780-890 | fixExtremelyHighSingleEdges(L8787, EXTREME_EDGE_FLOW_FACTOR=200) → removeLightEdges 系(L13010-13222, 严格 `<` vs In/Out 版 `<=` 且不重算 total) → compactLinearPaths(L12936, 迭代中改图快照+while) → removeSingleNtBubbles(L8603, 平局保 v2、addToPrevIDs 的 `id >= lastRealID` 疑似 bug 方向保留) → calcSubComponentsStats(L13320, avgCov < 1-0.5 清小组件) |
| dfs.rs::run_dfs2 | My_DFS.java | 双向 visitVertex2 + down-only 两阶段（最终深度来自第二阶段）+ finish time 降序拓扑序；顶点迭代序 JUNG HashSet 任意序 → 确定插入序（depth 对迭代序不敏感，down-only 平局 (depth,插入序) 稳定排序） |
| align.rs | jaligner 内嵌源 NeedlemanWunschGotoh{,Banded}.java / Alignment.java / AlignmentStats.java / ZipperAlignment.java + TA 的 NWalign 封装 | 全 DP f32（Java float，含 `==` 平局）；8 条 QUIRK 逐字保留：行末 vDiagonal=0 重置、首行宽松初始化、报告分 v[n-1] 与 traceback 起点可不同、带外 DIRECTION_INIT 提前截断、Zipper 右锚双 i 下标 bug 等；52 组黄金向量（fixtures/p4/align，Java driver 直调 jar 内 jaligner 类）tests/align_golden.rs 全等 |
| threading.rs | TA L11766-12687 | getReadStarts → readAndMapSingleRead → findPathInGraph → updatePathRecursively（递归 + memo 深拷贝 + 后继平局 `<=` 取迭代序末 + tied 直接改子对象）；三级比对切换（zipper / 短 NW 预检 / 带状 NW）；loc_in_node 精读：段计数被直接当碱基下标用（Java 原文如此）；L4336 环节 tests/threading_c0.rs 全等 |
| pair_paths.rs | PairPath.java + TA getSuffStats_wPairs L9294-9410 | PairPath equals/hashCode = 双路径整体相等（派生）；**combinePaths 合并为 jar 实证行为**——源码树该调用被注释，发布版 jar 字节码实际执行（invokestatic ×2），Rust 按 jar 语义移植（详见 setup.md 裁定）；Dijkstra isAncestral 可达性以闭包注入 |
| pog.rs | TA L1617-3321（create_DAG_from_OverlapLayout L1617 / construct_path_overlap_graph L3321 / break_cycles_in_path_overlap_graph L3147 / convert_path_DAG_to_SeqVertex_DAG L2146 / zipping 调度 L2220-2330 / zip_up/zip_down L2573/2621 / attempt_zip_merge L2671 / remove_containments L3460 / find_dispersed_repeat_nodes L1865 等）+ Path.java + PathWithOrig.java + TopologicalSort.java | **DFS_add_path_to_graph (L2842) 按 jar 新版移植**（源码树为旧版，c0 对拍据此调和）；JUNG/HashMap 迭代序不可复现 → POG 对拍按结构比较（路径内容多重集 + 内容对边集 / orig-id 节点与边多重集）tests/pog_c0.rs |
| paths.rs | TA getAllProbablePaths L9600-10084 + 前后编排（reorganizeReadPairings/triplet/addSandT/remove_identical_subseqs/remove_short_seqs） | 组件内路径搜索主体；DijkstraDistanceWoVer（JUNG Dijkstra 无顶点版本语义，PairPath gap 检查供参）；-L 200 -F 10000 -R 2 默认 → PATH_REINFORCEMENT_DISTANCE=2500；tests/paths_c0.rs |
| postprocess.rs | TA L15392 reduce_cdhit_like / L10608 twoPathsAreTooSimilar / L10748 getPrevCalcNumMismatches / L10683 findLastSharedNode / L9129 group_paths_into_genes / L8954 printFinalPaths / L1438 convert_to_orig_ids / L15663 get_pathName_string / L7000 assignCompatibleReadsToPaths / L11130 removeTheLesserSupportedPath | CD-HIT 式去冗余 + `_g{i}_i{j}` 基因分组 + MISO `path=[id:j-k]` 输出；**EM（期望最大化读支持估计）完整移植**；6 条 Java quirk 保留（单侧空路径不写缓存 / 已删 path_i 继续当过滤证据不 break / ALL_VS_ALL_MAX_DP_LEN=1000 临时上限 / per-id 只算 mismatches / DIFFS_WINDOW 字段存在但未用 / 输出序按 java_hashmap_order 仿真：List 31 进制 hashCode + spread 桶序） |
| bin/butterfly.rs | TA main + Butterfly CLI 参数面 | -N/-L/-F/-R/-C/-V/--NO_EM_REDUCE 等同名同参；组件单线程 + 大栈（默认 512MB 可配，镜像原版进程级并行语义） |

## trinity-cli

| Rust | 原版 | 验证 |
|---|---|---|
| args.rs::parse_args | `Trinity` GetOptions L660-760（主线 17 个同名 flag + 逗号列表 L1127 展开 + 校验 L1088-1130/L969-985/L1156 `--max_memory xG` 镜像） | 未知参数 "do not understand option"；PE/SE × SS 组合校验单测 |
| checkpoint.rs::run_with_checkpoint | Pipeliner.pm（exists → skip；成功 touch）| stderr `checkpoint found, skipping` 与原版同名同位 `.ok`（断点续跑演示见 benchmarks.md P5 节） |
| prep.rs::prep_reads | `Trinity` prep_seqs L2743-2879（fq 经 seqtk `-A -R 1\|2`、`split(//,SS)` RF→L='R'/R='F'、仅字符=='R' 才 revcomp、字节数校验、both.fa 先左后右） | both.fa 互喂抽查与原版记录多重集完全相等（xcheck-trinity [互喂]） |
| orchestrate.rs::run_trinity | `Trinity` 主流程 L1376-1990（归一化 L1441-1453 → prep → jellyfish count/dump → inchworm → iworm 重命名 L1727 → chrysalis（GraphFromFasta→BubbleUp→Bundle→RTT→F2DB→partition→QuantifyGraph）→ butterfly → harvest） | 端到端对拍 xcheck-trinity（fast 50000: 双向覆盖率 57.4%/70.7% ≥ 阈值 50% PASS）；jellyfish histo 及其 `.ok` 未复刻（无下游消费，见 concerns） |
| butterfly_pool.rs::run_butterfly_pool | ParaFly 等价（component_base_listing 读取 + 组件级线程池 + 失败组件收集非零退出） | rayon 池统一（T4）；组件并行 vs 单线程输出一致性见 benchmarks.md P5-T4 |
| harvest.rs::harvest | `print_butterfly_assemblies.pl` + `get_Trinity_gene_to_trans_map.pl`（收集 c*.graph.out.flower → `<out>.Trinity.fasta` + `.gene_trans_map`，outdir 上一级同名拼接） | xcheck-trinity 判定行：双方产出 Trinity.fasta + gene 数对照（58/58） |
| orchestrate.rs::memory_guard_error/warn_memory | 原版无此护栏（Rust 版扩展） | 单测 + benchmarks.md RSS 基准（0.40× max_memory） |

<!-- c2 已知差异归因（P4 终审）：[2957,799] len 2182 短异构体差异根因 = get_all_probable_paths 裸边保底播种的迭代序敏感性（JUNG HashSet 后继序 + 等支持度路径 sort+reverse 平局序，两者叠加使"哪条边被播种"两侧合法不同）；em/noem 两形态独立复现同一差异、黄金序列全覆盖、只多不丢——保守方向，非移植缺口。 -->

## 验证资产（非移植，供追溯）

| 资产 | 对拍对象 |
|---|---|
| xtask xcheck-kmer [1/3] | jellyfish count+dump（DS/SS 双模式） |
| xtask xcheck-kmer [2/3] | 原版 Inchworm/bin/fastaToKmerCoverageStats |
| xcheck-kmer [3/3] | 原版 insilico_read_normalization.pl（DS/SS × max_cov） |
| trinity-inchworm tests/smoke_vs_original.rs（重放） | 原版 `--monitor 2` 抓取的种子序重放 populate→剪枝→主循环：smoke fixture、sample_data 全量 kmers（288,470）、--reads 模式三组 stdout **BYTE-MATCH**——贪心核心/tie 打破/rand/清零/去重/格式化位级复刻，默认种子序的残余分歧全部归因种子平局序 |
| xtask xcheck-inchworm [1/4] | 原版 Inchworm/bin/inchworm 单线程（rc 不变 + header 去 aN 多重集；smoke 全等） |
| xcheck-inchworm [2/4] | 原版 --PARALLEL_IWORM 多线程对拍（双方同为竞态，严格全等重试至多 10 次；大输入降为差异率统计） |
| xcheck-inchworm [3/4] | 原版 --reads 模式（默认 prune_error_kmers=true 路径） |
| xcheck-inchworm [4/4] | 原版 Chrysalis/bin/GraphFromFasta 消化我们的输出（-k 24 -kk 48 主脚本默认，exit 0 + weld 图非空——P2 门） |
| fixtures/p1/* | seqtk-trinity、perl 5.38.2（rand 序列）、原版管线 golden |
| fixtures/p2/* | generateHash C harness、glibc srand(1) 100+50 值黄金、原版 inchworm 冒烟产物与种子序 TSV |
| trinity-chrysalis tests/gff_vs_original.rs | 原版 GraphFromFasta（-t 1 消除 OMP toasted 竞态，Trinity:2180 实参形态）：边多重集（(A,B) 有向对 + weldmers/total/min_len 全字段）**完全相等**（24 边）；行序差异在预期内（非下游契约） |
| trinity-chrysalis tests/bubble_vs_original.rs | 原版 BubbleUpClustering + CreateIwormFastaBundle：COMPONENT 块多重集**完全相等**（55 组件）；bundle 输出**逐字节相等**（110 行） |
| trinity-chrysalis tests/quantify_vs_original.rs | 原版 partition perl + QuantifyGraph：分区产物逐字节相等（55 组件）；graph.out **刀口边白名单**内逐行相等（c0 8 行 / c2 4 行 / c1 0 行仅第 3 列计数可差），reads 输出过滤刀口边提及行后逐行相等 |
| xtask xcheck-chrysalis [1-7] | [1] GFF 边多重集、[2] bubble+bundle、[3] rtt+sort 排后逐行、[4] F2DB 块行多重集、[5] partition 逐字节 + listing、[6] quantify 刀口白名单、[7] chrysalis-all 自一致（listing 逐字节 55 组件 + c0 graph.out 自一致）+ Butterfly 冒烟 8 转录本 |
| 刀口边契约 | QuantifyGraph `strncpy(&d[1], kmer_length)` 越界读堆尾巴一字节使熵第 25 字符在 {NUL,A,C,G,T} 间浮动——poly-A/GA 等熵 ≈1.0 的边跳过与否随堆历史翻转（真 UB，不可精确复刻）；本库取确定语义（NUL 补位），非刀口边与原版完全一致 |
| fixtures/p3/* | 原版六程序黄金：gff.welds.orig.txt / bubble.orig.out / bundle.orig.fa / rtt.orig.out / f2db.orig.txt / quantify/c{0,1,2}（含 c*.graph.tmp、reads.tmp、orig.graph.out、orig.reads.out）；iworm.fa 为黄金向量勿再生 |
| trinity-butterfly tests/（8 层对拍，P4 门） | ① checkpoint_a.rs：建图 vs jar DOT 检查点 A；② checkpoints_bchd.rs：剪轻边/压缩+DFS/SNP 坍缩/小组件 vs 检查点 B/C/H/D；③ align_golden.rs：52 组 jaligner 黄金向量全等；④ threading_c0.rs：read 穿线 vs jar L4336 中间态全等；⑤ pair_stats_c0.rs：getSuffStats_wPairs（含 jar 实证 combinePaths）输出一致；⑥ pog_c0.rs：POG/SeqVertex DAG 结构比较一致（jar 版 DFS_add_path）；⑦ paths_c0.rs：路径搜索输出一致；⑧ 端到端见下行 |
| xtask xcheck-butterfly | c0/c1/c2 × {em,noem} = 6 检查点，对拍**发布版 jar**（312f…，见 setup.md 裁定）：3 PASS + 3 PASS-WARN（c0/em 序列多重集全等仅 header/顺序差 = Java HashMap 迭代序；c2 两形态各多 1 条 [2182] 短异构体——jar 的路径搜索/过滤丢弃之，黄金 4/5 条序列全覆盖） |
| fixtures/p4/* | jar 产黄金：align/（jaligner 52 向量）+ c{0,1,2}/allprobPaths.{em,noem}.fasta（发布版 jar 固化，勿再生） |
| xtask xcheck-trinity | **第 3 层端到端**：fast(50000)/full 双侧全管线对拍 + eval（精确 rc 多重集交集 + ≥99% nw_gotoh 聚类 → 双向覆盖率）+ both.fa 互喂 + SS(RF) 合成小集；报告 docs/xcheck-trinity-report.md |
| cargo xtask eval-trinity | 双向覆盖率统计器（阈值 50%，校准依据：同实现自对拍带 78~94%、跨实现 57~71%，差异主体为两侧并行种子平局序） |
