// 黄金向量生成器: 直接链接原版 sequenceUtil.cpp，从 stdin 逐行读 kmer，
// 输出 TSV: kmer \t intval \t revcomp \t dsval \t entropy
// 编译: g++ -O2 -I$TRINITY_SRC/Inchworm/src 本文件 $TRINITY_SRC/Inchworm/src/sequenceUtil.cpp \
//       $TRINITY_SRC/Inchworm/src/stacktrace.cpp -o dump_kmer_golden
#include <cstdio>
#include <iostream>
#include <string>
#include "sequenceUtil.hpp"

int main() {
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        kmer_int_type_t v = kmer_to_intval(line);
        unsigned int k = (unsigned int)line.length();
        printf("%s\t%llu\t%llu\t%llu\t%.9g\n",
               line.c_str(),
               (unsigned long long)v,
               (unsigned long long)revcomp_val(v, k),
               (unsigned long long)get_DS_kmer_val(v, k),
               compute_entropy(v, k));
    }
    return 0;
}
