// 黄金向量生成器: 直接链接原版 sequenceUtil.cpp，从 stdin 逐行读序列，
// 输出 TSV: seq \t generateHash 的 u64 值。
// 空行是合法输入（generateHash("") = 0），不跳过——输入文件中间的空行即空串样本。
// 编译: g++ -O2 -I$TRINITY_SRC/Inchworm/src 本文件 $TRINITY_SRC/Inchworm/src/sequenceUtil.cpp \
//       $TRINITY_SRC/Inchworm/src/stacktrace.cpp -o dump_hash_golden
#include <cstdio>
#include <iostream>
#include <string>
#include "sequenceUtil.hpp"

int main() {
    std::string line;
    while (std::getline(std::cin, line)) {
        unsigned long long h = generateHash(line);
        printf("%s\t%llu\n", line.c_str(), h);
    }
    return 0;
}
