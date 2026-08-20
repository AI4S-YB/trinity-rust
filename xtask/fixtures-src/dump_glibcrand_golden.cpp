// 黄金向量生成器: glibc random()（TYPE_3，degree 31 / separation 3）。
// 独立小程序，不链原版 Trinity 源——oracle 就是本机 glibc 的 srand/rand 本身。
// 原版 inchworm 从不调 srandom → rand() 即 srand(1) 序列，故固定 srand(1)。
// 用法: dump_glibcrand_golden raw   → srand(1) 后 rand() 前 100 值
//       dump_glibcrand_golden mod2  → srand(1) 后 rand()%2 前 50 值
//   （两模式各自独立 srand(1)，序列从同一点起）
// 编译: g++ -O2 本文件 -o dump_glibcrand_golden
#include <cstdio>
#include <cstdlib>
#include <cstring>

int main(int argc, char **argv) {
    bool mod2 = argc > 1 && strcmp(argv[1], "mod2") == 0;
    srand(1);
    int n = mod2 ? 50 : 100;
    for (int i = 0; i < n; i++) {
        printf("%u\n", mod2 ? static_cast<unsigned>(rand() % 2)
                            : static_cast<unsigned>(rand()));
    }
    return 0;
}
