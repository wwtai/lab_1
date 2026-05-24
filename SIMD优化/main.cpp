#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <algorithm>
#include <arm_neon.h>
#include <stdint.h>
#include <pthread.h>
// 可以自行添加需要的头文件

typedef long long ll;

typedef struct {
    int* arr;
    int gnt;
    int gnt2;
    int pt;
    int width;
    int start_block;
    int end_block;
    int len;
    int wn;
} ThreadData; // 创造线程参数的结构

pthread_barrier_t barrier_start;
pthread_barrier_t barrier_done; //设置屏障
ThreadData global_tasks[8]; // 设置全局变量

int qpow(int a, int b, int p) {
    int res = 1;
    while (b) {
        if (b & 1) res = res * a % p;
        a = a * a % p;
        b >>= 1;
    }
    return res;
}

void* multiply_rows(void* arg) {
    int tid = *(int*)arg;
    while (1)
    {

        pthread_barrier_wait(&barrier_start);

        ThreadData data = global_tasks[tid];
        for (int block = data.start_block; block < data.end_block; block++) {
            int i = block * data.len;  // 块的起始索引
            int g = 1;
            for (int j = 0; j < data.width; j++) {
                int u = data.arr[i + j];
                int v = 1LL * g * data.arr[i + j + data.width] % data.pt;
                data.arr[i + j] = (u + v) % data.pt;
                data.arr[i + j + data.width] = (u - v + data.pt) % data.pt;
                g = 1LL * g * data.gnt % data.pt;
            }
        }

        pthread_barrier_wait(&barrier_done);
    }
    return NULL;
}

void* multiply_radix4(void* arg) {
    int tid = *(int*)arg;
    while (1)
    {

        pthread_barrier_wait(&barrier_start);

        ThreadData data = global_tasks[tid];
        for (int block = data.start_block; block < data.end_block; block++) {
            int i = block * data.len;  // 块的起始索引
            int w1 = 1, w2 = 1;
            int wn2 = data.wn;
            for (int j = 0; j < data.width; j++) {
                int a0 = data.arr[i + j];
                int a1 = data.arr[i + j + data.len / 4];
                int a2 = data.arr[i + j + data.len / 2];
                int a3 = data.arr[i + j + 3 * data.len / 4];
                int t0 = (a0 + 1LL * w1 * a1 % data.pt + 1LL * w2 * (a2 + 1LL * w1 * a3 % data.pt) % data.pt) % data.pt;
                int t2 = (a0 + 1LL * w1 * a1 % data.pt - 1LL * w2 * (a2 + 1LL * w1 * a3 % data.pt) % data.pt) % data.pt;
                int t1 = ((a0 - 1LL * w1 * a1 % data.pt + data.pt) % data.pt + 1LL * wn2 * (a2 - 1LL * w1 * a3 % data.pt + data.pt) % data.pt) % data.pt;
                int t3 = ((a0 - 1LL * w1 * a1 % data.pt + data.pt) % data.pt - 1LL * wn2 * (a2 - 1LL * w1 * a3 % data.pt + data.pt) % data.pt) % data.pt;

                data.arr[i + j] = (t0 + data.pt) % data.pt;
                data.arr[i + j + data.len / 4] = (t1 + data.pt) % data.pt;
                data.arr[i + j + data.len / 2] = (t2 + data.pt) % data.pt;
                data.arr[i + j + 3 * data.len / 4] = (t3 + data.pt) % data.pt;

                // 更新旋转因子
                w1 = 1LL * w1 * data.gnt2 % data.pt;
                w2 = 1LL * w2 * data.gnt % data.pt;
                wn2 = 1LL * wn2 * data.gnt % data.pt;
            }
        }
        pthread_barrier_wait(&barrier_done);
    }
    return NULL;
}

uint32_t compute_R2(uint32_t m) {
    __uint128_t R2 = ((__uint128_t)1 << 64);
    return (uint32_t)(R2 % m);
}

uint64_t inverse_mod(uint64_t a, uint64_t m) {
    int64_t t = 0, nt = 1;
    int64_t r = (int64_t)m, nr = (int64_t)a;
    int64_t q, tmp;

    while (nr != 0) {
        q = r / nr;
        tmp = nr; nr = r - q * nr; r = tmp;
        tmp = nt; nt = t - q * nt; t = tmp;
    }

    if (t < 0) t += m;
    return (uint64_t)t;
}

void fRead(int* a, int* b, int* n, int* p, int input_id) {
    // 数据输入函数
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strin = str1 + str2 + ".in";
    char data_path[strin.size() + 1];
    std::copy(strin.begin(), strin.end(), data_path);
    data_path[strin.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    fin >> *n >> *p;
    for (int i = 0; i < *n; i++) {
        fin >> a[i];
    }
    for (int i = 0; i < *n; i++) {
        fin >> b[i];
    }
}

void fRead_big(ll* a, ll* b, int* n, ll* p, int input_id) {
    // 数据输入函数
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strin = str1 + str2 + ".in";
    char data_path[strin.size() + 1];
    std::copy(strin.begin(), strin.end(), data_path);
    data_path[strin.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    fin >> *n;
    *p = 263882790666241;
    // std::cout << "input_id=" << input_id << ", 实际读到的 n=" << *n << ", p=" << *p << std::endl;
    for (int i = 0; i < *n; i++) {
        fin >> a[i];
    }
    for (int i = 0; i < *n; i++) {
        fin >> b[i];
    }
}

void fCheck(int* ab, int n, int input_id) {
    // 判断多项式乘法结果是否正确
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char data_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), data_path);
    data_path[strout.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    for (int i = 0; i < n * 2 - 1; i++) {
        int x;
        fin >> x;
        if (x != ab[i]) {
            std::cout << "多项式乘法结果错误" << std::endl;
            return;
        }
    }
    std::cout << "多项式乘法结果正确" << std::endl;
    return;
}

void fCheck_big(ll* ab, int n, int input_id) {
    // 判断多项式乘法结果是否正确
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char data_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), data_path);
    data_path[strout.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    for (int i = 0; i < n * 2 - 1; i++) {
        ll x;
        fin >> x;
        if (x != ab[i]) {
            std::cout << "x=" << x << std::endl;
            std::cout << "多项式乘法结果错误" << std::endl;
            return;
        }
    }
    std::cout << "多项式乘法结果正确" << std::endl;
    return;
}

void fWrite(int* ab, int n, int input_id) {
    // 数据输出函数, 可以用来输出最终结果, 也可用于调试时输出中间数组
    std::string str1 = "files/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()] = '\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++) {
        fout << ab[i] << '\n';
    }
}

void fWrite_big(ll* ab, int n, int input_id) {
    // 数据输出函数, 可以用来输出最终结果, 也可用于调试时输出中间数组
    std::string str1 = "files/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()] = '\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++) {
        fout << ab[i] << '\n';
    }
}

void poly_multiply(int* a, int* b, int* ab, int n, int p) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            ab[i + j] = (1LL * a[i] * b[j] % p + ab[i + j]) % p;
        }
    }
}

// 朴素NTT迭代算法的实现

int calculatef(int d, int z, int p) // 计算指数的函数
{
    int result = 1;
    while (z)
    {
        if (z & 1)
        {
            result = 1LL * result * d % p;
        }
        d = 1LL * d * d % p;
        z >>= 1;
    }
    return result;
}

int get_gn(int len, int p) // 获取g_n的函数
{
    int root = calculatef(3, (p - 1) / len, p);
    return root;
}

int reverse_get_gn(int len, int p) // 获取逆变换的值
{
    return calculatef(get_gn(len, p), p - 2, p); // 费马小定理可得
}

void reverse_index(int* a, int n) { // 实现关于索引二进制的倒转，从而方便迭代实现蝴蝶变换
    int j = 0;
    for (int i = 1; i < n - 1; i++)
    {
        int bin = n >> 1; // 类似倒着的二进制加法
        while (j >= bin)
        {
            j -= bin;
            bin >>= 1;
        }
        j += bin;
        if (i < j) // 只进行一次转换
        {
            std::swap(a[i], a[j]);
        }
    }
}

// 找到大于等于 n 的最小 2 的幂
int get_len(int n) {
    int len = 1;
    while (len < n) len <<= 1;
    return len;
}

void butterfly_multiply(int* a, int n, int p, bool isreverse) { // 蝴蝶变换迭代算法
    reverse_index(a, n);

    for (int len = 2; len <= n; len = len * 2)
    {
        int g_n = get_gn(len, p);
        if (isreverse)
        {
            g_n = reverse_get_gn(len, p);
        }
        for (int i = 0; i < n; i += len)
        {
            int g = 1;
            for (int j = 0; j < len / 2; j++)
            {
                int u = a[i + j];
                int v = 1LL * g * a[i + j + len / 2] % p;
                a[i + j] = (u + v) % p;
                a[i + j + len / 2] = (u - v + p) % p;
                g = 1LL * g * g_n % p;
            }
        }
    }
}
void ntt_multiply(int* a, int* b, int* ab, int n, int p) {
    int len = get_len(2 * n - 1);
    static int na[300000], nb[300000];
    for (int i = 0; i < n; i++)
    {
        na[i] = a[i] % p;
        nb[i] = b[i] % p;
    }
    for (int i = n; i < len; i++)
    {
        na[i] = 0;
        nb[i] = 0;
    }

    butterfly_multiply(na, len, p, false);
    butterfly_multiply(nb, len, p, false);

    for (int i = 0; i < len; i++)
    {
        na[i] = 1LL * na[i] * nb[i] % p;
    }

    butterfly_multiply(na, len, p, true);

    int inv = calculatef(len, p - 2, p);
    for (int i = 0; i < 2 * n - 1; ++i) {
        ab[i] = 1LL * na[i] * inv % p;
    }
}

void butterfly_multiply_openmp(int* a, int n, int p, bool isreverse) {
    reverse_index(a, n);

#pragma omp parallel
    {
        for (int len = 2; len <= n; len = len * 2)
        {
            int g_n = get_gn(len, p);
            if (isreverse)
            {
                g_n = reverse_get_gn(len, p);
            }

#pragma omp for schedule(static) nowait
            for (int i = 0; i < n; i += len)
            {
                int g = 1;
                for (int j = 0; j < len / 2; j++)
                {
                    int u = a[i + j];
                    int v = 1LL * g * a[i + j + len / 2] % p;
                    a[i + j] = (u + v) % p;
                    a[i + j + len / 2] = (u - v + p) % p;
                    g = 1LL * g * g_n % p;
                }
            }
        }
    }
}

void ntt_multiply_openmp(int* a, int* b, int* ab, int n, int p) {
    int len = get_len(2 * n - 1);
    static int na[300000], nb[300000];
    for (int i = 0; i < n; i++)
    {
        na[i] = a[i] % p;
        nb[i] = b[i] % p;
    }
    for (int i = n; i < len; i++)
    {
        na[i] = 0;
        nb[i] = 0;
    }

    butterfly_multiply_openmp(na, len, p, false);
    butterfly_multiply_openmp(nb, len, p, false);

    for (int i = 0; i < len; i++)
    {
        na[i] = 1LL * na[i] * nb[i] % p;
    }

    butterfly_multiply_openmp(na, len, p, true);

    int inv = calculatef(len, p - 2, p);
    for (int i = 0; i < 2 * n - 1; ++i) {
        ab[i] = 1LL * na[i] * inv % p;
    }
}

void ntt_radix4(int* a, int n, int p, bool invert) {
    reverse_index(a, n);
    int len = 4;
    for (; len <= n; len = len * 4)
    {
        int g_n = get_gn(len, p);
        int g_n0 = get_gn(len / 2, p);
        if (invert) g_n = reverse_get_gn(len, p);
        if (invert) g_n0 = reverse_get_gn(len / 2, p);
        int wn = 1;
        for (int k = 0; k < len / 4; k++)
        {
            wn = 1LL * wn * g_n % p;
        }
        for (int i = 0; i < n; i += len) {
            int w1 = 1, w2 = 1;
            int wn2 = wn;
            for (int j = 0; j < len / 4; j++) {
                int a0 = a[i + j];
                int a1 = a[i + j + len / 4];
                int a2 = a[i + j + len / 2];
                int a3 = a[i + j + 3 * len / 4];
                int t0 = (a0 + 1LL * w1 * a1 % p + 1LL * w2 * (a2 + 1LL * w1 * a3 % p) % p) % p;
                int t2 = (a0 + 1LL * w1 * a1 % p - 1LL * w2 * (a2 + 1LL * w1 * a3 % p) % p) % p;
                int t1 = ((a0 - 1LL * w1 * a1 % p + p) % p + 1LL * wn2 * (a2 - 1LL * w1 * a3 % p + p) % p) % p;
                int t3 = ((a0 - 1LL * w1 * a1 % p + p) % p - 1LL * wn2 * (a2 - 1LL * w1 * a3 % p + p) % p) % p;

                a[i + j] = (t0 + p) % p;
                a[i + j + len / 4] = (t1 + p) % p;
                a[i + j + len / 2] = (t2 + p) % p;
                a[i + j + 3 * len / 4] = (t3 + p) % p;

                // 更新旋转因子
                w1 = 1LL * w1 * g_n0 % p;
                w2 = 1LL * w2 * g_n % p;
                wn2 = 1LL * wn2 * g_n % p;
            }
        }
    }
    if (len == n * 2)
    {
        int g_n = get_gn(len / 2, p);
        if (invert)
        {
            g_n = reverse_get_gn(len / 2, p);
        }
        int g = 1;
        for (int j = 0; j < len / 4; j++)
        {
            int u = a[j];
            int v = 1LL * g * a[j + len / 4] % p;
            a[j] = (u + v) % p;
            a[j + len / 4] = (u - v + p) % p;
            g = 1LL * g * g_n % p;
        }
    }
}
void ntt_multiply_radix4(int* a, int* b, int* ab, int n, int p) {
    int len = get_len(2 * n - 1);
    static int na[300000], nb[300000];
    for (int i = 0; i < n; i++)
    {
        na[i] = a[i] % p;
        nb[i] = b[i] % p;
    }
    for (int i = n; i < len; i++)
    {
        na[i] = 0;
        nb[i] = 0;
    }

    ntt_radix4(na, len, p, false);
    ntt_radix4(nb, len, p, false);

    for (int i = 0; i < len; i++)
    {
        na[i] = 1LL * na[i] * nb[i] % p;
    }

    ntt_radix4(na, len, p, true);

    int inv = calculatef(len, p - 2, p);
    for (int i = 0; i < 2 * n - 1; ++i) {
        ab[i] = 1LL * na[i] * inv % p;
    }
}

void ntt_radix4_openmp(int* a, int n, int p, bool invert) {
    reverse_index(a, n);

    int len = 4;
    for (; len <= n; len = len * 4)
    {
        int g_n = get_gn(len, p);
        int g_n0 = get_gn(len / 2, p);
        if (invert) {
            g_n = reverse_get_gn(len, p);
            g_n0 = reverse_get_gn(len / 2, p);
        }

        int wn = 1;
        for (int k = 0; k < len / 4; k++) {
            wn = 1LL * wn * g_n % p;
        }

        int w1, w2, wn2;

#pragma omp parallel for schedule(static) private(w1, w2, wn2) shared(a, len, n, p, g_n, g_n0, wn)
        for (int i = 0; i < n; i += len) {
            w1 = 1;
            w2 = 1;
            wn2 = wn;

            for (int j = 0; j < len / 4; j++) {
                int a0 = a[i + j];
                int a1 = a[i + j + len / 4];
                int a2 = a[i + j + len / 2];
                int a3 = a[i + j + 3 * len / 4];

                int t0 = (a0 + 1LL * w1 * a1 % p + 1LL * w2 * (a2 + 1LL * w1 * a3 % p) % p) % p;
                int t2 = (a0 + 1LL * w1 * a1 % p - 1LL * w2 * (a2 + 1LL * w1 * a3 % p) % p) % p;
                int t1 = ((a0 - 1LL * w1 * a1 % p + p) % p + 1LL * wn2 * (a2 - 1LL * w1 * a3 % p + p) % p) % p;
                int t3 = ((a0 - 1LL * w1 * a1 % p + p) % p - 1LL * wn2 * (a2 - 1LL * w1 * a3 % p + p) % p) % p;

                a[i + j] = (t0 + p) % p;
                a[i + j + len / 4] = (t1 + p) % p;
                a[i + j + len / 2] = (t2 + p) % p;
                a[i + j + 3 * len / 4] = (t3 + p) % p;

                w1 = 1LL * w1 * g_n0 % p;
                w2 = 1LL * w2 * g_n % p;
                wn2 = 1LL * wn2 * g_n % p;
            }
        }
    }

    // 处理最后的 Radix-2 阶段
    if (len == n * 2) {
        int g_n = get_gn(len / 2, p);
        if (invert) {
            g_n = reverse_get_gn(len / 2, p);
        }

        int g;

#pragma omp parallel for schedule(static) private(g) shared(a, len, p, g_n)
        for (int j = 0; j < len / 4; j++) {
            g = 1;
            for (int k = 0; k < j; k++) {
                g = 1LL * g * g_n % p;
            }
            int u = a[j];
            int v = 1LL * g * a[j + len / 4] % p;
            a[j] = (u + v) % p;
            a[j + len / 4] = (u - v + p) % p;
        }
    }
}

void ntt_multiply_radix4_openmp(int* a, int* b, int* ab, int n, int p) {
    int len = get_len(2 * n - 1);
    static int na[300000], nb[300000];
    for (int i = 0; i < n; i++)
    {
        na[i] = a[i] % p;
        nb[i] = b[i] % p;
    }
    for (int i = n; i < len; i++)
    {
        na[i] = 0;
        nb[i] = 0;
    }

    ntt_radix4_openmp(na, len, p, false);
    ntt_radix4_openmp(nb, len, p, false);

    for (int i = 0; i < len; i++)
    {
        na[i] = 1LL * na[i] * nb[i] % p;
    }

    ntt_radix4_openmp(na, len, p, true);

    int inv = calculatef(len, p - 2, p);
    for (int i = 0; i < 2 * n - 1; ++i) {
        ab[i] = 1LL * na[i] * inv % p;
    }
}

static pthread_t threads[8];
static int initialized = 0;

void ntt_radix4_Pthread(int* a, int n, int p, bool invert) {
    reverse_index(a, n);
    int len = 4;
    ThreadData thread_data[4];
    int tp[4];
    if (!initialized) {
        pthread_barrier_init(&barrier_start, NULL, 5); // 4 workers + 1 main
        pthread_barrier_init(&barrier_done, NULL, 5);
        for (int i = 0; i < 4; i++) {
            tp[i] = i;

            pthread_create(&threads[i], NULL, multiply_radix4, &tp[i]);
        }
        initialized = 1;
    }
    for (; len <= n; len = len * 4)
    {
        int width = len / 4;
        int g_n = get_gn(len, p);
        int g_n0 = get_gn(len / 2, p);
        if (invert) g_n = reverse_get_gn(len, p);
        if (invert) g_n0 = reverse_get_gn(len / 2, p);
        int wn = 1;
        for (int k = 0; k < len / 4; k++)
        {
            wn = 1LL * wn * g_n % p;
        }
        int total_blocks = n / len;
        for (int t = 0; t < 4; t++)
        {
            int start_block = t * total_blocks / 4;
            int end_block = (t + 1) * total_blocks / 4;
            global_tasks[t].start_block = start_block;
            global_tasks[t].end_block = end_block;
            global_tasks[t].arr = a;
            global_tasks[t].gnt = g_n;
            global_tasks[t].gnt2 = g_n0;
            global_tasks[t].wn = wn;
            global_tasks[t].pt = p;
            global_tasks[t].width = width;
            global_tasks[t].len = len;
        }
        pthread_barrier_wait(&barrier_start);
        // 等待所有线程完成
        pthread_barrier_wait(&barrier_done);
    }
    if (len == n * 2)
    {
        int g_n = get_gn(len / 2, p);
        if (invert)
        {
            g_n = reverse_get_gn(len / 2, p);
        }
        int g = 1;
        for (int j = 0; j < len / 4; j++)
        {
            int u = a[j];
            int v = 1LL * g * a[j + len / 4] % p;
            a[j] = (u + v) % p;
            a[j + len / 4] = (u - v + p) % p;
            g = 1LL * g * g_n % p;
        }
    }
}

void ntt_multiply_radix4_Pthread(int* a, int* b, int* ab, int n, int p) {
    int len = get_len(2 * n - 1);
    static int na[300000], nb[300000];
    for (int i = 0; i < n; i++)
    {
        na[i] = a[i] % p;
        nb[i] = b[i] % p;
    }
    for (int i = n; i < len; i++)
    {
        na[i] = 0;
        nb[i] = 0;
    }

    ntt_radix4_Pthread(na, len, p, false);
    ntt_radix4_Pthread(nb, len, p, false);

    for (int i = 0; i < len; i++)
    {
        na[i] = 1LL * na[i] * nb[i] % p;
    }

    ntt_radix4_Pthread(na, len, p, true);

    int inv = calculatef(len, p - 2, p);
    for (int i = 0; i < 2 * n - 1; ++i) {
        ab[i] = 1LL * na[i] * inv % p;
    }
}

void butterfly_multiply_Pthread(int* a, int n, int p, bool isreverse) { // 蝴蝶变换迭代算法
    reverse_index(a, n);
    ThreadData thread_data[4];
    int tp[4];
    if (!initialized) {
        pthread_barrier_init(&barrier_start, NULL, 5); // 4 workers + 1 main
        pthread_barrier_init(&barrier_done, NULL, 5);
        for (int i = 0; i < 4; i++) {
            tp[i] = i;

            pthread_create(&threads[i], NULL, multiply_rows, &tp[i]);
        }
        initialized = 1;
    }

    for (int len = 2; len <= n; len = len * 2)
    {
        int g_n = get_gn(len, p);
        int width = len / 2;
        if (isreverse)
        {
            g_n = reverse_get_gn(len, p);
        }
        int total_blocks = n / len; // 按照块来进行处理，去除原先的行索引
        for (int t = 0; t < 4; t++)
        {
            int start_block = t * total_blocks / 4;
            int end_block = (t + 1) * total_blocks / 4;
            global_tasks[t].start_block = start_block;
            global_tasks[t].end_block = end_block;
            global_tasks[t].arr = a;
            global_tasks[t].gnt = g_n;
            global_tasks[t].pt = p;
            global_tasks[t].width = width;
            global_tasks[t].len = len;
        }
        pthread_barrier_wait(&barrier_start);
        // 等待所有线程完成
        pthread_barrier_wait(&barrier_done);
    }
}

void ntt_multiply_Pthread(int* a, int* b, int* ab, int n, int p) {
    int len = get_len(2 * n - 1);
    static int na[300000], nb[300000];
    for (int i = 0; i < n; i++)
    {
        na[i] = a[i] % p;
        nb[i] = b[i] % p;
    }
    for (int i = n; i < len; i++)
    {
        na[i] = 0;
        nb[i] = 0;
    }

    butterfly_multiply_Pthread(na, len, p, false);
    butterfly_multiply_Pthread(nb, len, p, false);

    for (int i = 0; i < len; i++)
    {
        na[i] = 1LL * na[i] * nb[i] % p;
    }

    butterfly_multiply_Pthread(na, len, p, true);

    int inv = calculatef(len, p - 2, p);
    for (int i = 0; i < 2 * n - 1; ++i) {
        ab[i] = 1LL * na[i] * inv % p;
    }
}

void butterfly_multiply_Pthread2(int* a, int n, int p, bool isreverse) { // 蝴蝶变换迭代算法
    reverse_index(a, n);
    ThreadData thread_data[2];
    int tp[2];
    if (!initialized) {
        pthread_barrier_init(&barrier_start, NULL, 3); // 4 workers + 1 main
        pthread_barrier_init(&barrier_done, NULL, 3);
        for (int i = 0; i < 2; i++) {
            tp[i] = i;

            pthread_create(&threads[i], NULL, multiply_rows, &tp[i]);
        }
        initialized = 1;
    }

    for (int len = 2; len <= n; len = len * 2)
    {
        int g_n = get_gn(len, p);
        int width = len / 2;
        if (isreverse)
        {
            g_n = reverse_get_gn(len, p);
        }
        int total_blocks = n / len; // 按照块来进行处理，去除原先的行索引
        for (int t = 0; t < 2; t++)
        {
            int start_block = t * total_blocks / 2;
            int end_block = (t + 1) * total_blocks / 2;
            global_tasks[t].start_block = start_block;
            global_tasks[t].end_block = end_block;
            global_tasks[t].arr = a;
            global_tasks[t].gnt = g_n;
            global_tasks[t].pt = p;
            global_tasks[t].width = width;
            global_tasks[t].len = len;
        }
        pthread_barrier_wait(&barrier_start);
        // 等待所有线程完成
        pthread_barrier_wait(&barrier_done);
    }
}

void ntt_multiply_Pthread2(int* a, int* b, int* ab, int n, int p) {
    int len = get_len(2 * n - 1);
    static int na[300000], nb[300000];
    for (int i = 0; i < n; i++)
    {
        na[i] = a[i] % p;
        nb[i] = b[i] % p;
    }
    for (int i = n; i < len; i++)
    {
        na[i] = 0;
        nb[i] = 0;
    }

    butterfly_multiply_Pthread2(na, len, p, false);
    butterfly_multiply_Pthread2(nb, len, p, false);

    for (int i = 0; i < len; i++)
    {
        na[i] = 1LL * na[i] * nb[i] % p;
    }

    butterfly_multiply_Pthread2(na, len, p, true);

    int inv = calculatef(len, p - 2, p);
    for (int i = 0; i < 2 * n - 1; ++i) {
        ab[i] = 1LL * na[i] * inv % p;
    }
}


const int MOD1 = 469762049;
const int MOD2 = 998244353;
const int MOD3 = 1004535809;

void reverse_index_big(ll* a, int n) { // 实现关于索引二进制的倒转，从而方便迭代实现蝴蝶变换
    int j = 0;
    for (int i = 1; i < n - 1; i++)
    {
        int bin = n >> 1; // 类似倒着的二进制加法
        while (j >= bin)
        {
            j -= bin;
            bin >>= 1;
        }
        j += bin;
        if (i < j) // 只进行一次转换
        {
            std::swap(a[i], a[j]);
        }
    }
}

void butterfly_multiply_big(ll* a, int n, int p, bool isreverse) { // 蝴蝶变换迭代算法
    reverse_index_big(a, n);

    for (int len = 2; len <= n; len = len * 2)
    {
        int g_n = get_gn(len, p);
        if (isreverse)
        {
            g_n = reverse_get_gn(len, p);
        }
        for (int i = 0; i < n; i += len)
        {
            int g = 1;
            for (int j = 0; j < len / 2; j++)
            {
                ll u = a[i + j];
                ll v = 1LL * g * a[i + j + len / 2] % p;
                a[i + j] = (u + v) % p;
                a[i + j + len / 2] = (u - v + p) % p;
                g = 1LL * g * g_n % p;
            }
        }
    }
}
void ntt_multiply_big(ll* a, ll* b, ll* ab, int n, int p) {
    int len = get_len(2 * n - 1);
    static ll na[300000], nb[300000];
    for (int i = 0; i < n; i++)
    {
        na[i] = a[i] % p;
        nb[i] = b[i] % p;
    }
    for (int i = n; i < len; i++)
    {
        na[i] = 0;
        nb[i] = 0;
    }

    butterfly_multiply_big(na, len, p, false);
    butterfly_multiply_big(nb, len, p, false);

    for (int i = 0; i < len; i++)
    {
        na[i] = 1LL * na[i] * nb[i] % p;
    }

    butterfly_multiply_big(na, len, p, true);

    int inv = calculatef(len, p - 2, p);
    for (int i = 0; i < 2 * n - 1; ++i) {
        ab[i] = 1LL * na[i] * inv % p;
    }
}

void combine_crt(ll* a, ll* b, ll* ab, int n, ll p)
{
    ll* a1 = new ll[2 * n - 1];
    ll* a2 = new ll[2 * n - 1];
    ll* a3 = new ll[2 * n - 1];
    ll* result = new ll[2 * n - 1];
    ntt_multiply_big(a, b, a1, n, MOD1);
    ntt_multiply_big(a, b, a2, n, MOD2);
    ntt_multiply_big(a, b, a3, n, MOD3);
    const ll M1 = MOD2 * (ll)MOD3;
    const ll M2 = MOD1 * (ll)MOD3;
    const ll M3 = MOD1 * (ll)MOD2;
    const ll M = (ll)MOD1 * MOD2 % ((ll)MOD1 * MOD2);
    ll inv1 = qpow(MOD1 % MOD2, MOD2 - 2, MOD2);
    ll inv2 = qpow(MOD2 % MOD1, MOD1 - 2, MOD1);
    ll inv3 = qpow(MOD3 % MOD1, MOD1 - 2, MOD1);

    for (int i = 0; i < 2 * n - 1; i++) {
        ll x = (ll)(a2[i] - a1[i] + MOD2) % MOD2;
        x = (x * inv1) % MOD2;
        ll t = a1[i] + x * MOD1;

        ll y = (ll)(a3[i] - t % MOD3 + MOD3) % MOD3;
        ll inv_M12 = qpow(M1 % MOD3, MOD3 - 2, MOD3);
        y = y * inv_M12 % MOD3;

        result[i] = t + y * M1;
    }
    if (p > 0) {
        for (int i = 0; i < 2 * n - 1; i++) {
            ab[i] = result[i] % p;
        }
    }
}

uint32_t inverse(uint32_t n, uint64_t mod) {
    // 对于 mod = 2^32，直接当作 uint64_t 处理
    int64_t t = 0, nt = 1;
    int64_t r = (int64_t)mod, nr = (int64_t)n;
    int64_t q, tmp;

    while (nr != 0) {
        q = r / nr;
        tmp = nr; nr = r - q * nr; r = tmp;
        tmp = nt; nt = t - q * nt; t = tmp;
    }

    // r 应该是 1（因为 n 是奇数，与 2^32 互质）
    if (t < 0) t += mod;
    return (uint32_t)t;
}


uint32x2_t reduce(uint64x2_t a, int p) {
    uint32_t n = uint32_t(p);
    uint32_t nr = inverse(p, 1ULL << 32);
    // std::cout<<nr<<std::endl;
    uint32x2_t v_n = vdup_n_u32(n);
    uint32x2_t v_nr = vdup_n_u32(nr);
    uint32x2_t low = vmovn_u64(a);
    uint32x2_t q = vmul_u32(low, v_nr);
    uint64x2_t m = vmull_u32(q, v_n);
    uint64x2_t yt = vsubq_u64(a, m);

    uint32x2_t y = vmovn_u64(vshrq_n_u64(yt, 32));

    uint32_t m_vals[2];
    vst1_u32(m_vals, y);
    // printf("y[0] = %llu, y[1] = %llu\n", m_vals[0], m_vals[1]);

    uint64x2_t vcmp = vcltq_u64(a, m);
    uint32x2_t vn_low = vmovn_u64(vandq_u64(vcmp, vdupq_n_u64(n)));
    uint32x2_t vresult = vadd_u32(y, vn_low);
    return vresult;
}


uint32x2_t reduce_mod(uint64x2_t t, int p) {

    uint64_t t_arr[2];
    vst1q_u64(t_arr, t);

    uint32_t result_arr[2];
    for (int i = 0; i < 2; i++) {
        result_arr[i] = t_arr[i] % p;
    }

    uint32x2_t result = vld1_u32(result_arr);

    // uint32_t R_mod = (1ULL << 32) % p;
    // uint32x2_t t_low = vmovn_u64(t);
    // uint32x2_t t_high = vshrn_n_u64(t, 32);
    // uint32x2_t r = vdup_n_u32(R_mod);

    // // 计算各部分
    // uint64x2_t low_part = vmull_u32(t_low, r);     // a_low * b
    // uint64x2_t high_part = vmull_u32(t_high, r);   // a_high * b

    // // 组合：low_part + (high_part << 32)
    // uint64x2_t use = vaddq_u64(low_part, vshlq_n_u64(high_part, 32));

    // return reduce(use,p);
    return result;
}

uint32x4_t reduce_2(uint32x4_t a, int p) {
    uint32x2_t v1 = vget_low_u32(a);
    uint32x2_t v2 = vget_high_u32(a);

    // 计算 R2 = R² mod m
    uint32_t R2 = compute_R2(p);

    uint32x2_t v_R2 = vdup_n_u32(R2);  // 将 R2 复制到两个通道
    uint64x2_t multiplied1 = vmull_u32(v1, v_R2);  // 向量乘法
    uint64x2_t multiplied2 = vmull_u32(v2, v_R2);
    uint32x2_t r1 = reduce(multiplied1, p);
    uint32x2_t r2 = reduce(multiplied2, p);
    uint32x4_t r = vcombine_u32(r1, r2);
    return r;
}

uint32x4_t montgomery_mult(uint32x4_t a, uint32x4_t b, int p) // Montgomery模乘化
{
    uint32_t n = uint32_t(p);
    uint32_t nr = inverse(p, 1ULL << 32);
    uint32x4_t ar = reduce_2(a, p);
    uint32x4_t br = reduce_2(b, p);
    uint32x2_t v_n = vdup_n_u32(n);
    uint32x2_t v_nr = vdup_n_u32(nr);
    uint32x2_t ar1 = vget_low_u32(ar);
    uint32x2_t ar2 = vget_high_u32(ar);
    uint32x2_t br1 = vget_low_u32(br);
    uint32x2_t br2 = vget_high_u32(br);
    uint64x2_t t1 = vmull_u32(ar1, br1);
    uint64x2_t t2 = vmull_u32(ar2, br2);

    uint64_t r32 = 1ULL << 32;
    uint32x2_t t1_mod = reduce_mod(t1, p);
    uint32x2_t t2_mod = reduce_mod(t2, p);
    uint64_t r_inv = inverse_mod(r32 % n, (uint64_t)n);
    uint32_t r_inv32 = (uint32_t)r_inv;  // 转换为 32 位
    uint32x2_t v_r_inv = vdup_n_u32(r_inv32);  // [r_inv, r_inv]
    uint64x2_t sum1 = vmull_u32(t1_mod, v_r_inv);
    uint64x2_t sum2 = vmull_u32(t2_mod, v_r_inv);
    uint32x2_t sum1_mod = reduce_mod(sum1, p);
    uint32x2_t sum2_mod = reduce_mod(sum2, p);
    uint32x4_t r = vcombine_u32(sum1_mod, sum2_mod);

    uint32x2_t vv1 = vget_low_u32(r);
    uint32x2_t vv2 = vget_high_u32(r);
    uint64x2_t v1_64 = vmovl_u32(vv1);  // 低2个32位 -> 2个64位
    uint64x2_t v2_64 = vmovl_u32(vv2);  // 高2个32位 -> 2个64位


    uint32x2_t rr1 = reduce(v1_64, p);
    uint32x2_t rr2 = reduce(v2_64, p);
    uint32x4_t r_original = vcombine_u32(rr1, rr2);

    return r_original;
}


uint32x4_t Montgomery_add(uint32x4_t a, uint32x4_t b, int p) {
    uint32x4_t addnum = vaddq_u32(a, b);
    uint32x4_t vp = vdupq_n_u32(p);
    uint32x4_t v_ge = vcgeq_u32(addnum, vp);
    addnum = vsubq_u32(addnum, vandq_u32(v_ge, vp));
    return addnum;
}

uint32x4_t Montgomery_sub(uint32x4_t a, uint32x4_t b, int p) {
    uint32x4_t sub = vsubq_u32(a, b);
    uint32x4_t lt = vcltq_u32(a, b);    // a < b 才加 p
    uint32x4_t vp = vdupq_n_u32(p);
    sub = vaddq_u32(sub, vandq_u32(lt, vp));
    return sub;
}

uint32x4_t vec_mul_mod(uint32x4_t a, uint32x4_t b, int p) {
    uint32_t a_arr[4], b_arr[4], r_arr[4];
    vst1q_u32(a_arr, a);
    vst1q_u32(b_arr, b);
    for (int i = 0; i < 4; i++) {
        r_arr[i] = (uint64_t)a_arr[i] * b_arr[i] % p;
    }
    return vld1q_u32(r_arr);
}

void butterfly_neon_multiply(int* a, int n, int p, bool isreverse) { // 蝴蝶变换迭代算法
    reverse_index(a, n);

    for (int len = 2; len <= n; len = len * 2)
    {
        int g_n = get_gn(len, p);
        if (isreverse)
        {
            g_n = reverse_get_gn(len, p);
        }

        int half_len = len / 2;
        int* gn_powers = new int[half_len];  // 预分配数组
        gn_powers[0] = 1;
        for (int j = 1; j < half_len; j++) {
            gn_powers[j] = 1LL * gn_powers[j - 1] * g_n % p;
        }

        for (int i = 0; i < n; i += len)
        {
            int g = 1;
            for (int j = 0; j < len / 2; j += 4)
            {
                if (j + 4 >= len / 2)
                {
                    for (int k = j; k < len / 2; k++)
                    {
                        int u = a[i + k];
                        int v = 1LL * gn_powers[k] * a[i + k + len / 2] % p;
                        a[i + k] = (u + v) % p;
                        a[i + k + len / 2] = (u - v + p) % p;
                    }
                }
                else
                {
                    uint32x4_t ui = vld1q_u32((const uint32_t*)(a + i + j));
                    uint32x4_t vg = vld1q_u32((uint32_t*)(gn_powers + j));
                    uint32x4_t vv = vld1q_u32((const uint32_t*)(a + i + j + len / 2));
                    uint32x4_t vi = montgomery_mult(vg, vv, p);

                    uint32x4_t uav = Montgomery_add(ui, vi, p);
                    uint32x4_t usv = Montgomery_sub(ui, vi, p);
                    vst1q_u32((uint32_t*)(&a[i + j]), uav);
                    vst1q_u32((uint32_t*)(&a[i + j + len / 2]), usv);
                }
            }
        }
    }
}


void ntt_neon_multiply(int* a, int* b, int* ab, int n, int p) {
    int len = get_len(2 * n - 1);
    static int na[300000], nb[300000];
    for (int i = 0; i < n; i++)
    {
        na[i] = a[i] % p;
        nb[i] = b[i] % p;
    }
    for (int i = n; i < len; i++)
    {
        na[i] = 0;
        nb[i] = 0;
    }
    butterfly_neon_multiply(na, len, p, false);
    butterfly_neon_multiply(nb, len, p, false);

    for (int i = 0; i < len; i += 4)
    {
        if (i + 4 >= len)
        {
            for (int k = i; k < len; k++)
            {
                na[k] = 1LL * na[k] * nb[k] % p;
            }
        }
        else
        {
            uint32x4_t a = vld1q_u32((const uint32_t*)(na + i));
            uint32x4_t b = vld1q_u32((const uint32_t*)(nb + i));
            uint32x4_t r = montgomery_mult(a, b, p);
            vst1q_u32((uint32_t*)(&na[i]), r);
        }
    }

    butterfly_neon_multiply(na, len, p, true);

    int inv = calculatef(len, p - 2, p);
    for (int i = 0; i < 2 * n - 1; i += 4) {
        if (i + 4 >= 2 * n - 1)
        {
            for (int k = i; k < 2 * n - 1; k++)
            {
                ab[k] = 1LL * na[k] * inv % p;
            }

        }
        else
        {
            uint32x4_t vinv = vdupq_n_u32(inv);
            uint32x4_t vna = vld1q_u32((const uint32_t*)(na + i));
            uint32x4_t vab = montgomery_mult(vna, vinv, p);
            vst1q_u32((uint32_t*)(&ab[i]), vab);
        }
    }
}


int a[300000], b[300000], ab[300000];
ll a_big[300000], b_big[300000], ab_big[300000];
int main(int argc, char* argv[])
{

    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1 的形式
    // 输入模数分别为 7340033 104857601 469762049 263882790666241
    // 第四个模数超过了整型表示范围, 如果实现此模数意义下的多项式乘法需要修改框架
    // 对第四个模数的输入数据不做必要要求, 如果要自行探索大模数 NTT, 请在完成前三个模数的基础代码及优化后实现大模数 NTT
    // 输入文件共五个, 第一个输入文件 n = 4, 其余四个文件分别对应四个模数, n = 131072
    // 在实现快速数论变化前, 后四个测试样例运行时间较久, 推荐调试正确性时只使用输入文件 1
    int test_begin = 0;
    int test_end = 3;
    for (int i = test_begin; i <= test_end; ++i) {
        // if(i==4)
        // {
        //     long double ans = 0;
        //     int n_; 
        //     ll p_n;
        //     fRead_big(a_big, b_big, &n_, &p_n, i);
        //     memset(ab_big, 0, sizeof(ab_big));
        //     auto Start = std::chrono::high_resolution_clock::now();
        //     std::cout<<a_big[0]<<" "<<b_big[0]<<std::endl;
        //     combine_crt(a_big, b_big, ab_big, n_, p_n);
        //     std::cout<<ab_big[0]<<std::endl;
        //     auto End = std::chrono::high_resolution_clock::now();
        //     std::chrono::duration<double, std::ratio<1, 1000>>elapsed = End - Start;
        //     ans += elapsed.count();
        //     fCheck_big(ab_big, n_, i);
        //     std::cout << "average latency for n = " << n_ << " p = " << p_n << " : " << ans << " (us) " << std::endl;
        //     // 可以使用 fWrite 函数将 ab 的输出结果打印到 files 文件夹下
        //     // 禁止使用 cout 一次性输出大量文件内容
        //     fWrite_big(ab_big, n_, i);
        // }
        // else
        // {
        long double ans = 0;
        int n_, p_;
        fRead(a, b, &n_, &p_, i);
        memset(ab, 0, sizeof(ab));
        auto Start = std::chrono::high_resolution_clock::now();
        // TODO : 将 poly_multiply 函数替换成你写的 ntt
        // poly_multiply(a, b, ab, n_, p_);
        // ntt_neon_multiply(a, b, ab, n_, p_);
        // ntt_multiply(a, b, ab, n_, p_);
        // ntt_multiply_Pthread(a, b, ab, n_, p_);
        // ntt_multiply_Pthread2(a, b, ab, n_, p_);
        // ntt_multiply_radix4(a, b, ab, n_, p_);
        ntt_multiply_radix4_Pthread(a, b, ab, n_, p_);
        // ntt_multiply_openmp(a, b, ab, n_, p_);
        // ntt_multiply_radix4_openmp(a, b, ab, n_, p_);
        auto End = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::ratio<1, 1000>>elapsed = End - Start;
        ans += elapsed.count();
        fCheck(ab, n_, i);
        std::cout << "average latency for n = " << n_ << " p = " << p_ << " : " << ans << " (us) " << std::endl;
        // 可以使用 fWrite 函数将 ab 的输出结果打印到 files 文件夹下
        // 禁止使用 cout 一次性输出大量文件内容
        fWrite(ab, n_, i);
        // }
    }

    return 0;

}