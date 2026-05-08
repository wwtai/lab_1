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
// 可以自行添加需要的头文件

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
                    uint32x4_t vi = vec_mul_mod(vg, vv, p);

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
            uint32x4_t r = vec_mul_mod(a, b, p);
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
            uint32x4_t vab = vec_mul_mod(vna, vinv, p);
            vst1q_u32((uint32_t*)(&ab[i]), vab);
        }
    }
}


int a[300000], b[300000], ab[300000];
int main(int argc, char* argv[])
{

    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1 的形式
    // 输入模数分别为 7340033 104857601 469762049 263882790666241
    // 第四个模数超过了整型表示范围, 如果实现此模数意义下的多项式乘法需要修改框架
    // 对第四个模数的输入数据不做必要要求, 如果要自行探索大模数 NTT, 请在完成前三个模数的基础代码及优化后实现大模数 NTT
    // 输入文件共五个, 第一个输入文件 n = 4, 其余四个文件分别对应四个模数, n = 131072
    // 在实现快速数论变化前, 后四个测试样例运行时间较久, 推荐调试正确性时只使用输入文件 1
    int test_begin = 0;
    int test_end = 4;
    for (int i = test_begin; i <= test_end; ++i) {
        long double ans = 0;
        int n_, p_;
        fRead(a, b, &n_, &p_, i);
        memset(ab, 0, sizeof(ab));
        auto Start = std::chrono::high_resolution_clock::now();
        // TODO : 将 poly_multiply 函数替换成你写的 ntt
        // poly_multiply(a, b, ab, n_, p_);
        // ntt_neon_multiply(a, b, ab, n_, p_);
        ntt_multiply(a, b, ab, n_, p_);
        auto End = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::ratio<1, 1000>>elapsed = End - Start;
        ans += elapsed.count();
        fCheck(ab, n_, i);
        std::cout << "average latency for n = " << n_ << " p = " << p_ << " : " << ans << " (us) " << std::endl;
        // 可以使用 fWrite 函数将 ab 的输出结果打印到 files 文件夹下
        // 禁止使用 cout 一次性输出大量文件内容
        fWrite(ab, n_, i);
    }

    return 0;

}