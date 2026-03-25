#include<iostream>
#include<chrono>  // 替代 Windows.h 的计时功能
// #include<Windows.h>  // 删除这行
using namespace std;
using namespace std::chrono;  // 使用 chrono 命名空间

int *normal(int n,int** A,int* b)
{
    int *sum=new int[n];
    for(int i=0;i<n;i++)
    {
        sum[i]=0;
        for(int j=0;j<n;j++)
        {
            sum[i]+=A[j][i]*b[j];
        }
    }
    return sum;
}

int *cache1(int n,int** A,int* b)
{
    int *sum=new int[n];
    for(int c=0;c<n;c+=4)
    {
        int li=(c+4<n)?c+4:n;
        for(int i=c;i<li;i++)
        {
            sum[i]=0;
        }
        for(int i=0;i<n;i++)
        {
            int bi=b[i];
            for(int j=c;j<li;j++)
            {
                sum[j]+=A[i][j]*bi;
            }
        }
    }
    return sum;
}

int *cache2(int n,int** A,int* b)
{
    int *sum=new int[n];
    for(int c=0;c<n;c+=8)
    {
        int li=(c+8<n)?c+8:n;
        for(int i=c;i<li;i++)
        {
            sum[i]=0;
        }
        for(int i=0;i<n;i++)
        {
            int bi=b[i];
            for(int j=c;j<li;j++)
            {
                sum[j]+=A[i][j]*bi;
            }
        }
    }
    return sum;
}

// 修改后的计时函数 - 使用 chrono
void printmatrixcache1(int *b,int n,int** A)
{
    auto start = high_resolution_clock::now();  // 开始计时
    int *sum = cache1(n,A,b);
    auto end = high_resolution_clock::now();    // 结束计时
    
    auto duration = duration_cast<microseconds>(end - start);
    cout << "cache1(4) Col: " << duration.count() / 1000.0 << " ms" << endl;
    
    delete[] sum;  // 记得释放内存
}

void printmatrixcache2(int *b,int n,int** A)
{
    auto start = high_resolution_clock::now();
    int *sum = cache2(n,A,b);
    auto end = high_resolution_clock::now();
    
    auto duration = duration_cast<microseconds>(end - start);
    cout << "cache2(8) Col: " << duration.count() / 1000.0 << " ms" << endl;
    
    delete[] sum;
}

void printmatrixnormal(int *b,int n,int** A)
{
    auto start = high_resolution_clock::now();
    int *sum = normal(n,A,b);
    auto end = high_resolution_clock::now();
    
    auto duration = duration_cast<microseconds>(end - start);
    cout << "normal Col: " << duration.count() / 1000.0 << " ms" << endl;
    
    delete[] sum;
}

int main()
{
    int n;
    cin>>n;
    int **A=new int*[n];
    int *b=new int[n];
    for(int i=0;i<n;i++)
    {
        int *t= new int[n];
        for(int j=0;j<n;j++)
        {
            t[j]=j;
            b[j]=j;
        }
        A[i]=t;
    }
    printmatrixnormal(b,n,A);
    printmatrixcache1(b,n,A);
    printmatrixcache2(b,n,A);
    
    // 释放内存
    for(int i=0;i<n;i++) {
        delete[] A[i];
    }
    delete[] A;
    delete[] b;
    
    return 0;
}