#include<iostream>
#include<Windows.h>
using namespace std;

int normaladd(int* arr,int n)
{
    int sum=0;
    for(int i=0;i<n;i++)
    {
        sum+=arr[i];
    }
    return sum;
}

int supadd(int *arr,int n)
{
    // int time=2;
    // while(n/time!=0)
    // {
    //     for(int i=0;i<n;i+=time)
    //     {
    //         arr[i]+=arr[i+time/2];
    //     }
    //     time=time*2;
    // }
    // return arr[0];
    int sum1=0;
    int sum2=0;
    for(int i=0;i<n;i+=2)
    {
        sum1+=arr[i];
        sum2+=arr[i+1];
    }
    return sum1+sum2;
}

int main()
{
    int n;
    cin>>n;
    int* arr1=new int[n];
    for(int i=0;i<n;i++)
    {
        arr1[i]=i;
    }
    int* arr2=new int[n];
    for(int i=0;i<n;i++)
    {
        arr2[i]=i;
    }
    long long head1, tail1 , freq1 ;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq1 );
    QueryPerformanceCounter((LARGE_INTEGER *)&head1);
    for(int i=0;i<10000;i++)
    {
        int sum1=supadd(arr1,n);
    }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail1 );
    cout << "supadd Col: " << ( tail1 - head1) * 1000.0 / freq1<<"ms"<<endl;
    long long head2, tail2 , freq2 ;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq2 );
    QueryPerformanceCounter((LARGE_INTEGER *)&head2);
    for(int i=0;i<10000;i++)
    {
        int sum2=normaladd(arr2,n);
    }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail2 );
    cout << "normaladd Col: " << ( tail2 - head2) * 1000.0 / freq2<<"ms"<<endl;
    return 0;
}