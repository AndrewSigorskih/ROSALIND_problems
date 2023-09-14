#include <iostream>
#include <fstream>
#include <utility>
#include <unordered_map>

const char* INFILENAME = "rosalind_3sum.txt";
const char* OUTFILENAME = "out.txt";
using ipair = std::pair<int, int>;
using std::make_pair;
const int zero = 0; //sign flip

ipair twoSum(int* arr, int n, int target)
{
    std::unordered_map<int, int> map;
    for (int i = 0; i < n; ++i)
    {
        if (map.find(arr[i]) != map.end())
        {
            return make_pair(map[arr[i]] + 1, i + 1);
        } else {
            map[target - arr[i]] = i;
        }
    }
    return make_pair(-1, -1);
}

void threeSum(int* arr, int n, std::ofstream& ost)
{
    for (int i = 0; i < n; ++i)
    {
        ipair res = twoSum(arr+i+1, n-i-1, zero-arr[i]);
        if ((res.first != -1) && (res.second != -1))
        {
            ost << (i + 1) << " " << (res.first + i + 1) << " " << (res.second + i + 1) << '\n';
            return;
        }
    }
    ost << -1 << '\n';
} 

int main()
{
    std::ifstream ist{INFILENAME};
    if (!ist) 
    {
        std::cout << "Cannot open input file!\n";
        exit(1);
    }
    std::ofstream ost{OUTFILENAME};
    if (!ost) 
    {
        std::cout << "Cannot open output file!\n";
        exit(1);
    }

    int n, k;
    ist >> k >> n;
    int* arr = new int[n];
    for(int i = 0; i < k; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            ist >> arr[j];
        }
        threeSum(arr, n, ost);
    }
    delete[] arr;
    return 0;
}