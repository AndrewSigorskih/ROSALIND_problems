#include <iostream>
#include <fstream>
#include <unordered_map>

const char* INFILENAME = "rosalind_2sum.txt";
const char* OUTFILENAME = "out.txt";
const int zero = 0; // sign flip

void twoSum(int* arr, int n, std::ofstream& ost)
{
    std::unordered_map<int, int> map;
    for (int i = 0; i < n; ++i)
    {
        if (map.find(arr[i]) != map.end())
        {
            ost << (map[arr[i]] + 1) << " " << (i + 1) << '\n';
            return;
        } else {
            map[zero - arr[i]] = i;
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
        twoSum(arr, n, ost);
    }
    delete[] arr;
    return 0;
}