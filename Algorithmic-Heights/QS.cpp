#include <iostream>
#include <fstream>
#include <cstdlib>

const char* INFILENAME = "rosalind_qs.txt";
const char* OUTFILENAME = "out.txt";

int compare (const void* a, const void* b)
{
  return ( *(int*)a - *(int*)b );
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

    int n;
    ist >> n;
    int* arr = new int[n];

    for (int i = 0; i < n; ++i)
    {
        ist >> arr[i];
    }

    std::qsort(arr,
                n,
                sizeof(int),
                compare);

    for (int i = 0; i < n; ++i)
    {
        ost << arr[i];
        if (i < n-1)
            ost << " ";
    }
    ost << '\n';

    delete[] arr;
    return 0;
}