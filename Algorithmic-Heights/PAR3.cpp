#include <iostream>
#include <fstream>
#include <vector>
using std::vector;
using std::swap;
const char* INFILENAME = "rosalind_par3.txt";
const char* OUTFILENAME = "out.txt";

void threeWayPartition(vector<int>& arr)
{
    int pivot = arr[0];
    int left = 0, i = 0;
    int right = arr.size() - 1, j = arr.size() - 1;
    while (left <= right)
    {
        while ((left <= right) && (arr[left] <= pivot))
        {
            if (arr[left] < pivot)
            {
                swap(arr[i], arr[left]);
                i++;
            }
            left++;
        }
        swap(arr[left], arr[right]);

        while ((left <= right) && (arr[right]>= pivot))
        {
            if (arr[right] > pivot)
            {
                swap(arr[right], arr[j]);
                j--;
            }
            right--;
        }
        swap(arr[left], arr[right]);
    }
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
    vector<int> arr(n);
    for (int i = 0; i < n; ++i)
    {
        ist >> arr[i];
    }
    threeWayPartition(arr);
    for (int i = 0; i < n; ++i)
    {
        ost << arr[i];
        if (i < n-1)
            ost << " ";
    }
    ost << '\n';
    return 0;
}