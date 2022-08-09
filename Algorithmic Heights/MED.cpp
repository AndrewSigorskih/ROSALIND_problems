#include <iostream>
#include <fstream>
#include <algorithm>
using std::swap;
using std::sort;
const char* INFILENAME = "rosalind_med.txt";
const char* OUTFILENAME = "out.txt";

// reference:
// https://rcoh.me/posts/linear-time-median-finding/

int partition(int arr[], int l, int r, int k); 
  
// called on arrays of fixed size (5) =>
// linear time 
int findMedian(int arr[], int n) 
{ 
    sort(arr, arr+n);
    return arr[n/2];
} 
  
// Returns k'th smallest element in arr[l..r] in worst case 
// linear time.
int kthSmallest(int arr[], int l, int r, int k) 
{ 
    if (k > 0 && k <= r - l + 1) 
    { 
        int n = r-l+1;
  
        // Divide arr[] in groups of size 5, calculate median 
        // of every group and store it in median[] array. 
        int i, median[(n+4)/5]; // There will be floor((n+4)/5) groups; 
        for (i=0; i<n/5; i++) 
            median[i] = findMedian(arr+l+i*5, 5); 
        if (i*5 < n) //For last group with less than 5 elements 
        { 
            median[i] = findMedian(arr+l+i*5, n%5);  
            i++; 
        }     
  
        int medOfMed = (i == 1)? median[i-1]: 
                                 kthSmallest(median, 0, i-1, i/2); 
  
        int pos = partition(arr, l, r, medOfMed); 
  
        if (pos-l == k-1) 
            return arr[pos]; 
        if (pos-l > k-1)
            return kthSmallest(arr, l, pos-1, k); 
  
        return kthSmallest(arr, pos+1, r, k-pos+l-1); 
    } 
    return -1; 
} 
  
// It searches for x in arr[l..r], and partitions the array  
// around x. 
int partition(int arr[], int l, int r, int x) 
{ 
    // Search for x in arr[l..r] and move it to end 
    int i; 
    for (i = l; i < r; i++) 
        if (arr[i] == x) 
           break; 
    swap(arr[i], arr[r]); 
  
    // Standard partition algorithm 
    i = l; 
    for (int j = l; j <= r - 1; j++) 
    { 
        if (arr[j] <= x) 
        { 
            swap(arr[i], arr[j]); 
            i++; 
        } 
    } 
    swap(arr[i], arr[r]); 
    return i; 
} 

int main()
{
    std::ifstream ist{INFILENAME};
    if (!ist) 
    {
        std::cout << "Cannot open input file!\n";
        exit(1);
    }

    int n, k;
    ist >> n;
    int* arr = new int[n];

    for (int i= 0; i < n; ++i)
    {
        ist >> arr[i];
    }
    ist >> k;

    std::ofstream ost{OUTFILENAME};
    if (!ost) 
    {
        std::cout << "Cannot open output file!\n";
        exit(1);
    }

    ost << kthSmallest(arr, 0, n-1, k) << '\n';
    delete[] arr;
    return 0;
}