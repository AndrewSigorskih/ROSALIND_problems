#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <sstream>
#include <algorithm>
using namespace std;

void read_array(ifstream& ist, int* arr, int n)
{
    string input;
    getline(ist, input);
    istringstream iss(input);
    int temp;
    for (int i = 0; i < n; i++)
    {
        iss >> temp;
        *(arr + i) = temp;
    }
}

void printHeap(int arr[], int n)
{ 
    for (int i = 0; i < n; ++i)
        cout << arr[i] << " ";
    cout << "\n";
}

int main()
{
    ifstream ist{"rosalind_ps.txt"};
	if (!ist) 
	{
		cout << "Cannot open input file!\n";
		exit(1);
	}

    int n, k;
    ist >> n;
    ist.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    int* arr = new int[n];
    read_array(ist, arr, n);

    ist >> k;
    //std::ps for now _)
    partial_sort(arr, arr+k, arr+n);
 
    printHeap(arr, k);

    delete [] arr;
    return 0;
}