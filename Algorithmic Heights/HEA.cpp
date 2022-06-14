#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <sstream>
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

void heapify(int* arr, int n, int i)
{
    int largest = i; 
    int l = 2 * i + 1;
    int r = 2 * i + 2;
    
    if ((l < n) && (arr[l] > arr[largest]))
        largest = l;
 
    if ((r < n) && (arr[r] > arr[largest]))
        largest = r;
 
    if (largest != i) 
    {
        swap(arr[i], arr[largest]);
        heapify(arr, n, largest);
    }
}

void buildHeap(int* arr, int n)
{
    // last non-leaf node
    int startIdx = ((n - 1) / 2);

    for (int i = startIdx; i >= 0; i--) 
    {
        heapify(arr, n, i);
    }
}

int main()
{
    ifstream ist{"rosalind_hea.txt"};
	if (!ist) 
	{
		cout << "Cannot open input file!\n";
		exit(1);
	}

    int n;
    ist >> n;
    ist.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    int* arr = new int[n];
    read_array(ist, arr, n);

    buildHeap(arr, n);
 
    printHeap(arr, n);

    delete [] arr;
    return 0;
}
