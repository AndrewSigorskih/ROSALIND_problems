#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
using namespace std;

void read_int_array(ifstream& ist, vector<int>& arr)
{
    string input;
    getline(ist, input);
    istringstream iss(input);
    int temp;
    while(iss >> temp)
    {
    arr.push_back(temp);
    }
}

int binary_search(vector<int>& arr, int start, int end, int x)
{   
    if (end >= start)
    {
        int mid  = start + (end - start) / 2;
        if (arr[mid] == x)
        {
            return mid + 1;
        }
        if (arr[mid] > x)
        {
            return binary_search(arr, start, end-1, x);
        }
        return binary_search(arr, start+1, end, x);
    }
    return -1;
}

int main()
{
    ifstream ist{"rosalind_bins.txt"};
	if (!ist) 
	{
		cout << "Cannot open input file!\n";
		exit(1);
	}
    int n, m;
    vector<int> arr, list;

    ist >> n;
    ist >> m;
    ist.get();
    read_int_array(ist, arr);
    read_int_array(ist, list);
    
    for(int i = 0; i < m; i++)
    {   
        cout << binary_search(arr, 0, arr.size()-1, list[i]) << " ";
    }
    cout << "\n";
    return 0;
}