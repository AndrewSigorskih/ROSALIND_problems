#include <iostream>
#include <fstream>
#include <limits>
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

int insertionSort(vector<int>& arr)
{
    int inversions = 0;
    for (int i = 1; i < arr.size(); i++)
    {
        int j = i;
        while((j > 0) && (arr[j-1] > arr[j]))
        {
            int temp = arr[j-1];
            arr[j-1] = arr[j];
            arr[j] = temp;
            inversions++;
            j--;
        }
    }
    return inversions;
}


int main()
{
    ifstream ist{"rosalind_ins.txt"};
	if (!ist) 
	{
		cout << "Cannot open input file!\n";
		exit(1);
	}

    int k;
    ist >> k;
    ist.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    vector<int> arr;
    read_int_array(ist, arr);

    cout << insertionSort(arr) << endl;

    return 0;
}