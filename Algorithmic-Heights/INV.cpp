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

vector<int> merge(vector<int>& left, vector<int>& right, int* invCount)
{
    vector<int> result;
    int i = 0, j = 0;
    while ((i < left.size()) && (j < right.size()))
    {
        if (left[i] <= right[j])
        {
            result.push_back(left[i]);
            i++;
        } else {
            result.push_back(right[j]);
            j++;
            *invCount += (left.size() - i);
        }
    }

    while (i < left.size())
    {
        result.push_back(left[i]);
        i++;
    }
    while (j < right.size())
    {
        result.push_back(right[j]);
        j++;
    }
    return result;
}

vector<int> mergeSort(vector<int>& arr, int* invCount)
{
    if (arr.size() <= 1)
        return arr;
    vector<int> left, right;
    for(int i = 0; i < arr.size(); i++)
    {
        if (i < arr.size() / 2)
            left.push_back(arr[i]);
        else
            right.push_back(arr[i]);
    }
    left = mergeSort(left, invCount);
    right = mergeSort(right, invCount);
    return merge(left, right, invCount);
}

int main()
{
    ifstream ist{"rosalind_inv.txt"};
	if (!ist) 
	{
		cout << "Cannot open input file!\n";
		exit(1);
	}

    int n;
    vector<int> arr;

    ist >> n;
    ist.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    read_int_array(ist, arr);

    int invCount = 0;
    vector<int> result = mergeSort(arr, &invCount);
    /*for (int i = 0; i < result.size(); i++)
    {
        cout << result[i];
        if (i < result.size()-1)
            cout << " ";
    }
    cout << endl;*/
    cout << invCount << endl;
    return 0;
}