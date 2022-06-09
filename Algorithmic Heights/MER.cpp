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

vector<int> merge(vector<int>& left, vector<int>& right)
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

int main()
{
    ifstream ist{"rosalind_mer.txt"};
	if (!ist) 
	{
		cout << "Cannot open input file!\n";
		exit(1);
	}

    int n,m;
    vector<int> left, right;

    ist >> n;
    ist.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    read_int_array(ist, left);

    ist >> m;
    ist.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    read_int_array(ist, right);

    vector<int> result = merge(left, right);
    for (int i = 0; i < result.size(); i++)
    {
        cout << result[i];
        if (i < result.size()-1)
            cout << " ";
    }
    cout << endl;

    return 0;
}