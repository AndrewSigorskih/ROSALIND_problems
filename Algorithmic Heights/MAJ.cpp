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

//Boyer–Moore majority vote algorithm
//linear time and space complexity
//step 1: find the candidate for Majority
int findCandidate(vector<int>& a)
{
    int maj_index = 0, count = 1;
    for (int i = 1; i < a.size(); i++) 
    {
        if (a[maj_index] == a[i])
            count++;
        else
            count--;
        if (count == 0) 
        {
            maj_index = i;
            count = 1;
        }
    }
    return a[maj_index];
}

//step 2: check if the candidate occurs more than n/2 times
bool isMajority(vector<int>& a, int& cand)
{
    int count = 0;
    for (int i = 0; i < a.size(); i++)
 
        if (a[i] == cand)
            count++;
 
    if (count > a.size() / 2)
        return 1;
    else
        return 0;
}

void printMajority(vector<int>& a)
{
    int cand = findCandidate(a);
 
    if (isMajority(a, cand))
        cout << cand;
    else
        cout << -1;
}

int main()
{
    ifstream ist{"rosalind_maj.txt"};
	if (!ist) 
	{
		cout << "Cannot open input file!\n";
		exit(1);
	}

    int k, n;
    ist >> k >> n;
    ist.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    for (int i = 0; i < k; i++)
    {
        vector<int> arr;
        read_int_array(ist, arr);
        printMajority(arr);
        if (i < (k-1))
            cout << " ";
    }   
    cout << endl;

    return 0;
}