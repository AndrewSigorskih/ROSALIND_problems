#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <sstream>
#include <vector>
#include <map>
using namespace std;

void read_edges(ifstream& ist, map<int, int>& m)
{
    string input;
    while(getline(ist, input)) 
    {
        istringstream iss(input);
        int temp;
        while(iss >> temp)
        {
            m[temp]++;
        }
    }
}

int main()
{
    ifstream ist{"rosalind_deg.txt"};
    if (!ist) 
    {
	cout << "Cannot open input file!\n";
	exit(1);
    }
    ist.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    map<int, int> m;
    read_edges(ist, m);

    for ( const auto &myPair : m ) 
    {
        cout << myPair.second << " ";
    }
    cout << "\n";

    return 0;
}