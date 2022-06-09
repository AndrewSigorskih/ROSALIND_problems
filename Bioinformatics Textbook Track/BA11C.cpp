#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

using namespace std;

struct A
{
    static map<char,int> create_map()
        {
          map<char,int> m;
          m['X'] = 4;
          m['Z'] = 5;
          m['A'] = 71;
          m['R'] = 156;
          m['N'] = 114;
          m['D'] = 115;
          m['C'] = 103;
          m['E'] = 129;
          m['Q'] = 128;
          m['G'] = 57;
          m['H'] = 137;
          m['I'] = 113;
          m['L'] = 113;
          m['K'] = 128;
          m['M'] = 131;
          m['F'] = 147;
          m['P'] = 97;
          m['S'] = 87;
          m['T'] = 101;
          m['W'] = 186;
          m['Y'] = 163;
          m['V'] = 99;
          return m;
        }
    static const map<char,int> myMap;

};

const map<char,int> A:: myMap =  A::create_map();

void print_vector(vector<int>& vec)
{
    for (int i = 0; i < vec.size(); i++)
    {
        if (i == 0)
            cout << "" << vec[i];
        else
            cout << " " << vec[i];
    }
    cout << "\n";
}

int main() 
{
    ifstream ist{"rosalind_ba11c.txt"};
	if (!ist) 
	{
		cout << "Cannot open input file!\n";
		exit(1);
	}

    string input;
    getline(ist, input);
    
    int curmass = 0;
    vector<int> masses;

    for (int i = 0; i < input.size(); i++)
    {
        curmass += A::myMap.at(input[i]);
        masses.push_back(curmass);
    }
    //print_vector(masses);
    vector<int> result(curmass, 0);
    
    for (int i = 0; i < masses.size(); i++)
        result[masses[i]-1] = 1;

    print_vector(result);
    return 0;
}