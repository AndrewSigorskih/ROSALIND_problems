#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

using namespace std;

struct A
{
    static map<int,char> create_map()
        {
            map<int,char> m;
            m[4] = 'X';
            m[5] = 'Z';
            m[57] = 'G';
            m[71] = 'A';
            m[87] = 'S';
            m[97] = 'P';
            m[99] = 'V';
            m[101] = 'T';
            m[103] = 'C';
            m[113] = 'L'; // sorry ILE
            m[114] = 'N';
            m[115] = 'D';
            m[128] = 'Q'; // sorry LYS
            m[129] = 'E';
            m[131] = 'M';
            m[137] = 'H';
            m[147] = 'F';
            m[156] = 'R';
            m[163] = 'Y';
            m[186] = 'W';
          return m;
        }
    static const map<int,char> myMap;
};

const map<int,char> A:: myMap =  A::create_map();

int main()
{
    ifstream ist{"rosalind_ba11d.txt"};
	if (!ist) 
	{
		cout << "Cannot open input file!\n";
		exit(1);
	}

    string input;
    getline(ist, input);

    int i = 0;
    for (char c: input) 
    {
        switch(c)
        {
            case '0':
                i++;
                break;
            case '1':
                cout << A::myMap.at(i+1);
                i = 0;
                break;
            default:
                break;
        }
    }
    cout << "\n";

    return 0;
}