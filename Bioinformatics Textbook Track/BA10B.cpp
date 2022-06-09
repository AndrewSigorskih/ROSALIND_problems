#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <sstream>
#include <iterator>
#include <map>

using namespace std;

//get line and split by any number of whitespaces
vector<string> parseLine(std::istream& str)
{
	vector<std::string>   result;
    string                buf;
    getline(str, buf);

	istringstream iss(buf);
    copy(istream_iterator<string>(iss), istream_iterator<string>(),
              back_inserter(result));

	return result;
} 
/*
void printVector(vector<string>& vec)
{
	for (string s: vec)
    	cout << s << " ";
	cout << "\n";
}
*/
map<string, map<string, double>> readMatrix(vector<string>& states,
									vector<string>& header, std::istream& ist)
{
	map<string, map<string, double>> m;
	for (int i = 0; i < states.size(); i++)
	{
		vector<string> buf = parseLine(ist);

		for (int j = 1; j <= header.size(); j++)
		{
			m[buf[0]][header[j-1]] = double{ stod(buf[j]) };
		}
		buf.clear();
	}
	return m;
}

double calcPathProb(string& str, string& path, map<string, map<string, double>>& m)
{
    double prob = 1.0;
    for (int i = 0; i < str.size(); i++)
    {
        prob *= m[path.substr(i,1)][str.substr(i,1)];
        
    }
	return prob;
}

int main()
{
	ifstream ist{"rosalind_ba10b.txt"};
	if (!ist) 
	{
		cout << "Cannot open input file!\n";
		exit(1);
	}

    string str;
    getline(ist, str);

    ist.ignore(numeric_limits<streamsize>::max(), '\n');

    vector<string> alphabet = parseLine(ist);

    ist.ignore(numeric_limits<streamsize>::max(), '\n');

	string path;
	getline(ist, path);

	ist.ignore(numeric_limits<streamsize>::max(), '\n');

	vector<string> states = parseLine(ist);

	ist.ignore(numeric_limits<streamsize>::max(), '\n');

	vector<string> header = parseLine(ist);
	map<string, map<string, double>> emissions = readMatrix(states, header, ist);

	cout.precision(12);
	cout << calcPathProb(str, path, emissions) <<"\n";
	return 0;
}