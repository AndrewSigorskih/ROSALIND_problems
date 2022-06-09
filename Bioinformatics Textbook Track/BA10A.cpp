#include <iostream>
#include <fstream>
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
map<string, map<string, double>> readMatrix(vector<string>& alphabet,
									vector<string>& header, std::istream& ist)
{
	map<string, map<string, double>> m;
	for (int i = 0; i < alphabet.size(); i++)
	{
		vector<string> buf = parseLine(ist);

		for (int j = 1; j <= alphabet.size(); j++)
		{
			m[buf[0]][header[j-1]] = double{ stod(buf[j]) };
		}
		buf.clear();
	}
	return m;
}

double calcPathProb(string& path, map<string, map<string, double>> m, 
					int size)
{
	double prob = 1.0 / size;
	for (int i = 1; i < path.size(); i++)
	{
		prob *= m[path.substr(i-1,1)][path.substr(i,1)];
	}
	return prob;
}

int main()
{
	ifstream ist{"rosalind_ba10a.txt"};
	if (!ist) 
	{
		cout << "Cannot open input file!\n";
		exit(1);
	}

	string path;
	getline(ist, path);

	ist.ignore(2048, '\n');

	vector<string> alphabet = parseLine(ist);

	ist.ignore(2048, '\n');

	vector<string> header = parseLine(ist);
	map<string, map<string, double>> m = readMatrix(alphabet, header, ist);

	cout.precision(12);
	cout << calcPathProb(path, m, alphabet.size()) <<"\n";
	return 0;
}
