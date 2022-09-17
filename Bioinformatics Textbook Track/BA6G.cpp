#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <stdlib.h> // atoi

using std::string;
using std::vector;
const char* INFILENAME = "rosalind_ba6g.txt";
const char* OUTFILENAME = "out.txt";

vector<int> parseLine(const string& buf)
{
    vector<string> tmp;
    std::istringstream iss(buf);
    std::copy(std::istream_iterator<string>(iss), 
              std::istream_iterator<string>(),
              std::back_inserter(tmp));
    vector<int> result(tmp.size());
    for (int i = 0; i < tmp.size(); ++i)
        result[i] = atoi(tmp[i].c_str());
	return result;
}

vector<int> CycleToChromosome(const vector<int>& nodes)
{
    vector<int> result(nodes.size()/2);
    for (int i = 0; i < result.size(); ++i)
    {
        if (nodes[2*i] < nodes[2*i + 1])
        {
            result[i] = nodes[2*i + 1] / 2;
        } else {
            result[i] = -(nodes[2*i]) / 2;
        }
    }
    return result;
}

int main()
{   
    std::ifstream ist{INFILENAME};
    if (!ist) 
    {
        std::cerr << "Cannot open input file!\n";
        exit(1);
    }
    std::ofstream ost{OUTFILENAME};
    if (!ost) 
    {
        std::cerr << "Cannot open output file!\n";
        exit(1);
    }

    string input;
    getline(ist, input);
    input = input.substr(1, input.size()-2); // trim brackets
    vector<int> nodes = parseLine(input);
    vector<int> result = CycleToChromosome(nodes);
    ost << "(";
    for (int i = 0; i < result.size(); ++i)
    {
        if (result[i] > 0)
            ost << "+" << result[i];
        else
            ost << result[i];
        if (i < result.size()-1)
            ost << " ";
    }
    ost << ")\n";
    return 0;
}