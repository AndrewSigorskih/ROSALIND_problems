#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <stdlib.h> // atoi

using std::string;
using std::vector;
const char* INFILENAME = "rosalind_ba6f.txt";
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

vector<int> ChromosomeToCycle(const vector<int>& chromosome)
{
    vector<int> result(chromosome.size()*2);
    for (int i = 0; i < chromosome.size(); ++i)
    {
        int val = chromosome[i];
        if (val > 0)
        {
            result[2*i] = 2*val - 1;
            result[2*i+1] = 2*val;
        } else {
            result[2*i] = -2*val;
            result[2*i+1] = -2*val - 1;
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
    vector<int> chromosome = parseLine(input);
    for (auto i: ChromosomeToCycle(chromosome))
        ost << i << " ";
    ost << '\n';
    return 0;
}