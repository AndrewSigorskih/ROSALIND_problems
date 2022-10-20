#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

using std::vector;
using std::string;
const char* INFILENAME = "rosalind_ba6b.txt";
const char* OUTFILENAME = "out.txt";

void parsePermutation(string const& input, vector<int>& result)
{
    std::stringstream ss(input);
    int tmp;
    while (ss >> tmp)
        result.push_back(tmp);
}

int countBreakpoints(vector<int> const& permutation)
{
    int result = 0;
    for (int i = 1; i < permutation.size(); ++i)
        if (permutation[i-1] + 1  != permutation[i])
            ++result;
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
    vector<int> permutation = {0};
    getline(ist, input);
    parsePermutation(input.substr(1, input.size()-2), permutation);
    permutation.push_back(permutation.size());
    ost << countBreakpoints(permutation) << '\n';
    return 0;
}