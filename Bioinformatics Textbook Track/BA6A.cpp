#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

using std::vector;
using std::string;
const char* INFILENAME = "rosalind_ba6a.txt";
const char* OUTFILENAME = "out.txt";

void parsePermutation(string const& input, vector<int>& result)
{
    std::stringstream ss(input);
    int tmp;
    while (ss >> tmp)
        result.push_back(tmp);
}

void printPermutation(vector<int> const& permutation, std::ofstream& ost)
{
    ost << "(";
    for (int i = 0; i < permutation.size()-1; ++i)
        ost << permutation[i] << " ";
    ost << permutation[permutation.size()-1] << ")\n";
}

void reversePermutation(vector<int>& permutation, int ind)
{
    auto start = permutation.begin() + (ind - 1);
    auto end = permutation.begin();
    for (int i = (ind-1); i < permutation.size(); ++i)
    {
        if ((permutation[i] == ind) || (permutation[i] == -ind))
        {
            end += i;
            break;
        }
    }
    for (auto it = start; it <= end; ++it)
        (*it) *= (-1);
    std::reverse(start, end+1);
}

void greedySorting(vector<int>& permutation, std::ofstream& ost)
{
    int ind = 1;
    while (ind <= permutation.size())
    {
        if (permutation[ind-1] == ind)
        {
            ++ind;
        } else {
            reversePermutation(permutation, ind);
            printPermutation(permutation, ost);
        }     
    }
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
    ost << std::showpos; // sign of int
    string input;
    vector<int> permutation;
    getline(ist, input);
    parsePermutation(input.substr(1, input.size()-2), permutation);
    greedySorting(permutation, ost);
    return 0;
}