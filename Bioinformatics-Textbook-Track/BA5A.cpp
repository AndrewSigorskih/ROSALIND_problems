#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using std::string;
using std::vector;
const char* INFILENAME = "rosalind_ba5a.txt";
const char* OUTFILENAME = "out.txt";

vector<int> parseCSVline(std::ifstream& ist)
{
    string buf;
    vector<int> result;
    getline(ist, buf);
    std::stringstream ss(buf);
    for (int i; ss >> i;) 
    {
        result.push_back(i);    
        if (ss.peek() == ',')
            ss.ignore();
    }
    return result;
}

int DPChange(int const M, vector<int> const& coins)
{
    vector<int> tab(M+1, INT32_MAX);
    tab[0] = 0;
    for (int m = 1; m <= M; ++m)
        for (int i = 0; i < coins.size(); ++i)
            if ((m >= coins[i]) && (tab[m-coins[i]] + 1 < tab[m]))
                    tab[m] = tab[m-coins[i]] + 1;
    return tab[M];
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

    int M;
    vector<int> coins;
    ist >> M;
    ist.ignore();
    coins = parseCSVline(ist);
    ost << DPChange(M, coins) << '\n';
    return 0;
}