#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm> // min

using std::vector;
using std::string;
using matrix = vector<vector<int>>;
const char* INFILENAME = "rosalind_ba5g.txt";
const char* OUTFILENAME = "out.txt";

int editDistance(string const& s1, string const& s2)
{
    int m = s1.size();
    int n = s2.size();
    matrix tab(m+1, vector<int>(n+1));
    // fill table
    for (int i = 0; i <= m; ++i)
        tab[i][0] = i;
    for (int j = 0; j <= n; ++j)
        tab[0][j] = j;
    for (int i = 1; i <= m; ++i)
    {
        for (int j = 1; j <= n; ++j)
        {
            int cost = 0;
            if (s1[i-1] != s2[j-1])
                cost = 1;
            tab[i][j] = std::min({
                tab[i-1][j] + 1,
                tab[i][j-1] + 1,
                tab[i-1][j-1] + cost
            });
        }
    }
    return tab[m][n];
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

    string s1, s2;
    ist >> s1 >> s2;

    ost << editDistance(s1, s2) << '\n';

    return 0;
}