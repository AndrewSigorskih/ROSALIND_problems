#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using std::vector;
using std::string;
using matrix = vector<vector<int>>;
const char* INFILENAME = "rosalind_ba5c.txt";
const char* OUTFILENAME = "out.txt";

inline int max(int a, int b)
{
    return ((a > b) ? (a) : (b));
}

string lcs(string const& s1, string const& s2)
{   
    int m = s1.size();
    int n = s2.size();
    matrix tab(m+1, vector<int>(n+1));
    // fill dp table
    for (int i = 0; i <= m; ++i)
    {
        for (int j = 0; j <= n; ++j)
        {
            if ((i == 0) || (j == 0))
                tab[i][j] = 0;
            else if (s1[i-1] == s2[j-1])
                tab[i][j] = tab[i-1][j-1] + 1;
            else
                tab[i][j] = max(tab[i-1][j], tab[i][j-1]);
        }
    }
    // backtrack lcs
    int index = tab[m][n];
    //char result[index];
    char* tmp = new char[index + 1];
    tmp[index] = '\0';
    int i = m, j = n;
    while (i > 0 && j > 0)
    {   // same symbols -> part of LCS
        if (s1[i-1] == s2[j-1])
        {
            tmp[index-1] = s1[i-1];
            i--;
            j--;
            index--; 
        } // if not same -> go larger way
        else if (tab[i-1][j] > tab[i][j-1])
            i--;
        else
            j--;
    }
    string result(tmp);
    delete[] tmp;
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

    string s1, s2;
    ist >> s1 >> s2;

    ost << lcs(s1, s2) << '\n';

    return 0;
}
