#include <iostream>
#include <fstream>
#include <vector>

using std::vector;
using matrix = vector<vector<int>>;
const char* INFILENAME = "rosalind_ba5b.txt";
const char* OUTFILENAME = "out.txt";

inline int max(int a, int b)
{
    return ((a > b) ? (a) : (b));
}

int ManhattanTourist(matrix const& down, matrix const& right)
{
    int n = down.size();
    int m = right[0].size();
    matrix tab(n+1, vector<int>(m+1));

    for (int i = 1; i <= n; ++i)
        tab[i][0] = tab[i-1][0] + down[i-1][0];
    for (int j = 1; j <= m; ++j)
        tab[0][j] = tab[0][j-1] + right[0][j-1];
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= m; ++j)
            tab[i][j] = max(
                tab[i-1][j] + down[i-1][j],
                tab[i][j-1] + right[i][j-1]
            );
    return tab[n][m];
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

    int n, m;
    ist >> n >> m;
    matrix down(n, vector<int>(m+1));
    matrix right(n+1, vector<int>(m));

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < (m+1); ++j)
            ist >> down[i][j];
    ist.ignore(256, '-');
    for (int i = 0; i < (n+1); ++i)
        for (int j = 0; j < m; ++j)
            ist >> right[i][j];
    ost << ManhattanTourist(down, right) << '\n';
    return 0;
}