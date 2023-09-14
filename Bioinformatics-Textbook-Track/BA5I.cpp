#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm> // max element, reverse

enum Direction
{
    DIAG,
    UP,
    LEFT
};

int const MATCH = 1;
int const GAP_PEN = -2;
using std::vector;
using std::string;
using array2d = vector<vector<int>>;
using ptrs2d = vector<vector<Direction>>;
const char* INFILENAME = "rosalind_ba5i.txt";
const char* OUTFILENAME = "out.txt";

class Aligner
{
    string s1, s2;
    array2d tab;
    ptrs2d pntrs;
public:
    Aligner(std::istream&);
    void align(std::ofstream&);
};

Aligner::Aligner(std::istream& ist)
{
    ist >> s1 >> s2;
    tab = array2d(s1.size()+1, vector<int>(s2.size()+1));
    pntrs = ptrs2d(s1.size()+1, vector<Direction>(s2.size()+1));
}

void Aligner::align(std::ofstream& ost)
{
    int m, n;
    m = s1.size();
    n = s2.size();
    for (int i = 0; i < m+1; ++i)
    {
        tab[i][0] = 0;
        pntrs[i][0] = LEFT;
    }
    for (int j = 0; j < n+1; ++j)
    {
        tab[0][j] = j * GAP_PEN;
        pntrs[0][j] = UP;
    }
    for (int i = 1; i < m+1; ++i)
    {
        for (int j = 1; j < n+1; ++j)
        {
            vector<int> scores = {
                tab[i-1][j-1] + (s1[i-1] == s2[j-1] ? MATCH : GAP_PEN),
                tab[i][j-1] + GAP_PEN,
                tab[i-1][j] + GAP_PEN
            };
            auto max_iter = std::max_element(scores.begin(), scores.end());
            tab[i][j] = *max_iter;
            pntrs[i][j] = static_cast<Direction>(max_iter - scores.begin());
        }
    }
    // find max in last row
    string s1_al, s2_al;
    int max_score = tab[m][0], j = 0, i = m;
    for (int k = 1; k < n+1; ++k)
    {
        if (tab[m][k] >= max_score)
        {
            max_score = tab[m][k];
            j = k;
        }
    }
    // backtrack
    while ((i>0) && (j > 0))
    {
        switch (pntrs[i][j])
        {
            case DIAG:
                s1_al += s1[--i];
                s2_al += s2[--j];
                break;
            case UP:
                s1_al += "-";
                s2_al += s2[--j];
                break;
            case LEFT:
                s1_al += s1[--i];
                s2_al += "-";
                break;
        }
    }
    
    std::reverse(s1_al.begin(), s1_al.end());
    std::reverse(s2_al.begin(), s2_al.end());
    ost << max_score << "\n";
    ost << s1_al << "\n" << s2_al << "\n";
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
    
    Aligner alner(ist);
    alner.align(ost);

    return 0;
}