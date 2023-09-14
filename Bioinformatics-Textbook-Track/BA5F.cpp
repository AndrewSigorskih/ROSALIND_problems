#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm> // max element, reverse
#include <map>

enum Direction
{
    DIAG,
    UP,
    LEFT,
    STOP
};

int const GAP_PEN = -5;
using std::vector;
using std::string;
using array2d = vector<vector<int>>;
using ptrs2d = vector<vector<Direction>>;
using Matrix = std::map<string, int>;
const char* INFILENAME = "rosalind_ba5f.txt";
const char* OUTFILENAME = "out.txt";
const char* MATRIXFILENAME = "PAM250";

void findMaxElement(array2d& tab, int& row, int& col, int& max)
{
    max = tab[0][0];
    for (int i = tab.size()-1; i >= 0; --i)
    {
        for (int j = tab[0].size()-1; j >= 0; --j)
        {
            if (tab[i][j] >= max)
            {
                max = tab[i][j];
                row = i;
                col = j;
            }
        }
    }
}

class Aligner
{
    string s1, s2;
    Matrix mat;
    array2d tab;
    ptrs2d pntrs;
    vector<string> parseLine(std::istream&);
    void readMatrix();
public:
    Aligner(std::istream&);
    void align(std::ofstream&);
};

Aligner::Aligner(std::istream& ist)
{
    readMatrix();
    ist >> s1 >> s2;
    tab = array2d(s1.size()+1, vector<int>(s2.size()+1));
    pntrs = ptrs2d(s1.size()+1, vector<Direction>(s2.size()+1));
}

vector<string> Aligner::parseLine(std::istream& ist)
{
    string buf;
    vector<string> result;
    std::getline(ist, buf);
    std::istringstream iss(buf);
    std::copy(std::istream_iterator<string>(iss), 
              std::istream_iterator<string>(),
              std::back_inserter(result));
	return result;
}

void Aligner::readMatrix()
{
    std::ifstream matrix_file{MATRIXFILENAME};
    if (!matrix_file) 
    {
        std::cerr << "Cannot open substitution matrix file!\n";
        exit(1);
    }
    vector<string> header = parseLine(matrix_file);
    for (int i = 0; i < header.size(); ++i)
    {
        vector<string> buf = parseLine(matrix_file);
        for (int j = 1; j <= header.size(); ++j)
        {
            mat[buf[0] + header[j-1]] = stoi(buf[j]);
        }
        buf.clear();
    }
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
        tab[0][j] = 0;
        pntrs[0][j] = UP;
    }
    for (int i = 1; i < m+1; ++i)
    {
        for (int j = 1; j < n+1; ++j)
        {
            vector<int> scores = {
                tab[i-1][j-1] + mat[s1.substr(i-1,1) + s2.substr(j-1, 1)],
                tab[i][j-1] + GAP_PEN,
                tab[i-1][j] + GAP_PEN,
                0 // end of local alignment
            };
            auto max_iter = std::max_element(scores.begin(), scores.end());
            tab[i][j] = *max_iter;
            pntrs[i][j] = static_cast<Direction>(max_iter - scores.begin());
        }
    }
    // jump to highest scoring local alignment
    int i, j, max_score;
    string s1_al, s2_al;
    findMaxElement(tab, i, j, max_score);
    // backtrack
    bool to_stop = false;
    while ((i>0) || (j > 0))
    {
        if (to_stop)
            break;
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
            case STOP:
                to_stop = true;
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