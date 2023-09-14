#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm> // max element
#include <map>

enum Direction
{
    LOWER = 0,
    MIDDLE,
    UPPER
};

using std::vector;
using std::string;
using Matrix = std::map<string, int>;

template <typename T>
using array2d = vector<vector<T>>;
template <typename T>
using array3d = vector<array2d<T>>;

const int GAPOPEN = -11;
const int GAPEXTEND = -1;
const char* INFILENAME = "rosalind_ba5j.txt";
const char* OUTFILENAME = "out.txt";
const char* MATRIXFILENAME = "BLOSUM62";

void findMax(vector<int>& scores, std::pair<int, Direction>& result)
{
    auto max_iter = std::max_element(scores.begin(), scores.end());
    result.first = *max_iter;
    result.second = static_cast<Direction>(max_iter - scores.begin());
}

class Aligner
{
    string s1, s2;
    Matrix mat;
    array3d<int> tab;
    array3d<Direction> pntrs;
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
    tab = array3d<int>(3, 
                    array2d<int>(s1.size()+1, vector<int>(s2.size()+1))
                        );
    pntrs = array3d<Direction>(3, 
                    array2d<Direction>(s1.size()+1, vector<Direction>(s2.size()+1))
                                );
                            
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
    // prepare tensors
    for (int i = 1; i < m+1; ++i)
    {
        tab[0][i][0] = GAPOPEN;
        tab[1][i][0] = GAPOPEN;
        tab[2][i][0] = INT16_MIN;
    }
    for (int j = 1; j < n+1; ++j)
    {
        tab[0][0][j] = INT16_MIN;
        tab[1][0][j] = GAPOPEN;
        tab[2][0][j] = GAPOPEN;
    }
    //  fill table and keep backtrack
    for (int i = 1; i < m+1; ++i)
    {
        for (int j = 1; j < n+1; ++j)
        {
            vector<int> scores;
            std::pair<int, Direction> tmp;
            // upper
            scores = { tab[2][i][j-1] + GAPEXTEND, tab[1][i][j-1] + GAPOPEN };
            findMax(scores, tmp);
            tab[2][i][j] = tmp.first;
            pntrs[2][i][j] = tmp.second;
            // lower
            scores = { tab[0][i-1][j] + GAPEXTEND, tab[1][i-1][j] + GAPOPEN };
            findMax(scores, tmp);
            tab[0][i][j] = tmp.first;
            pntrs[0][i][j] = tmp.second;
            // middle
            scores = { 
                    tab[0][i][j], 
                    tab[1][i-1][j-1] + mat[s1.substr(i-1,1) + s2.substr(j-1, 1)], 
                    tab[2][i][j] 
                    };
            findMax(scores, tmp);
            tab[1][i][j] = tmp.first;
            pntrs[1][i][j] = tmp.second;
        }
    }
    // backtrack
    int i = m, j = n;
    string s1_al = s1, s2_al = s2;
    std::pair<int, Direction> tmp;
    vector<int> scores = {
        tab[0][i][j],
        tab[1][i][j],
        tab[2][i][j]
    };
    findMax(scores, tmp);
    Direction back_pntr = tmp.second;
    while ((i != 0) && (j != 0))
    {
        switch (back_pntr)
        {
            case LOWER:
                // lower matrix
                if (pntrs[0][i][j] == MIDDLE)
                    back_pntr = MIDDLE;
                --i;
                s2_al.insert(j, "-");
                break;
            case MIDDLE:
                // middle matrix
                if (pntrs[1][i][j] == LOWER)
                { 
                    back_pntr = LOWER;
                } else if (pntrs[1][i][j] == UPPER) {
                    back_pntr = UPPER;
                } else {
                    --i;
                    --j;
                }
                break;
            case UPPER:
                // upper matrix
                if (pntrs[2][i][j] == MIDDLE)
                    back_pntr = MIDDLE;
                --j;
                s1_al.insert(i, "-");
                break;
        }
    }
    for (; i > 0; --i)
        s1_al.insert(0, "-");
    for (; j > 0; --j)
        s2_al.insert(0, "-");
    // save result
    ost << tmp.first << "\n";
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