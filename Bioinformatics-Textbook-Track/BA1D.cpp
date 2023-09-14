#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using std::string;
using std::vector;
const char* INFILENAME = "rosalind_ba1d.txt";
const char* OUTFILENAME = "out.txt";

vector<int> computeSuffPref(const string& pat)
{
    vector<int> result(pat.size(), 0);
    int len = 0;
    result[0] = 0;
    int i = 1;
    while (i < pat.size())
    {
        if (pat[i] == pat[len])
        {   // len++; result[i] = len; i++;
            result[i++] = ++len;
        } else { // (pat[i] != pat[len])
            if (len != 0)
            {
                len = result[len - 1];
            } else { // len == 0)
                result[i] = 0;
                ++i;
            }
        }
    }
    return result;
}

vector<int> KMPsearch(const string& pat, const string& text)
{
    vector<int> result;
    int M = pat.size();
    int N = text.size();

    vector<int> table = computeSuffPref(pat);

    int i = 0; // index for text
    int j = 0; // index for pat

    while ((N - i) >= (M - j)) {
        if (pat[j] == text[i]) 
        {
            j++;
            i++;
        }
 
        if (j == M) 
        {
            //printf("Found pattern at index %d ", i - j);
            result.push_back(i - j);
            j = table[j - 1];
        // mismatch after j matches
        } else if ((i < N) && (pat[j] != text[i])) { 
            // Do not match table[0..table[j-1]] characters,
            // they will match anyway
            if (j != 0)
                j = table[j - 1];
            else
                i = i + 1;
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
    string pat, genome;
    std::getline(ist, pat);
    std::getline(ist, genome);

    vector<int> result = KMPsearch(pat, genome);
    for(auto c: result)
    {
        ost << c << " ";
    }
    ost << '\n';
    return 0;
}