#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using std::string;
using std::vector;
const char* INFILENAME = "rosalind_ba1f.txt";
const char* OUTFILENAME = "out.txt";

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

    string genome;
    getline(ist, genome);
    vector<int> skews(genome.size());
    int C_cnt = 0, G_cnt = 0, min = INT32_MAX;

    for (int i = 0; i < genome.size(); ++i)
    {
        if (genome[i] == 'C')
            ++C_cnt;
        else if (genome[i] == 'G')
            ++G_cnt;
        skews[i] = (G_cnt - C_cnt);
        if (skews[i] < min)
            min = skews[i];
    }
    for (int i = 0; i < genome.size(); ++i)
        if (skews[i] == min)
            ost << (i+1) << " ";
    ost << '\n';
    return 0;
}