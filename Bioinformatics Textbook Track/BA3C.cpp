#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using std::vector;
using std::string;
const char* INFILENAME = "rosalind_ba3c.txt";
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

    vector<string> reads;
    reads.reserve(64); // avoiding too many reallocations
    string tmp;
    while(getline(ist, tmp))
        reads.push_back(tmp);
    // create graph overlap
    int k = reads[0].size();
    for (int i = 0; i < reads.size(); ++i)
    {
        for (int j = 0; j < reads.size(); ++j)
        {
            if (i == j)
                continue;
            if (reads[i].substr(1, k-1) == reads[j].substr(0, k-1))
                ost << reads[i] << " -> " << reads[j] << '\n';
        }
    }
    return 0;
}