#include <iostream>
#include <fstream>
#include <unordered_map>
#include <set>
using std::set;
using umap = std::unordered_map<int, u_long>;
const char* INFILENAME = "rosalind_ba4d.txt";
const char* OUTFILENAME = "out.txt";
set<int> MASSES;
umap MEMO;

void loadMasses()
{
    // using set here due to the fact
    // equal-mass amino acids should be counted once
    std::ifstream source{"masses.lst"};
    if (!source) 
    {
        std::cerr << "Cannot open file with AA masses!\n";
        exit(1);
    }
    char aa;
    double mass;
    while (source >> aa >> mass)
    {
        if ((aa == 'X') || (aa == 'Z'))
            continue;
        MASSES.insert(static_cast<int>(mass));
    }
}

u_long countPeptides(int mass)
{
    if (mass < 0)
        return 0;
    auto it = MEMO.find(mass);
    if (it != MEMO.end())
        return it->second;
    u_long result = 0;
    for (auto x: MASSES)
        result += countPeptides(mass - x);
    MEMO[mass] = result;
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

    loadMasses();
    MEMO[0] = 1;

    int mass;
    ist >> mass;
    ost << countPeptides(mass) << '\n';
    return 0;
}