#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

using std::vector;
using std::string;
using Profile = vector<vector<double>>;
const int PSEUDOCOUNT = 0;
const char* BASES = "ACGT";
const char* INFILENAME = "rosalind_ba2d.txt";
const char* OUTFILENAME = "out.txt";

inline int getIndex(char const x)
{
    switch (x)
    {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return -1;
    }
}

class motifSearcher
{
    int k, t;
    vector<string> dna;
    Profile profile;
    vector<string> best_motifs;
    int score(vector<string> const&);
    string profMostProbable(string const&);
    void createProfile(vector<string> const&);
public:
    motifSearcher(std::ifstream&);
    void greedyMotifSearch(std::ofstream&);
};

motifSearcher::motifSearcher(std::ifstream& ist)
{
    ist >> k >> t;
    ist.ignore(1, '\n');
    string tmp;
    while (getline(ist, tmp))
        dna.push_back(tmp);
    profile = Profile(4);
}

int motifSearcher::score(vector<string> const& motifs)
{
    int score = 0;
    for (int i = 0; i < motifs[0].size(); ++i)
    {
        vector<int> basecounts(4);
        for (int j = 0; j < motifs.size(); ++j)
            ++basecounts[getIndex(motifs[j][i])];
        char c = BASES[
            std::max_element(basecounts.begin(), basecounts.end()) - basecounts.begin()
        ];
        for (int j = 0; j < motifs.size(); ++j)
            if (motifs[j][i] != c)
                ++score;
    }
    return score;
}

void motifSearcher::createProfile(vector<string> const& motifs)
{
    for (int i = 0; i < profile.size(); ++i)
    {
        profile[i].clear();
        profile[i].reserve(motifs.size());
    }
        
    for (int c = 0; c < profile.size(); ++c)
    {
        for (int j = 0; j < motifs[0].size(); ++j)
        {
            double SUM = static_cast<double>(PSEUDOCOUNT);
            for (int i = 0; i < motifs.size(); ++i)
            {
                if (motifs[i][j] == BASES[c])
                    ++SUM;
            }
            profile[c].push_back(SUM / motifs.size());
        }
    }
}

string motifSearcher::profMostProbable(string const& text)
{
    string result;
    double max_prob = -1.0;
    for (int i = 0; i <= (text.size()-k); ++i)
    {
        double cur_prob = 1.0;
        string cur_kmer = text.substr(i, k);
        for (int j = 0; j < k; ++j)
        {
            cur_prob *= profile[getIndex(cur_kmer[j])][j];
        }
        if (cur_prob > max_prob)
        {
            max_prob = cur_prob;
            result = cur_kmer;
        }
    }
    return result;
}

void motifSearcher::greedyMotifSearch(std::ofstream& ost)
{
    for (int i = 0; i < dna.size(); ++i)
        best_motifs.push_back(dna[i].substr(0, k));
    for (int i = 0; i <= (dna[0].size()-k); ++i)
    {
        vector<string> motifs = {dna[0].substr(i, k)};
        for (int j = 1; j < dna.size(); ++j)
        {
            createProfile(motifs);
            motifs.push_back(profMostProbable(dna[j]));
        }
        if (score(motifs) < score(best_motifs))
            best_motifs = motifs;
    }
    for (auto item: best_motifs)
        ost << item << '\n';
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

    motifSearcher processor(ist);
    processor.greedyMotifSearch(ost);
    return 0;
}