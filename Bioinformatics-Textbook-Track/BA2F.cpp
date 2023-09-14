#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <random>

using std::vector;
using std::string;
using Profile = vector<vector<double>>;
const int N_ITER = 1000;
const int PSEUDOCOUNT = 1;
const char* BASES = "ACGT";
const char* INFILENAME = "rosalind_ba2f.txt";
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
    int k;
    vector<string> dna;
    Profile profile;
    vector<string> best_motifs;
    vector<string> estimateMotifs();
    int score(vector<string> const&);
    string profMostProbable(string const&);
    void createProfile(vector<string> const&);
    
public:
    motifSearcher(std::ifstream&);
    void randomizedMotifSearch(std::ofstream&);
};

motifSearcher::motifSearcher(std::ifstream& ist)
{
    int t; // useless variable == dna.size()
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

vector<string> motifSearcher::estimateMotifs()
{
    vector<string> result;
    for (auto text: dna)
    {
        result.push_back(profMostProbable(text));
    }
    return result;
}

void motifSearcher::randomizedMotifSearch(std::ofstream& ost)
{
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, dna[0].size() - k - 1);

    int best_score = INT32_MAX;

    for (int iter = 1; iter <= N_ITER; ++iter)
    {  // randomly select k-mers Motifs
        vector<string> motifs;
        vector<string> leadingMotifs;
        motifs.reserve(dna.size());
        for (int i = 0; i < dna.size(); ++i)
        {
            int randind = distr(gen);
            motifs.push_back(dna[i].substr(randind, k));
        }
        leadingMotifs = motifs;
        while (true)
        {
            createProfile(motifs);
            motifs = estimateMotifs();
            if (score(motifs) < score(leadingMotifs))
                leadingMotifs = motifs;
            else    
                break; // search untill minimizing stops
        }
        
        int cur_score = score(leadingMotifs);
        if (cur_score < best_score)
        { // update global min
            best_score = cur_score;
            best_motifs = leadingMotifs;
        }
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
    processor.randomizedMotifSearch(ost);
    return 0;
}