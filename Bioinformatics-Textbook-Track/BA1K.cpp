#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath> // pow
#include <map>

using std::map;
using std::string;
using std::vector;
const char* INFILENAME = "rosalind_ba1k.txt";
const char* OUTFILENAME = "out.txt";

class kmerProcessor
{
    int _k;
    //int ind; // store pos in freq array 
    string const alphabet = "ACGT";
    map<string, int> kmers; // kmers from text
    vector<int> freq_arr; // frequency array
    void processAllKmersRec(string const&, int const);
public:
    kmerProcessor(string const& text, int k)
    {
        _k = k;
        //ind = 0;
        //freq_arr = vector<int>(pow(4, k), 0);
        freq_arr.reserve(pow(4, k));
        for (int i = 0; i <= (text.size()-k); ++i)
            kmers[text.substr(i, k)]++;
    }
    void printKmers(std::ofstream&);
    void processAllKmers(std::ofstream&);
};

void kmerProcessor::printKmers(std::ofstream& ost)
{
    for (auto const& item: kmers)
        ost << item.first << " ";
    ost << "\n";
}

void kmerProcessor::processAllKmersRec(string const& prefix, int const k)
{
    if (k == 0) // resursion stop: process current kmer
    {   
        int count = 0; // how many kmers from text match current one
        for (auto const& item: kmers)
        {
            if (prefix == item.first) // if kmer in text
            {
                count += item.second; // store its freq
                break; // only looking for exact matches here 
            }
        }
        freq_arr.push_back(count);
        return;
    }
    for (int i = 0; i < alphabet.size(); ++i)
    {   // far all possible combinations in ith pos
        string newPrefix;
        newPrefix = prefix + alphabet[i];
        processAllKmersRec(newPrefix, k-1);
    }
}

void kmerProcessor::processAllKmers(std::ofstream& ost)
{
    processAllKmersRec("", _k);
    for (int i = 0; i < freq_arr.size()-1; ++i)
    {
        ost << freq_arr[i] << " ";
    } 
    ost << freq_arr[freq_arr.size()-1] << '\n';
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

    string text;
    int k;
    getline(ist, text);
    ist >> k;

    kmerProcessor processor(text, k);
    processor.processAllKmers(ost);
    return 0;
}