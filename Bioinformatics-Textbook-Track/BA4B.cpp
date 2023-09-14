#include <iostream>
#include <fstream>
#include <string>
#include <map>

using std::map;
using std::string;
const int CODON_NUM = 64;
const char* CODON_TABLE = "dna_codons.lst";
const char* INFILENAME = "rosalind_ba4b.txt";
const char* OUTFILENAME = "out.txt";

const map<string, string> Complement = {
    {"A", "T"},
    {"T", "A"},
    {"C", "G"},
    {"G", "C"}
};

string RevCompl(string const& dna)
{
    string result;
    result.reserve(dna.size());
    for (int i = dna.size()-1; i >= 0; --i)
    {
        result += Complement.at(dna.substr(i, 1));
    }
    return result;
}

void readTableCodonToAA(std::ifstream& tablefile, map<string, string>& table)
{
    for (int i = 0; i < CODON_NUM; ++i)
    {
        string key, value;
        tablefile >> key >> value;
        table[key] = value;
    }
}

string translate(string const& rna, map<string, string> const& codon_table)
{
    string result;
    result.reserve(rna.size() / 3);
    for (int i = 0; i <= (rna.size()-3); i += 3)
    {
        string codon = rna.substr(i, 3);
        auto it = codon_table.find(codon);
        if (it == codon_table.end())
        { // unknown codon
            std::cerr << "Warning: unknown codon: " << codon << " ";
            std::cerr << "at position " << i << ", exiting";
            exit(1);
        }
        if (it->second == "Stop")
            break;
        result += it->second;
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
    std::ifstream codonTableFile{CODON_TABLE};
    if (!ist) 
    {
        std::cerr << "Cannot open input codon table file!\n";
        exit(1);
    }
    std::ofstream ost{OUTFILENAME};
    if (!ost) 
    {
        std::cerr << "Cannot open output file!\n";
        exit(1);
    }

    map<string, string> codonTable;
    readTableCodonToAA(codonTableFile, codonTable);
    string dna, peptide;
    getline(ist, dna);
    getline(ist, peptide);

    int k = peptide.size() * 3;
    for (int i = 0; i <= dna.size()-k; ++i)
    {
        string pat = dna.substr(i, k);
        string revc = RevCompl(pat);
        if ((translate(pat, codonTable) == peptide) || (translate(revc, codonTable) == peptide))
        {
            ost << pat << '\n';
        }
    }

    return 0;
}