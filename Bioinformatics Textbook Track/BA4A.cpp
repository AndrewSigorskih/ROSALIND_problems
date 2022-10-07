#include <iostream>
#include <fstream>
#include <string>
#include <map>

using std::map;
using std::string;
const int CODON_NUM = 64;
const char* CODON_TABLE = "rna_codons.lst";
const char* INFILENAME = "rosalind_ba4a.txt";
const char* OUTFILENAME = "out.txt";

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
    string rna;
    getline(ist, rna);
    ost << translate(rna, codonTable) << '\n';

    return 0;
}