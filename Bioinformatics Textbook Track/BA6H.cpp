#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <vector>
#include <set>

using std::set;
using std::string;
using std::vector;
const char* INFILENAME = "rosalind_ba6h.txt";
const char* OUTFILENAME = "out.txt";

vector<int> parseChromosome(string const& buf)
{
    vector<string> tmp;
    std::istringstream iss(buf);
    std::copy(std::istream_iterator<string>(iss), 
              std::istream_iterator<string>(),
              std::back_inserter(tmp));
    vector<int> result(tmp.size());
    for (int i = 0; i < tmp.size(); ++i)
        result[i] = atoi(tmp[i].c_str());
	return result;
}

vector<string> parseGenome(string const& genome)
{
    vector<string> result;
    int i = 0, j = 0;
    while (i < genome.size())
    {
        if (genome[i] == '(')
        {
            j = 0;
            while ((genome[i+j] != ')') && (i+j < genome.size()))
                ++j;
            result.push_back(genome.substr(i+1, j-1));
            i += (j+1);
        } else {
            // wrong input format
            std::cerr << "Error: input doesn't match bracket formula!\n";
            exit(1);
        }
    }
    return result;
}

vector<int> ChromosomeToCycle(vector<int> const& chromosome)
{
    vector<int> result(chromosome.size()*2);
    for (int i = 0; i < chromosome.size(); ++i)
    {
        int val = chromosome[i];
        if (val > 0)
        {
            result[2*i] = 2*val - 1;
            result[2*i+1] = 2*val;
        } else {
            result[2*i] = -2*val;
            result[2*i+1] = -2*val - 1;
        }
    }
    return result;
}

set<string> ColoredEdges(string const& genome)
{
    set<string> result;
    for (string chromosome: parseGenome(genome))
    {
        vector<int> nodes = ChromosomeToCycle(parseChromosome(chromosome));
        for (int i = 0; i < nodes.size()/2; ++i)
        { // for each "pair", get second element and group it with 
          // first element from next pair. Assume array to be cyclic.
            result.insert(
                "(" + std::to_string(nodes[2*i + 1]) + ", " +
                std::to_string(nodes[(2*i + 2) % nodes.size()]) + ")"
                );
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

    string genome;
    getline(ist, genome);
    set<string> result = ColoredEdges(genome);
    int i = 0;
    for (auto item: result)
    {
        ost << item;
        if (i < result.size() - 1)
            ost << ", ";
        ++i;
    }
    
    return 0;
}