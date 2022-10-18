#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <map>

struct AdjNode
{
    int end;
    char aa;
};

using std::map;
using std::vector;
using std::string;
using adjmap = map<int, vector<AdjNode>>;
const char* INFILENAME = "rosalind_ba11a.txt";
const char* MASSESFILENAME = "masses.lst";
const char* OUTFILENAME = "out.txt";

class SpectrumGraph
{
    adjmap nodes;
    vector<int> spectrum;
    map<int, char> masses;
    void parseSpectrum(std::ifstream&);
    void loadMasses();
public:
    SpectrumGraph(std::ifstream&);
    void printGraph(std::ofstream&);
};

SpectrumGraph::SpectrumGraph(std::ifstream& ist)
{
    loadMasses();
    parseSpectrum(ist);
    for (int i = 0; i < spectrum.size(); ++i)
    {
        for (int j = i+1; j < spectrum.size(); ++j)
        {
            int diff = spectrum[j] - spectrum[i];
            auto it = masses.find(diff);
            if (it != masses.end()) 
            {
                nodes[spectrum[i]].push_back({spectrum[j], it->second});
            }
        }
    }
}

void SpectrumGraph::parseSpectrum(std::ifstream& ist)
{
    spectrum.clear();
    spectrum.push_back(0); // first aa diff
    string temp;
    getline(ist, temp);
    std::stringstream ss(temp);

    int val;
    while (ss >> val)
        spectrum.push_back(val);
}

void SpectrumGraph::loadMasses()
{
    std::ifstream source{MASSESFILENAME};
    if (!source) 
    {
        std::cerr << "Cannot open file with AA masses!\n";
        exit(1);
    }
    char aa;
    double mass;
    while (source >> aa >> mass)
    {   // possible disign flaw: 
        // same-mass amino acids overwrite each other
        masses[static_cast<int>(mass)] = aa;
    }
}

void SpectrumGraph::printGraph(std::ofstream& ost)
{
    for (auto node: nodes)
    {
        for (AdjNode edge: node.second)
        {
            ost << node.first << "->" << edge.end << ":" << edge.aa << "\n";
        }
    }
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
    
    SpectrumGraph graph(ist);
    graph.printGraph(ost);

    return 0;
}