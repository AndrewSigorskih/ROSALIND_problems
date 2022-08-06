#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <iterator>

const char* INFILENAME = "rosalind_ddeg.txt";
const char* OUTFILENAME = "out.txt";
bool isDirected = false;

struct adjNode
{
    unsigned name;
    double weight;
};
using adjList = std::vector<std::list<adjNode>>;


class Graph
{
    unsigned edgesNum = 0;
    adjList nodes;

public:
    unsigned vNum() const { return nodes.size(); }
    unsigned eNum() const { return this->edgesNum; }

    Graph(unsigned verticesNum) 
    {
        this->nodes = adjList(verticesNum);
    }

    void addEdge(unsigned start, unsigned end, double weight)
    {
        adjNode node {end, weight};
        this->nodes[start].push_back(node);
        ++(this->edgesNum);
        if (!isDirected)
        {
            adjNode node {start, weight};
            this->nodes[end].push_back(node);
        }
    }

    unsigned get_ith_num_neighbors(unsigned v) const
    { return this->nodes[v].size(); } // mostly for debugging

    unsigned getNeiDegSum(unsigned v) const
    {
        unsigned res = 0;
        auto it  = nodes[v].begin();
        while (it != nodes[v].end())
        {
            unsigned nodeName = (*it).name;
            res += this->nodes[nodeName].size();
            it = std::next(it);
        }
        return res;
    }
};

int main()
{
    std::ifstream ist{INFILENAME};
    if (!ist) 
    {
	std::cout << "Cannot open input file!\n";
	exit(1);
    }

    double W = 0;
    int x, y;
    
    ist >> x >> y;
    Graph graph(x);
    for (int i=0; i < y; ++i)
    {
        int x, y;
        ist >> x >> y;
        graph.addEdge(x-1, y-1, W);
    }

    std::ofstream ost{OUTFILENAME};
    if (!ost) 
    {
	std::cout << "Cannot open output file!\n";
	exit(1);
    }

    for (int i=0; i < x; ++i)
    {
        ost << graph.getNeiDegSum(i);
        if ( i < x-1)
            ost << " ";
    }
    ost << '\n';
    
    return 0;
}
