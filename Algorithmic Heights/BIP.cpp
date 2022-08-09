#include <iostream>
#include <fstream>
#include <vector>
#include <list>
using std::vector;
using adjList = std::vector<std::list<uint>>;
const char* INFILENAME = "rosalind_bip.txt";
const char* OUTFILENAME = "out.txt";
const bool isDirected = false;

class Graph
{
    uint V, E;
    adjList nodes;
    vector<bool> visited;
    vector<bool> color;

public:
    Graph(uint vNum) : V(vNum), E(0), nodes(vNum), 
                    visited(vNum, false), color(vNum, false)
    { }
    uint vNum() const { return this->V; }
    uint eNum() const { return this->E; }
    bool isVisited(uint v) const { return visited[v]; }
    void addEdge(uint, uint);
    bool isBipartite(uint);
};

void Graph::addEdge(uint start, uint end)
{
    ++E;
    nodes[start].push_back(end);
    if (!isDirected)
    {
        nodes[end].push_back(start);
    }
}

bool Graph::isBipartite(uint v)
{
    for (uint u: nodes[v])
    {
        if (!visited[u])
        {
            visited[u] = true;
            color[u] = !color[v];
            if (!this->isBipartite(u))
                return false;
        } else if (color[u] == color[v])
            return false;
    }
    return true;
}

int check_graph(std::ifstream& ist)
{
    bool bip = true;
    uint x, y;
    ist >> x >> y;
    Graph graph(x);
    for (int i = 0; i < y; ++i)
    {
        uint x, y;
        ist >> x >> y;
        graph.addEdge(x-1, y-1);
    }
    for (uint v = 0; v < graph.vNum(); ++v)
    {
        if (!graph.isVisited(v))
        {
            bip &= graph.isBipartite(v);
            if (!bip)
            {
                return -1;
            }
        }
    }
    return 1;
}

int main()
{
    std::ifstream ist{INFILENAME};
    if (!ist) 
    {
        std::cout << "Cannot open input file!\n";
        exit(1);
    }

    std::ofstream ost{OUTFILENAME};
    if (!ost) 
    {
        std::cout << "Cannot open output file!\n";
        exit(1);
    }

    int k;
    ist >> k;
    for (int i = 1; i <= k; ++i)
    {
        ost << check_graph(ist);
        if (i < k)
            ost << " ";
    }
    ost << '\n';

    return 0;
}