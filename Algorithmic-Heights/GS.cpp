#include <iostream>
#include <fstream>
#include <vector>
#include <list>
using std::vector;
using adjList = std::vector<std::list<uint>>;
const char* INFILENAME = "rosalind_gs.txt";
const char* OUTFILENAME = "out.txt";
const bool isDirected = true;
// mother vertex of a graph:
// https://www.geeksforgeeks.org/find-a-mother-vertex-in-a-graph/
class Graph
{
    uint V, E;
    adjList nodes;
    vector<bool> visited;
    void DFSUtil(uint v);

public:
    Graph(uint vNum) : V(vNum), E(0), nodes(vNum), 
                    visited(vNum, false)
    { }
    uint vNum() const { return this->V; }
    uint eNum() const { return this->E; }
    bool isVisited(uint v) const { return visited[v]; }
    void wash() { std::fill(visited.begin(), visited.end(), false); }
    void addEdge(uint, uint);
    int findMother();
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

void Graph::DFSUtil(uint v)
{
    visited[v] = true;
    for (auto u: nodes[v])
        if (!visited[u])
            DFSUtil(u);
}

int Graph::findMother()
{
    uint v = 0;
    //this->wash();
    for(uint i = 0; i < V; ++i)
    {
        if (!visited[i])
        {
            DFSUtil(i);
            v = i;
        }
    }
    this->wash();
    DFSUtil(v);
    for (uint i = 0; i < V; ++i)
        if (!visited[i])
            return -1;
    return v + 1;
}

int check_graph(std::ifstream& ist)
{
    uint x, y;
    ist >> x >> y;
    Graph graph(x);
    for (int i = 0; i < y; ++i)
    {
        uint x, y;
        ist >> x >> y;
        graph.addEdge(x-1, y-1);
    }
    return graph.findMother();
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