#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <stack>
using std::vector;
using adjList = std::vector<std::list<uint>>;
const char* INFILENAME = "rosalind_scc.txt";
const char* OUTFILENAME = "out.txt";
const bool isDirected = true;
// Kosaraju’s algorithm:
// https://www.geeksforgeeks.org/strongly-connected-components/
class Graph
{
    uint V, E;
    adjList nodes;
    vector<bool> visited;
    void fillOrder(uint v, std::stack<uint>& stack);
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
    Graph getTranspose();
    int countSCCs();
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
    {
        if (!visited[u])
        {
            DFSUtil(u);
        }
    }
}

Graph Graph::getTranspose()
{
    Graph graph(V);
    for (uint v = 0; v < V; ++v)
    {
        for (auto u: nodes[v])
        {
            graph.addEdge(u, v);
        }
    }
    return graph;
}

void Graph::fillOrder(uint v, std::stack<uint> &stack)
{
    visited[v] = true;
    for (auto u: nodes[v])
    {
        if (!visited[u])
        {
            fillOrder(u, stack);
        }
    }
    stack.push(v);
}

int Graph::countSCCs()
{
    int res = 0;
    std::stack<uint> stack;
    this->wash();

    for(uint i = 0; i < V; ++i)
        if (!visited[i])
            fillOrder(i, stack);
    
    Graph rev = this->getTranspose();
    this->wash();

    while(!stack.empty())
    {
        uint v = stack.top();
        stack.pop();
        if (!rev.isVisited(v))
        {
            rev.DFSUtil(v);
            ++res;
        }
    }
    return res;
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

    uint x, y;
    ist >> x >> y;
    Graph graph(x);
    for (int i = 0; i < y; ++i)
    {
        uint x, y;
        ist >> x >> y;
        graph.addEdge(x-1, y-1);
    }

    ost << graph.countSCCs() << '\n';
    return 0;
}