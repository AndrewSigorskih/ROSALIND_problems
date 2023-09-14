#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <stack>

const char* INFILENAME = "rosalind_hdag.txt";
const char* OUTFILENAME = "out.txt";
const bool isDirected = true;
using std::vector;
using std::stack;
using adjList = std::vector<std::list<uint>>;

// for DAG this problem is solvable in O(n + m)
// https://stackoverflow.com/questions/16124844/algorithm-for-finding-a-hamiltonian-path-in-a-dag

class Graph
{
    uint V, E;
    adjList nodes;
    std::vector<bool> visited;
    void topoSortUtil(uint v, stack<uint>& stack);

public:
    Graph(uint vNum) : V(vNum), E(0), nodes(vNum),
                        visited(vNum, false)
    {  }
    uint vNum() const { return V; }
    uint eNum() const { return E; }
    bool isVisited(uint v) const { return visited[v]; }
    void wash() { std::fill(visited.begin(), visited.end(), false); }
    void addEdge(uint, uint);
    vector<uint> topoSort();
    bool isNeighbor(uint v, uint next);
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

void Graph::topoSortUtil(uint v, stack<uint>& stack)
{
    this->visited[v] = true;
    for (auto i: this->nodes[v])
    {
        if (!this->visited[i])
        {
            this->topoSortUtil(i, stack);
        }
    }
    stack.push(v);
}

vector<uint> Graph::topoSort()
{
    stack<uint> stack;
    for(uint v = 0; v < this->vNum(); ++v)
    {
        if (!this->visited[v])
        {
            this->topoSortUtil(v, stack);
        }
    }

    vector<uint> res;
    res.reserve(this->vNum());
    while(!stack.empty())
    {
        int v = stack.top();
        stack.pop();
        res.push_back(v);
    }
    return res;
}

bool Graph::isNeighbor(uint v, uint next)
{
    bool isNei = false;
    for (auto u: nodes[v])
    {
        if (next == u)
        {
            isNei = true;
            break;
        }
    }
    return isNei;
}

void check_graph(std::ifstream& ist, std::ofstream& ost)
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

    bool isHam = true;
    vector<uint> path = graph.topoSort();
    for (uint i = 0; i < path.size()-1; ++i)
    {
        isHam &= graph.isNeighbor(path[i], path[i+1]);
        if (!isHam) 
        {
            ost << -1 << '\n';
            return;
        }
    }
    ost << 1;
    for (auto i: path)
    {
        ost << " " << i+1;
    }
    ost << '\n';
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
        check_graph(ist, ost);
    }
    return 0;
}
