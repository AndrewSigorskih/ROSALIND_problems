#include <iostream>
#include <fstream>
#include <vector>
//#include <list>
#include <unordered_map>
#include <stack>
#include <climits> // INT_MAX
// for DAG shortest paths can be found in O(V+E) time
// topo sort + Bellman-Ford in linearized order
// https://www.geeksforgeeks.org/shortest-path-for-directed-acyclic-graphs/
// BEWARE! Edges are given multiple times in this dataset 
// hence hashtable instead of a list.
// only the last weight given for each edge should be considered valid.
const char* INFILENAME = "rosalind_sdag.txt";
const char* OUTFILENAME = "out.txt";
const bool isDirected = true;
using std::vector;
using std::stack;

struct adjNode
{
    uint name;
    int weight; 
};
//using adjList = vector<std::list<adjNode>>;
using adjList = vector<std::unordered_map<uint, adjNode>>;
class Graph
{
    uint V, E;
    adjList nodes;
    vector<bool> visited;
    void topoSortUtil(uint v, stack<uint>& stack);
public:
    Graph(uint vNum) : V(vNum), E(0), nodes(vNum),
                        visited(vNum, false)
    {  }
    uint vNum() const { return V; }
    uint eNum() const { return E; }
    bool isVisited(uint v) const { return visited[v]; }
    void wash() { std::fill(visited.begin(), visited.end(), false); }
    void addEdge(uint, uint, int);
    void shortestPath(uint, std::ofstream&);
};

void Graph::addEdge(uint start, uint end, int weight)
{
    ++E;
    adjNode node{end, weight};
    auto it = nodes[start].find(end);
    if (it == nodes[start].end())
    {
        nodes[start][end] = node;
    } else {
        nodes[start].erase(it);
        nodes[start][end] = node;
    }
    //directed-only for now
}

void Graph::topoSortUtil(uint v, stack<uint>& stack)
{
    visited[v] = true;
    for (auto it: nodes[v])
    {   
        adjNode i = it.second;
        if (!visited[i.name])
        {
            topoSortUtil(i.name, stack);
        }
    }
    stack.push(v);
}

void Graph::shortestPath(uint s, std::ofstream& ost)
{
    stack<uint> stack;
    vector<int> path(V, INT_MAX);
    //this->wash();
    for (uint i = 0; i < V; ++i)
        if (!visited[i])
            topoSortUtil(i, stack);
    
    path[s] = 0;
    while(!stack.empty())
    {
        uint u = stack.top();
        stack.pop();
        if (path[u] != INT_MAX)
        {
            for (auto it: nodes[u])
            {
                adjNode i = it.second;
                if (path[i.name] > path[u] + i.weight)
                {
                    path[i.name] = path[u] + i.weight;
                }
            }
        }
    }
    for (uint i = 0; i < V; ++i)
    {
        (path[i] == INT_MAX)? ost << "x ": ost << path[i] << " ";
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
    uint x, y;
    ist >> x >> y;
    Graph graph(x);
    for (int i = 0; i < y; ++i)
    {
        uint x, y;
        int w;
        ist >> x >> y >> w;
        graph.addEdge(x-1, y-1, w);
    }
    graph.shortestPath(0, ost);
    return 0;
}