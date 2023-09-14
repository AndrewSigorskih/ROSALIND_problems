#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <climits>
#include <queue>
using std::list;
using std::vector;
using std::priority_queue;
const char* INFILENAME = "rosalind_dij.txt";
const char* OUTFILENAME = "out.txt";
const bool isDirected = true;
// https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-using-priority_queue-stl/
struct adjNode
{
    uint name;
    uint weight; // dijkstra operates positive edge weights only
};
using adjList = vector<list<adjNode>>;

bool Compare(const adjNode& a, const adjNode& b)
{
    if (a.weight != b.weight)
        return (a.weight < b.weight);
    else
        return (a.name < b.name);
}

class Graph
{
    uint V, E;
    adjList nodes;
    vector<uint> path;
public:
    Graph(uint vNum) : V(vNum), E(0), nodes(vNum),
                        path(vNum, UINT_MAX)
    {  }
    uint vNum() const { return V; }
    uint eNum() const { return E; }
    uint pathTo(uint v) const { return path[v]; }
    void initPath() { std::fill(path.begin(), path.end(), UINT_MAX); }
    void addEdge(uint, uint, uint);
    void Dijkstra(uint); 
};

void Graph::addEdge(uint start, uint end, uint weight)
{
    ++E;
    adjNode node {end, weight};
    nodes[start].push_back(node);
    if (!isDirected)
    {
        adjNode node {start, weight};
        nodes[end].push_back(node);
    }
}

void Graph::Dijkstra(uint s)
{   
    //initPath();
    priority_queue<adjNode, vector<adjNode>, decltype(&Compare)> queue(Compare);
    queue.push(adjNode{s, 0});
    path[s] = 0;
    while(!queue.empty())
    {
        uint u = queue.top().name;
        queue.pop();
        for (auto v: nodes[u])
        {
            if (path[v.name] > path[u] + v.weight)
            {
                path[v.name] = path[u] + v.weight;
                queue.push(adjNode{v.name, path[v.name]});
            }
        }
    }
}

void check_graph(std::ifstream& ist, std::ofstream& ost)
{
    uint x, y;
    ist >> x >> y;
    Graph graph(x);

    for (int i = 0; i < y; ++i)
    {
        uint x, y, w;
        ist >> x >> y >> w;
        graph.addEdge(x-1, y-1, w);

    }
    graph.Dijkstra(0);
    for (uint i = 0; i < graph.vNum(); ++i)
    {
        uint dist = graph.pathTo(i);
        if (dist < UINT_MAX)
            ost << dist;
        else
            ost << -1;
        if (i < graph.vNum() -1) 
            ost << " ";
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
    check_graph(ist, ost);
    return 0;
}
