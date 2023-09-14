#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <iterator>

const char* INFILENAME = "rosalind_bfs.txt";
const char* OUTFILENAME = "out.txt";
bool isDirected = true;

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

    std::vector<int> BFS(unsigned s); 
};

std::vector<int> Graph::BFS(unsigned s)
{
    std::vector<int> res;
    res.resize(this->vNum(), -1);
    res[s] = 0;

    std::vector<bool> visited;
    visited.resize(this->vNum(), false);

    std::list<unsigned> queue;
    visited[s] = true;
    queue.push_back(s);
 
    while(!queue.empty())
    {
        s = queue.front();
        queue.pop_front();

        for (auto adjecent: this->nodes[s])
        {
            unsigned name = adjecent.name;
            if (!visited[name])
            {
                visited[name] = true;
                // parents path + 1 level, no W for now
                res[name] = res[s] + 1; 
                queue.push_back(name);
            }
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

    std::vector<int> res = graph.BFS(0);

    for (auto i: res)
    {
        ost << i << " ";
    }
    ost << '\n';
    
    return 0;
}