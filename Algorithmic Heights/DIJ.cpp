#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <utility> //pair
#include <set> // imitating proirity queue
#include <climits> // UINT_MAX

const char* INFILENAME = "rosalind_dij.txt";
const char* OUTFILENAME = "out.txt";
bool isDirected = true;

struct adjNode
{
    unsigned name;
    unsigned weight; // dijkstra operates positive edge weights only
};
using adjList = std::vector<std::list<adjNode>>;

class Graph
{
    unsigned edgesNum = 0;
    adjList nodes;
    std::vector<bool> visited;
    std::vector<unsigned> path;

public:
    unsigned vNum() const { return nodes.size(); }
    unsigned eNum() const { return this->edgesNum; }

    Graph(unsigned verticesNum) 
    {
        this->nodes = adjList(verticesNum);
        this->visited = std::vector<bool>(verticesNum, false);
        this->path = std::vector<unsigned>(verticesNum, UINT_MAX);
    }

    void addEdge(unsigned start, unsigned end, unsigned weight)
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

    void Dijkstra(unsigned s); 
    void printPath(std::ofstream& ofs);
};

void Graph::Dijkstra(unsigned s)
{
    //init
    this->path[s] = 0;
    std::set< std::pair<unsigned, unsigned> > queue;
    queue.insert(std::make_pair(0, s));

    while (!queue.empty())
    {
        std::pair<unsigned, unsigned> tmp = *(queue.begin());
        queue.erase(queue.begin());
        unsigned u = tmp.second;
        //iterate all neighbours
        for (auto it = this->nodes[u].begin(); it != this->nodes[u].end(); ++it)
        {
            unsigned v, weight;
            v = (*it).name;
            weight = (*it).weight;

            if (this->path[v] > this->path[u] + weight)
            {
                if (this->path[v] != UINT_MAX)
                    queue.erase(queue.find(std::make_pair(this->path[v], v)));
                this->path[v] = this->path[u] + weight;
                queue.insert(std::make_pair(this->path[v], v));
            }
        }
    }
    
}

void Graph::printPath(std::ofstream& ofs)
{
    for (int i = 0; i < this->vNum(); ++i)
    {
        auto val = this->path[i];
        if (val == UINT_MAX) { ofs << -1; }
        else { ofs << val; }
        if (i < (this->vNum()-1))
            ofs << " ";
    }
    //ofs << '\n';
}

int main()
{
    std::ifstream ist{INFILENAME};
    if (!ist) 
    {
        std::cout << "Cannot open input file!\n";
        exit(1);
    }

    unsigned x, y;
    
    ist >> x >> y;
    Graph graph(x);
    for (int i=0; i < y; ++i)
    {
        unsigned x, y, w;
        ist >> x >> y >> w;
        graph.addEdge(x-1, y-1, w);
    }

    graph.Dijkstra(0);

    std::ofstream ost{OUTFILENAME};
    if (!ost) 
    {
        std::cout << "Cannot open output file!\n";
        exit(1);
    }

    graph.printPath(ost);
    
    return 0;
}