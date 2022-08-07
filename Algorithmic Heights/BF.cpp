#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <climits> // INT_MAX

const char* INFILENAME = "rosalind_bf.txt";
const char* OUTFILENAME = "out.txt";
bool isDirected = true;

struct Edge {
    int start, end, weight;
};

class Graph
{
    int V;
    int E = 0;
    std::vector<Edge> edges;
    std::vector<int> path;

public:
    Graph(int vNum)
    {
        this->V = vNum;
        this->E = 0;
        this->edges = std::vector<Edge>(0);
        this->path = std::vector<int>(this->V, INT_MAX);
    }

    int vNum() const { return this->V; }
    int eNum() const { return this->E; }

    void addEdge(int, int, int);
    void printPath(std::ofstream& ofs);
    void BF(int s);
};

void Graph::addEdge(int start, int end, int w)
{
    this->edges.push_back(Edge{start, end, w});
    ++(this->E);
}

void Graph::printPath(std::ofstream& ofs)
{
    for (int i = 0; i < this->vNum(); ++i)
    {
        auto val = this->path[i];
        if (val == INT_MAX) { ofs << "x"; }
        else { ofs << val; }
        if (i < (this->vNum()-1))
            ofs << " ";
    }
    ofs << '\n';
}

void Graph::BF(int s)
{
    //init
    this->path[s] = 0;
    // Relax all edges |V| - 1 times.
    for (int i = 1; i <= V - 1; ++i)
    {
        for (int j = 0; j < this->E; ++j)
        {
            int u = this->edges[j].start;
            int v = this->edges[j].end;
            int w = this->edges[j].weight;
            if ((this->path[u] != INT_MAX) 
                && (this->path[v] > this->path[u] + w))
                this->path[v] = this->path[u] + w;
        }
    }
}


int main()
{
    std::ifstream ist{INFILENAME};
    if (!ist) 
    {
        std::cout << "Cannot open input file!\n";
        exit(1);
    }

    int x, y;
    
    ist >> x >> y;
    Graph graph(x);
    for (int i=0; i < y; ++i)
    {
        int x, y, w;
        ist >> x >> y >> w;
        graph.addEdge(x-1, y-1, w);
    }

    graph.BF(0);

    std::ofstream ost{OUTFILENAME};
    if (!ost) 
    {
        std::cout << "Cannot open output file!\n";
        exit(1);
    }

    graph.printPath(ost);

    return 0;
}
