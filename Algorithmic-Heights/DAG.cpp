#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <stack>

const char* INFILENAME = "rosalind_dag.txt";
const char* OUTFILENAME = "out.txt";
bool isDirected = true;
using adjList = std::vector<std::list<uint>>;

class Graph
{
    uint E;
    adjList nodes;
    std::vector<bool> visited;
    std::vector<bool> on_stack;

public:
    Graph(uint vNum)
    {
        this->E  = 0;
        this->nodes = adjList(vNum);
        this->visited = std::vector<bool>(vNum, false);
        this->on_stack = std::vector<bool>(vNum, false);
    }
    uint vNum() const { return this->nodes.size(); }
    uint eNum() const { return this->E; }

    void addEdge(uint start, uint end)
    {
        this->nodes[start].push_back(end);
        ++(this->E);
        if (!isDirected)
        {
            this->nodes[end].push_back(start);
        }
    }
    uint isAcyclic();
};

uint Graph::isAcyclic()
{
    std::stack<uint> stack;
    for (uint w = 0; w < this->vNum(); ++w)
    {   // in case of several components
        if (this->visited[w])
            continue;
        
        stack.push(w);
        while(!stack.empty())
        {
            uint s = stack.top();
            if (!this->visited[s])
            {
                this->visited[s] = true;
                this->on_stack[s] = true;
            } else {
                this->on_stack[s] = false;
                stack.pop();
            }

            for (auto v: this->nodes[s])
            {
                if (!this->visited[v])
                {
                    stack.push(v);
                } else if (this->on_stack[v]) {
                    return -1;
                }
            }
        }
    }
    return 1;
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
    return graph.isAcyclic();
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