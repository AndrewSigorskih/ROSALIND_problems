#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <stack>

const char* INFILENAME = "rosalind_ts.txt";
const char* OUTFILENAME = "out.txt";
bool isDirected = true;
using adjList = std::vector<std::list<uint>>;

// for some reason testing system accepts only
// sortings obtained by recursive solutions

class Graph
{
    uint E;
    adjList nodes;
    std::vector<bool> visited;

public:
    Graph(uint vNum)
    {
        this->E  = 0;
        this->nodes = adjList(vNum);
        this->visited = std::vector<bool>(vNum, false);
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
    void topoSortUtil(uint v, std::stack<uint>& stack);
    std::vector<uint> topoSort();
};

void Graph::topoSortUtil(uint v, std::stack<uint>& stack)
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

std::vector<uint> Graph::topoSort()
{
    std::stack<uint> stack;
    for(uint v = 0; v < this->vNum(); ++v)
    {
        if (!this->visited[v])
        {
            this->topoSortUtil(v, stack);
        }
    }
    std::vector<uint> res;
    res.reserve(this->vNum());
    while(!stack.empty())
    {
        int v = stack.top();
        stack.pop();
        res.push_back(v);
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
    
    std::vector<uint> res = graph.topoSort();
    for (auto i: res)
    {
        ost << (i+1) << " ";
    }
    ost << '\n';
    return 0;
}


/*
std::vector<uint> Graph::topoSort()
{
    std::stack<uint> stack;
    std::stack<uint> out;

    for (int i = 0; i < this->vNum(); ++i)
    {
        if (this->visited[i])
            continue;

        stack.push(i);
        while (!stack.empty())
        {
            uint v = stack.top();
            stack.pop();
            if (this->visited[v])
            {
                out.push(v);
                continue;
            }
            this->visited[v] = true;
            stack.push(v);
            for (auto j: this->nodes[v])
            {
                if (!this->visited[j])
                {
                    stack.push(j);
                }
            }
        }
    }

    std::vector<uint> res;
    res.reserve(this->vNum());
    while(!out.empty())
    {
        int v = out.top();
        out.pop();
        res.push_back(v);
    }
    return res;
}
*/