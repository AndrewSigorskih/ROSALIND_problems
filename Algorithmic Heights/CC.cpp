#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <stack>

const char* INFILENAME = "rosalind_cc.txt";
const char* OUTFILENAME = "out.txt";
bool isDirected = false;
//using uint = unsigned int
using adjList = std::vector<std::list<uint>>;

class Graph
{
    uint E;
    adjList nodes;
    //std::vector<Color> visited;
    std::vector<bool> visited;

public:
    Graph(uint vNum)
    {
        this->nodes = adjList(vNum);
        //std::vector<Color> vec(vNum, WHITE);
        //this->visited = vec;
        this->visited = std::vector<bool>(vNum, false);
        this->E  = 0;
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

    void DFS(uint);
    uint CCnum();
};

void Graph::DFS(uint s)
{
    std::stack<uint> stack;
    stack.push(s);
    
    while(!stack.empty())
    {
        uint v = stack.top();
        this->visited[v] = true;
        stack.pop();
        for (auto it = this->nodes[v].begin();
            it != this->nodes[v].end(); ++it)
            {
                if (!this->visited[*it])
                {
                    stack.push(*it);
                }
            }
    }
}

uint Graph::CCnum()
{
    uint res = 0;
    std::fill(this->visited.begin(),
                this->visited.end(),
                false);
    for (uint v = 0; v < this->vNum(); ++v)
    {
        if (!this->visited[v])
        {
            this->DFS(v);
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
    for (int i=0; i < y; ++i)
    {
        uint x, y;
        ist >> x >> y;
        graph.addEdge(x-1, y-1);
    }

    ost << graph.CCnum() << '\n';

    return 0;
}

//enum Color {WHITE = 0, GREY, BLACK};
/*
void Graph::DFS(uint s)
{
    std::stack<uint> stack;
    stack.push(s);
    
    while(!stack.empty())
    {
        uint v = stack.top();
        //std::cout << "top node " << v + 1<< ' ';
        if (this->visited[v] == Color::WHITE)
        {
            //std::cout << " of color white, painted grey\n";
            this->visited[v] = Color::GREY;
            for (auto it = this->nodes[v].begin();
                      it != this->nodes[v].end(); ++it)
            {
                if (this->visited[*it] == Color::WHITE)
                {
                    //std::cout << "pushed white node " << *it + 1 << '\n';
                    stack.push(*it);
                }
            }
        } else if (this->visited[v] == Color::GREY) {
            //std::cout << " of color grey, painted black\n";
            this->visited[v] = Color::BLACK;
            stack.pop();
        } else {
            //std::cout << " of color black, popped\n";
            stack.pop();
        }
    }
}
*/