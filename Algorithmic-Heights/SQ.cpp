#include <iostream>
#include <fstream>
#include <vector>
using std::vector;
using adjMat = vector<vector<uint>>;

const char* INFILENAME = "rosalind_sq.txt";
const char* OUTFILENAME = "out.txt";
const int CycleSize = 4;
const bool isDirected = false;

class Graph
{
    adjMat mat;
    uint V, E;
    vector<bool> marked;
public:
    Graph(uint vNum) : V(vNum), E(0), mat(vNum), marked(vNum, false)
    {
        for (int i = 0; i < vNum; ++i)
        {
            mat[i].resize(vNum);
            std::fill(mat[i].begin(), mat[i].end(), 0);
        }
    }
    int vNum() const { return this->V; }
    int eNum() const { return this->E; }
    void addEdge(uint, uint);

    void DFS(int, uint, uint, int&);
    int checkCycles(int);
};

void Graph::addEdge(uint start, uint end)
{
    ++E;
    mat[start][end] = 1;
    if (!isDirected)
    {
        mat[end][start] = 1;
    }
}

void Graph::DFS(int n, uint vert, uint start, int& count)
{
    this->marked[vert] = true;
    if (n == 0)
    {
        this->marked[vert] = false;
        if (this->mat[vert][start])
        {
            ++count;
            return;
        } else {
            return;
        }
    }

    for (uint i = 0; i < this->V; ++i)
    {
        if ((!this->marked[i]) && (this->mat[vert][i]))
        {
            this->DFS( n-1, i, start, count);
        }
    }
    this->marked[vert] = false;
}

int Graph::checkCycles(int n)
{
    std::fill(marked.begin(), marked.end(), false);
    int count = 0;
    for (uint i = 0; i < this->V - (n-1); ++i)
    {
        this->DFS(n-1, i, i, count);
        this->marked[i] = true;
        if (count > 0)
            return 1;
    }
    return -1;
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
    return graph.checkCycles(CycleSize);
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
