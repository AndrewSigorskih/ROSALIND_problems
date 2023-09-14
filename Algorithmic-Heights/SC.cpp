#include <iostream>
#include <fstream>
#include <vector>
#include <list>
using std::vector;
using adjMat = vector<vector<uint>>;
const char* INFILENAME = "rosalind_sc.txt";
const char* OUTFILENAME = "out.txt";
const bool isDirected = true;
// Warshall's algorithm
// https://stackoverflow.com/questions/17860272/computation-of-path-matrix-from-the-adjacency-matrix
// https://www.geeksforgeeks.org/check-if-a-graph-is-strongly-unilaterally-or-weakly-connected/
class Graph
{
    uint V, E;
    adjMat mat;
    vector<bool> visited;

public:
    Graph(uint vNum) : V(vNum), E(0), mat(vNum), 
                    visited(vNum, false)
    {
        for (uint i = 0; i < vNum; ++i)
        {
            mat[i].resize(vNum);
            std::fill(mat[i].begin(), mat[i].end(), 0);
        }
    }
    uint vNum() const { return this->V; }
    uint eNum() const { return this->E; }
    bool isVisited(uint v) const { return visited[v]; }
    void wash() { std::fill(visited.begin(), visited.end(), false); }
    void addEdge(uint, uint);
    adjMat calcPathMatrix();
    int isSemiConnected();
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

adjMat Graph::calcPathMatrix()
{
    adjMat path(V);
    for (uint i = 0; i < V; ++i)
    {
        path[i].resize(V);
        for (uint j = 0; j < V; ++j) 
        {
			path[i][j] = mat[i][j];
		}
    }
    // O(n ^ 3) :(
    for (uint k = 0; k < V; ++k)
    {
        for (uint i = 0; i < V; ++i)
        {
            for (uint j = 0; j < V; ++j)
            {
                path[i][j] = ((path[i][j]) || 
                              ((path[i][k]) && (path[k][j])));
            }
        }
    }
    
    return path;
}

int Graph::isSemiConnected()
{
    adjMat path = calcPathMatrix();
    for (uint i = 0; i < V; ++i)
    {
        for (uint j = i; j < V; ++j)
        {
            if ((j != i) && !(path[i][j] || path[j][i]))
            {
                return -1;
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
    return graph.isSemiConnected();
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