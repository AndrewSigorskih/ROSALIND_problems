#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <stack>
using std::vector;
using std::stack;
using adjList = vector<std::list<int>>;
const char* INFILENAME = "rosalind_2sat.txt";
const char* OUTFILENAME = "out.txt";
const bool isDirected = true;
// https://cp-algorithms.com/graph/2SAT.html
// https://www.geeksforgeeks.org/2-satisfiability-2-sat-problem/
class Graph
{
    int V, E;
    adjList adj;
    adjList adjInv;
    vector<bool> visited;
    vector<bool> visitedInv;
    vector<int> scc;
    void dfs1(int, stack<int>&);
    void dfs2(int, int);

public:
    Graph(int vNum) : V(vNum), E(0), adj(vNum), adjInv(vNum),
                      visited(vNum), visitedInv(vNum), scc(vNum)
    { }
    int vNum() const { return this->V; }
    int eNum() const { return this->E; }
    void addEdge(int, int);
    void addEdgeInv(int, int);
    void printGraph();
    void is2SAT(std::ofstream&);
};

void Graph::addEdge(int start, int end)
{
    ++E;
    adj[start].push_back(end);
}

void Graph::addEdgeInv(int start, int end)
{
    adjInv[start].push_back(end);
}

void Graph::printGraph()
{
    std:: cout << "V: " << V << " E: " << E << '\n';
    for (int i = 0; i < adj.size(); ++i)
    {
        std::cout << i << ": ";
        for (auto j: adj[i])
        {
            std::cout << j << ", ";
        }
        std::cout << '\n';
    }
    std::cout << '\n';
    for (int i = 0; i < adjInv.size(); ++i)
    {
        std::cout << i << ": ";
        for (auto j: adjInv[i])
        {
            std::cout << j << " ";
        }
        std::cout << '\n';
    }
}

void Graph::dfs1(int u, stack<int>& s)
{
    if (visited[u])
        return;
    visited[u] = true;
    for (auto i: adj[u])
        dfs1(i, s);
    s.push(u);
}

void Graph::dfs2(int u, int counter)
{
    if (visitedInv[u])
        return;
    visitedInv[u] = true;
    for (auto i: adjInv[u])
        dfs2(i, counter);
    scc[u] = counter;
}

void Graph::is2SAT(std::ofstream& ost)
{
    stack<int> s;
    int counter = 1; // number of SCCs
    // STEP 1 of Kosaraju's Algorithm
    for (int i = 0; i < V; ++i)
        if (!visited[i])
            dfs1(i, s);
    // STEP 2 of Kosaraju's Algorithm
    while(!s.empty())
    {
        int n = s.top();
        s.pop();
        if (!visitedInv[n])
        {
            dfs2(n, counter);
            ++counter;
        }
    }
    // analyze results
    vector<bool> result(V/2);
    for (int i = 0; i < V; i += 2)
    {
        if (scc[i] == scc[i+1]) // expression is unsatisfiable
        {
            ost << "0\n";
            return;
        }
        result[i/2] = (scc[i] > scc[i+1]);
    }
    // we good
    ost << "1 ";
    for (int i =0; i < result.size(); ++i)
    {
        if (result[i])
            ost << (i + 1) << " ";
        else
            ost << -(i+1) << " ";
    }
    ost << '\n';

}

void check_graph(std::ifstream& ist, std::ofstream& ost)
{
    int x, y;
    ist >> x >> y;
    Graph graph(2*x);
    for (int i = 0; i < y; ++i)
    {
        int a, b;
        ist >> a >> b;
        if ((a > 0) && (b > 0))
        {
            graph.addEdge(2*abs(a)-1, 2*abs(b)-2);
            graph.addEdge(2*abs(b)-1, 2*abs(a)-2);
            graph.addEdgeInv(2*abs(b)-2, 2*abs(a)-1);
            graph.addEdgeInv(2*abs(a)-2, 2*abs(b)-1);
        } else if ((a < 0) && (b > 0)) {
            graph.addEdge(2*abs(a)-2, 2*abs(b)-2);
            graph.addEdge(2*abs(b)-1, 2*abs(a)-1);
            graph.addEdgeInv(2*abs(b)-2, 2*abs(a)-2);
            graph.addEdgeInv(2*abs(a)-1, 2*abs(b)-1);
        } else if ((a > 0) && (b < 0)) {
            graph.addEdge(2*abs(a)-1, 2*abs(b)-1);
            graph.addEdge(2*abs(b)-2, 2*abs(a)-2);
            graph.addEdgeInv(2*abs(b)-1, 2*abs(a)-1);
            graph.addEdgeInv(2*abs(a)-2, 2*abs(b)-2);
        } else /*(a < 0) && (b < 0) */ {
            graph.addEdge(2*abs(a)-2, 2*abs(b)-1);
            graph.addEdge(2*abs(b)-2, 2*abs(a)-1);
            graph.addEdgeInv(2*abs(b)-1, 2*abs(a)-2);
            graph.addEdgeInv(2*abs(a)-1, 2*abs(b)-2);
        }
    }
    //graph.printGraph();
    //std::cout << "----------------------------------\n";

    graph.is2SAT(ost);
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
        check_graph(ist, ost);
    }
    ost << '\n';
    return 0;
}