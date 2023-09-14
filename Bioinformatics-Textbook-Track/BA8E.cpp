#include <iostream>
#include <iomanip> // precision
#include <fstream>
#include <vector>
#include <unordered_map>
#include <utility>
#include <algorithm> // min element
/* Hierarchical Clustering (basically UPGMA, but clusters only, no distances) */
using std::vector;
using std::unordered_map;
const char* INFILENAME = "rosalind_ba8e.txt";
const char* OUTFILENAME = "out.txt";

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;  
    }
};
template <class Key, class Value, class HashFn>
std::pair<Key, Value> findMinValuePair(
    std::unordered_map<Key, Value, HashFn> &x)
{
    return *std::min_element(x.begin(), x.end(),
                             [](const std::pair<Key, Value> &p1,
                                const std::pair<Key, Value> &p2)
                             {
                                 return p1.second < p2.second;
                             });
}

struct Node
{
    vector<int> neighbors;
    int clusterSize;
    Node() : clusterSize(1) {};
    Node(int size) : clusterSize(size) {};
};

using adjlist = vector<Node>;
using key = std::pair<int, int>;
using distMatrix = unordered_map<key, double, pair_hash>;

inline key orderPair(int const a, int const b)
{
    return (a < b) ? key(a, b) : key(b, a);
}

class Graph
{
    int nodesNum, leavesNum;
    adjlist nodes;
    distMatrix matrix;
    vector<int> clusters;
    void addNode(int i, int j, double w);
    void addEdge(int start, int end, double weight);
    double calcDistance(int, int, int);
    void getDescendants(int idx, vector<int>& cluster);
public:
    Graph(std::ifstream&);
    void UPGMA_Clusters();
    friend std::ofstream& operator << (std::ofstream& ost, Graph& obj);
};

Graph::Graph(std::ifstream& ist)
{
    ist >> nodesNum;
    leavesNum = nodesNum;
    nodes.reserve(nodesNum);
    clusters.reserve(nodesNum);
    for (int i = 0; i < nodesNum; ++i)
    {
        nodes.push_back(Node());
        clusters.push_back(i);
        for (int j = 0; j < nodesNum; ++j)
        {
            double val;
            ist >> val;
            if (i >= j)
                continue;
            matrix[orderPair(i, j)] = val;
        }
    }
}

inline double Graph::calcDistance(int cluster, int x, int y)
{   // (|X| * dxw + |Y| * dyw) / (|X| + |Y|)
    return (nodes[x].clusterSize * matrix[orderPair(x, cluster)] + \
             nodes[y].clusterSize * matrix[orderPair(y, cluster)] ) / \
             (nodes[x].clusterSize + nodes[y].clusterSize);
}

inline void Graph::addEdge(int start, int end, double weight)
{
    nodes[start].neighbors.push_back(end);
}

void Graph::addNode(int i, int j, double w)
{
    // merge Ci and Cj into a new cluster Cnew with |Ci| + |Cj| elements
    Node newNode(nodes[i].clusterSize + nodes[j].clusterSize);
    // add a new node labeled by cluster Cnew to T
    nodes.push_back(newNode);
    // connect node Cnew to Ci and Cj by directed edges
    addEdge(nodesNum, i, w);
    addEdge(nodesNum, j, w);
}

void Graph::UPGMA_Clusters()
{
    while (matrix.size() > 1)
    {
        // find the two closest clusters Ci and Cj
        auto min_el = findMinValuePair(matrix);
        int i = min_el.first.first, j = min_el.first.second;
        double w = min_el.second / 2;
        // create new node and add to tree
        addNode(i, j, w);
        // add a row/column to D for Cnew by computing D(Cnew, C) for each C in Clusters
        for (auto cluster: clusters)
        {
            if ((cluster != i) && (cluster != j))
            {
                matrix[orderPair(cluster, nodesNum)] = calcDistance(cluster, i, j);
            }
        }
        // remove the rows and columns of D corresponding to Ci and Cj
        for (auto cluster: clusters)
        {
            if ((cluster != i) && (cluster != j))
            {
                matrix.erase(orderPair(i, cluster));
                matrix.erase(orderPair(j, cluster));
            }
        }
        matrix.erase(orderPair(i, j));
        // remove Ci and Cj from Clusters and add new node to clusters
        clusters.erase(std::remove(clusters.begin(), clusters.end(), i), clusters.end());
        clusters.erase(std::remove(clusters.begin(), clusters.end(), j), clusters.end());
        clusters.push_back(nodesNum);
        // increment index
        ++nodesNum;
    }
    // add root
    int i = clusters[0], j = clusters[1];
    double w = matrix[orderPair(i, j)] / 2;
    addNode(i, j, w);
}

void Graph::getDescendants(int idx, vector<int>& cluster)
{
    if (nodes[idx].clusterSize == 1)
    {
        cluster.push_back(idx);
        return;
    }
    for (auto nei: nodes[idx].neighbors)
        getDescendants(nei, cluster);
}

std::ofstream& operator << (std::ofstream& ost, Graph& obj)
{
    for (int i = obj.leavesNum; i <= obj.nodesNum; ++i)
    {
        vector<int> cluster;
        cluster.reserve(obj.nodes[i].clusterSize);
        obj.getDescendants(i, cluster);
        for (auto item: cluster)
            ost << item + 1 << " ";
        ost << '\n';
    }
    return ost;
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
    ost << std::fixed;
    ost << std::setprecision(3);

    Graph graph(ist);
    graph.UPGMA_Clusters();
    ost << graph;
    return 0;
}