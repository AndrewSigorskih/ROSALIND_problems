#include "BA7.h"
const char* INFILENAME = "rosalind_ba7d.txt";
const char* OUTFILENAME = "out.txt";

// actual graph declaration starts here
struct Node
{
    vector<Edge> neighbors;
    int clusterSize;
    Node() : clusterSize(1) {};
    Node(int size) : clusterSize(size) {};
};

using adjlist = vector<Node>;

class Graph
{
    int nodesNum;
    adjlist nodes;
    distMatrix matrix;
    vector<int> clusters;
    unordered_map<int, double> ages; // "ages" of the nodes
    void addNode(int i, int j, double w);
    void addEdge(int start, int end, double weight);
    double calcDistance(int, int, int);
public:
    Graph(std::ifstream&);
    void UPGMA();
    void printGraph(std::ofstream&);
};

Graph::Graph(std::ifstream& ist)
{
    ist >> nodesNum;
    nodes.reserve(nodesNum);
    clusters.reserve(nodesNum);
    for (int i = 0; i < nodesNum; ++i)
    {
        nodes.push_back(Node());
        clusters.push_back(i);
        ages[i] = 0;
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

void Graph::addEdge(int start, int end, double weight)
{
    nodes[start].neighbors.push_back(Edge{end, weight});
    nodes[end].neighbors.push_back(Edge{start, weight});
}

inline double Graph::calcDistance(int cluster, int x, int y)
{   // (|X| * dxw + |Y| * dyw) / (|X| + |Y|)
    return (nodes[x].clusterSize * matrix[orderPair(x, cluster)] + \
             nodes[y].clusterSize * matrix[orderPair(y, cluster)] ) / \
             (nodes[x].clusterSize + nodes[y].clusterSize);
}

void Graph::addNode(int i, int j, double w)
{
    // merge Ci and Cj into a new cluster Cnew with |Ci| + |Cj| elements
    Node newNode(nodes[i].clusterSize + nodes[j].clusterSize);
    // add a new node labeled by cluster Cnew to T
    nodes.push_back(newNode);
    // connect node Cnew to Ci and Cj by directed edges
    addEdge(i, nodesNum, w - ages[i]);
    addEdge(j, nodesNum, w - ages[j]);
}

void Graph::UPGMA()
{
    while (matrix.size() > 1)
    {
        // find the two closest clusters Ci and Cj
        auto min_el = findMinValuePair(matrix);
        int i = min_el.first.first, j = min_el.first.second;
        double w = min_el.second / 2;
        // create new node and add to tree
        addNode(i, j, w);
        // update "ages" of nodes
        ages[nodesNum] = w;
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

void Graph::printGraph(std::ofstream& ost)
{
    for (int i = 0; i < nodes.size(); ++i)
    {
        for (auto nei: nodes[i].neighbors)
        {
            ost << i << "->" << nei.end << ":" << nei.weight << '\n';
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

    std::ofstream ost{OUTFILENAME};
    if (!ost) 
    {
        std::cout << "Cannot open output file!\n";
        exit(1);
    }
    ost << std::fixed;
    ost << std::setprecision(3);

    Graph graph(ist);
    graph.UPGMA();
    graph.printGraph(ost);
    return 0;
}