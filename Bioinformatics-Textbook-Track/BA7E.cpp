#include "BA7.h"
const char* INFILENAME = "rosalind_ba7e.txt";
const char* OUTFILENAME = "out.txt";

struct Node
{
    vector<Edge> neighbors;
    int idx;
    Node(int num) : idx(num) {};
};

using adjlist = vector<Node>;

class Graph
{
    int nodesNum;
    adjlist nodes;
    vector<int> labels;
    distMatrix D_matrix;
    distMatrix Q_matrix;
    void compute_Q(int const);
    double RowSum(int const ind);
    void addEdge(int start, int end, double weight);
    void addNode(int i, int j, double w_i, double w_j);
public:
    Graph(std::ifstream&);
    void NJ();
    friend std::ofstream& operator << (std::ofstream& ost, Graph const& obj);
};

Graph::Graph(std::ifstream& ist)
{
    ist >> nodesNum;
    nodes.reserve(2 * nodesNum - 2);
    labels.reserve(nodesNum);
    for (int i = 0; i < nodesNum; ++i)
    {
        nodes.push_back(Node(i));
        labels.push_back(i);
        for (int j = 0; j < nodesNum; ++j)
        {
            double val;
            ist >> val;
            if (i >= j)
                continue;
            D_matrix[orderPair(i, j)] = val;
        }
    }
}

void Graph::addEdge(int start, int end, double weight)
{
    nodes[start].neighbors.push_back(Edge{end, weight});
    nodes[end].neighbors.push_back(Edge{start, weight});
}

void Graph::addNode(int i, int j, double w_i, double w_j)
{
    // create new node
    Node newNode(nodesNum);
    // add a new node labeled by Cnew to T
    nodes.push_back(newNode);
    // connect node Cnew to Ci and Cj by directed edges
    addEdge(i, nodesNum, w_i);
    addEdge(j, nodesNum, w_j);
}

double Graph::RowSum(int const ind)
{
    double res = 0.0;
    for (auto k: labels)
    {
        if (k != ind)
            res += D_matrix[orderPair(ind, k)];
    }
    return res;
}

void Graph::compute_Q(int const n)
{
    Q_matrix.clear();
    for (auto i : labels)
    {
        for (auto j: labels)
        {
            if (i >= j)
                continue;
            Q_matrix[orderPair(i, j)] = \
                            (n-2) * D_matrix[orderPair(i, j)] - RowSum(i) - RowSum(j);
        }
    }
}

void Graph::NJ()
{
    int n = labels.size();
    while (n > 2)
    {   // recalc Q_matrix
        compute_Q(n);
        // find min_element in Q_matrix
        auto min_el = findMinValuePair(Q_matrix);
        int i = min_el.first.first, j = min_el.first.second;
        double w = D_matrix[orderPair(j, i)];
        // get distances from i and j to new node
        double w_i = (D_matrix[orderPair(i, j)] +\
                (1.0 / (n-2)) * (RowSum(i) - RowSum(j))) / 2;
        double w_j = w - w_i;
        // add new node to tree
        addNode(i, j, w_i, w_j);
        // recalc distance to new node
        for (auto k: labels)
        {
            if ((k != i) && (k != j))
                D_matrix[orderPair(k, nodesNum)] = \
                    0.5 * (D_matrix[orderPair(i, k)] + D_matrix[orderPair(j, k)] - w);
            
        }
        // remove nodes i and j from matrix and labels
        for (auto k: labels)
        {
            if ((k != i) && (k != j))
            {
                D_matrix.erase(orderPair(i, k));
                D_matrix.erase(orderPair(j, k));
            }
        }
        D_matrix.erase(orderPair(i, j));
        labels.erase(std::remove(labels.begin(), labels.end(), i), labels.end());
        labels.erase(std::remove(labels.begin(), labels.end(), j), labels.end());
        // add new node to labels and increment indices
        labels.push_back(nodesNum);
        ++nodesNum;
        n = labels.size();
    }
    // add root
    int i = labels[0], j = labels[1];
    double w = D_matrix[orderPair(i, j)];
    addEdge(i, j, w);
}

std::ofstream& operator << (std::ofstream& ost, Graph const& obj)
{
    for (int i = 0; i < obj.nodes.size(); ++i)
    {
        for (auto nei: obj.nodes[i].neighbors)
        {
            ost << i << "->" << nei.end << ":" << nei.weight << '\n';
        }
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
    ost << std::setprecision(2);

    Graph graph(ist);
    graph.NJ();
    ost << graph;
    return 0;
}