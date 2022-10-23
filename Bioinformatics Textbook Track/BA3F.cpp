#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <sstream>

using std::string;
using std::vector;
using std::stack;
const char* INFILENAME = "rosalind_ba3f.txt";
const char* OUTFILENAME = "out.txt";
using adjlist = vector<vector<int>>;

// Hierholzer’s Algorithm for directed graphs
// source: https://www.geeksforgeeks.org/hierholzers-algorithm-directed-graph/

class Graph
{
    adjlist nodes;
public:
    Graph(std::ifstream&);
    void findEulerCircuit(std::ofstream&);
};

Graph::Graph(std::ifstream& ist)
{   // graph of unknown size 
    string tmp;
    while(getline(ist, tmp))
    {
        int node, nei;
        nodes.push_back(vector<int>());
        std::stringstream ss(tmp);
        ss >> node; // in "real" cases nodes do NOT come in order!!!
        if (node > (nodes.size()-1))
        { // create space for ALL missing nodes
            nodes.reserve(node);
            for (int i = nodes.size(); i <= node; ++i)
                nodes.push_back(vector<int>());
        }
        ss.ignore(4, '>');
        ss >> nei;
        nodes[node].push_back(nei);
        while (ss.peek() == ',')
        { // if more than one neighbor
            ss.ignore(1);
            ss >> nei;
            nodes[node].push_back(nei);
        }
    }
}

void Graph::findEulerCircuit(std::ofstream& ost)
{
    vector<int> result;
    result.reserve(nodes.size());
    stack<int> curr_path;

    int curr_v = 0;
    curr_path.push(curr_v);

    while (!curr_path.empty())
    {   
        curr_v = curr_path.top();
        if (nodes[curr_v].size()) // curr node still has edges
        {
            int next_v = nodes[curr_v].back();
            nodes[curr_v].pop_back();
            curr_path.push(next_v);
        } else {
            result.push_back(curr_path.top());
            curr_path.pop();
        }
    }
    for (int i = result.size()-1; i >= 0; --i)
    {
        ost << result[i];
        if (i)
            ost << "->";
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

    Graph graph(ist);
    graph.findEulerCircuit(ost);
    return 0;
}