#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>
#include <algorithm>
#include <unordered_set>
/* FarthestFirstTraversal */
using std::string;
using std::vector;
using set = std::unordered_set<int>;
const char* INFILENAME = "rosalind_ba8a.txt";
const char* OUTFILENAME = "out.txt";

class kCenters
{
    int k;
    set centers; // store centers' indexes, not copy of data
    vector<vector<double>> data;
    double dist(vector<double> const&, vector<double> const&);
    double minDistToCenters(int idx);
public:
    kCenters(std::ifstream&);
    void FarthestFirstTraversal();
    void printCenters(std::ofstream&);
};

kCenters::kCenters(std::ifstream& ist)
{
    int m;
    ist >> k >> m;
    string tmp;
    ist.ignore(INT32_MAX, '\n');
    while (getline(ist, tmp))
    {
        std::stringstream ss(tmp);
        vector<double> point(m);
        for (int i = 0; i < m; ++i)
            ss >> point[i];
        data.push_back(point);
    }
}

double kCenters::dist(vector<double> const& v, vector<double> const& w)
{   // Euclidean distance
    double res = 0.0;
    for (int i = 0; i < v.size(); ++i)
    {
        res +=pow(v[i] - w[i], 2);
    }
    return sqrt(res);
}

double kCenters::minDistToCenters(int idx)
{
    double distance = std::numeric_limits<double>::max();
    for (int i: centers)
    {
        double cur_dist = dist(data[idx], data[i]);
        if (cur_dist < distance)
            distance = cur_dist;
    }
    return distance;
}

void kCenters::FarthestFirstTraversal()
{
    centers.insert(0);
    while (centers.size() < k)
    {
        int new_center = -1;
        double max_dist = 0.0;
        for (int i = 0; i < data.size(); ++i)
        {
            if (centers.find(i) != centers.end())
                continue; // skip comparing centers to themselves
            double cur_dist = minDistToCenters(i);
            if (cur_dist > max_dist)
            {
                max_dist = cur_dist;
                new_center = i;
            }
        }
        centers.insert(new_center);
    }
}

void kCenters::printCenters(std::ofstream& ost)
{
    // faster than using tree-based set in main loop
    vector<int> result;
    result.reserve(centers.size());
    for (auto i: centers)
        result.push_back(i);
    std::sort(result.begin(), result.end());

    for (auto i : result)
    {
        for (int j = 0; j < data[i].size(); ++j)
            ost << data[i][j] << " ";
        ost << '\n';
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
    ost << std::setprecision(1);

    kCenters solver(ist);
    solver.FarthestFirstTraversal();
    solver.printCenters(ost);
    
    return 0;
}