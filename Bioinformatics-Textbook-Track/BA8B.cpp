#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>
/* Squared Error Distortion */
using std::string;
using std::vector;
const char* INFILENAME = "rosalind_ba8b.txt";
const char* OUTFILENAME = "out.txt";

class kCenters
{
    int k;
    vector<vector<double>> centers;
    vector<vector<double>> data;
    double minDistToCenters(int);
    double dist(vector<double> const& v, vector<double> const& w);
public:
    kCenters(std::ifstream&);
    double SquaredErrorDistortion();
};

kCenters::kCenters(std::ifstream& ist)
{
    int m;
    ist >> k >> m;
    for (int i = 0; i < k; ++i)
    {
        vector<double> point(m);
        for (int i = 0; i < m; ++i)
            ist >> point[i];
        centers.push_back(point);
    }
    ist.ignore(INT32_MAX, '\n');
    ist.ignore(INT32_MAX, '\n');
    string tmp;
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
    for (int i = 0; i < centers.size(); ++i)
    {
        double cur_dist = dist(data[idx], centers[i]);
        if (cur_dist < distance)
            distance = cur_dist;
    }
    return distance;
}

double kCenters::SquaredErrorDistortion()
{
    double distortion = 0;
    for (int i = 0; i < data.size(); ++i)
        distortion += pow(minDistToCenters(i), 2);
    return distortion * (1.0 / data.size());
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

    kCenters solver(ist);
    ost << solver.SquaredErrorDistortion() << '\n';
    return 0;
}