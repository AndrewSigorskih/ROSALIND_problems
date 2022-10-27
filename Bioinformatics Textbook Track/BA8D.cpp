#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
/* Soft k-Means Clustering Algorithm */
using std::string;
using std::vector;
const int STEPS = 100;
const char* INFILENAME = "rosalind_ba8d.txt";
const char* OUTFILENAME = "out.txt";
using matrix = vector<vector<double>>;

class softkMeans
{
    int k, m;
    double beta; // stiffness
    matrix data;
    matrix centers;
    matrix hidden;
    void E_step();
    void M_step();
    double stiffExp(int j, int i);
    double dist(vector<double> const& v, vector<double> const& w);
public:
    softkMeans(std::ifstream& ist);
    void Estimate();
    friend std::ofstream& operator << (std::ofstream& ost, softkMeans const& obj);
};

softkMeans::softkMeans(std::ifstream& ist)
{
    ist >> k >> m >> beta;
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
    // first k points = init centers
    centers.reserve(k);
    for (int i = 0; i < k; ++i)
        centers.push_back(data[i]);
    // k * n responsibility matrix
    hidden = matrix(k, vector<double>(data.size()));
}

double softkMeans::dist(vector<double> const& v, vector<double> const& w)
{   // Euclidean distance
    double res = 0.0;
    for (int i = 0; i < v.size(); ++i)
    {
        res +=pow(v[i] - w[i], 2);
    }
    return sqrt(res);
}

inline double softkMeans::stiffExp(int j, int i)
{
    return exp(-beta * dist(data[j], centers[i]));
}

void softkMeans::E_step()
{
    for (int j = 0; j < data.size(); ++j)
    {
        double cumsum = 0.0;
        vector<double> tmp(k);
        for (int i = 0; i < k; ++i)
        {
            tmp[i] = stiffExp(j, i);
            cumsum += tmp[i];
        }
        for (int i = 0; i < k; ++i)
            hidden[i][j] = tmp[i] / cumsum;
    }
}

void softkMeans::M_step()
{
    for (int i = 0; i < k; ++i)
    {
        vector<double> new_point(m);
        for (int pos = 0; pos < m; ++pos)
        {
            double divisor = 0.0, denominator = 0.0;
            for (int j = 0; j < data.size(); ++j)
            {
                divisor += hidden[i][j] * data[j][pos];
                denominator += hidden[i][j];
            }
            centers[i][pos] = divisor / denominator;
        }
    }
}

void softkMeans::Estimate()
{
    for (int i = 0; i < STEPS; ++i)
    {
        // E-step
        // assign each data point a “responsibility” value for each cluster
        E_step();
        // M-step
        // compute new centers
        M_step();
    }
}

std::ofstream& operator << (std::ofstream& ost, softkMeans const& obj)
{
    for (int i = 0; i < obj.k; ++i)
    {
        for (int j = 0; j < obj.centers[i].size(); ++j)
        {
            ost << obj.centers[i][j] << " ";
        }
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

    softkMeans solver(ist);
    solver.Estimate();
    ost << solver;
    return 0;
}