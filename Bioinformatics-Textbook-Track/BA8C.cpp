#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
/* Lloyd Algorithm for k-Means Clustering */
using std::string;
using std::vector;
const char* INFILENAME = "rosalind_ba8c.txt";
const char* OUTFILENAME = "out.txt";

class kMeans
{
    int k;
    vector<int> clusterIds; // which center is closest to point
    vector<vector<double>> centers; // cluster centers positions
    vector<vector<double>> data; // our data
    vector<int> pointsCount; // how many points in a cluster
    int closestCenter(int const); 
    double dist(vector<double> const& v, vector<double> const& w);
public:
    kMeans(std::ifstream&);
    void Estimate();
    friend std::ofstream& operator << (std::ofstream& ost, kMeans const& obj);
};

std::ofstream& operator << (std::ofstream& ost, kMeans const& obj)
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

kMeans::kMeans(std::ifstream& ist)
{
    int m; 
    ist >> k >> m;
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
    centers.reserve(k);
    for (int i = 0; i < k; ++i)
    {
        vector<double> point(data[i]);
        centers.push_back(point);
    }
    clusterIds = vector<int>(data.size());
    pointsCount = vector<int>(k);
}

double kMeans::dist(vector<double> const& v, vector<double> const& w)
{   // Euclidean distance
    double res = 0.0;
    for (int i = 0; i < v.size(); ++i)
    {
        res +=pow(v[i] - w[i], 2);
    }
    return sqrt(res);
}

int kMeans::closestCenter(int const idx)
{
    int res = 0;
    double distance = dist(data[idx], centers[res]);
    for (int i = 1; i < centers.size(); ++i)
    {
        double cur_dist = dist(data[idx], centers[i]);
        if (cur_dist < distance)
        {
            distance = cur_dist;
            res = i;
        }
    }
    return res;
}

void kMeans::Estimate()
{
    while (true)
    {
        // Centers to Clusters
        // assign each data point to its nearest center
        std::fill(pointsCount.begin(), pointsCount.end(), 0);
        for (int i = 0; i < data.size(); ++i)
        {
            clusterIds[i] = closestCenter(i);
            pointsCount[clusterIds[i]] += 1;
        }
        // Clusters to Centers
        // for each cluster, recompute its center
        vector<vector<double>> new_centers(k, vector<double>(data[0].size()));
        for (int i = 0; i < data.size(); ++i)
        {
            int clust_id = clusterIds[i];
            for (int j = 0; j < data[i].size(); ++j)
            {
                new_centers[clust_id][j] += data[i][j] / pointsCount[clust_id];
            }
        }
        // check if centers changed
        bool not_changed = true;
        for (int j = 0; j < k; ++j)
        {
            if (centers[j] != new_centers[j])
            {
                not_changed = false;
                break;
            }
        }
        // if all not changed -> stop; else store new centers and continue
        if (not_changed)
        {
            break;
        } else {
            centers = new_centers;
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

    kMeans solver(ist);
    solver.Estimate();
    ost << solver;
    return 0;
}