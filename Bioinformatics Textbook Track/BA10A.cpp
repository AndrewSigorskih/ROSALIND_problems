#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <limits>
#include <vector>
#include <unordered_map>

using std::string;
using std::vector;
long int STR_MAX = std::numeric_limits<std::streamsize>::max();
using Matrix = std::unordered_map<string, std::unordered_map<string, double>>;

class HMM
{
    string Path;
    vector<string> alphabet;
    Matrix Transition;

    vector<string> parseLine(std::istream& ist);
    void readMatrix(std::istream& ist);

public:
    HMM(std::istream& ist);
    double calcPathProb();
};

HMM::HMM(std::istream& ist)
{
    std::getline(ist, Path);
    ist.ignore(STR_MAX, '\n'); // ignoring "------"
    alphabet = parseLine(ist);
    ist.ignore(STR_MAX, '\n'); // ignoring "------"
    readMatrix(ist);
}

vector<string> HMM::parseLine(std::istream& ist)
{
    string buf;
    vector<string> result;
    //ist >> buf;
    std::getline(ist, buf);
    std::istringstream iss(buf);
    std::copy(std::istream_iterator<string>(iss), 
              std::istream_iterator<string>(),
              std::back_inserter(result));
	return result;
}

void HMM::readMatrix(std::istream& ist)
{
    vector<string> header = parseLine(ist);
    for (int i = 0; i < header.size(); ++i)
    {
        vector<string> buf = parseLine(ist);
        for (int j = 1; j <= header.size(); ++j)
        {
            Transition[buf[0]][header[j-1]] = stod(buf[j]);
        }
        buf.clear();
    }
}

double HMM::calcPathProb()
{
    double prob = 1.0 / alphabet.size();
    for (int i = 1; i < Path.size(); ++i)
    {
        prob *= Transition[Path.substr(i-1,1)][Path.substr(i,1)];
    }
    return prob;
}

int main()
{
    std::ifstream ist{"rosalind_ba10a.txt"};
    if (!ist) 
    {
        std::cout << "Cannot open input file!\n";
        exit(1);
    }
    std::ofstream ost{"out.txt"};
    if (!ost) 
    {
        std::cout << "Cannot open output file!\n";
        exit(1);
    }
    HMM model(ist);
    ost.precision(12);
    ost << model.calcPathProb() << '\n';

    return 0;
}
