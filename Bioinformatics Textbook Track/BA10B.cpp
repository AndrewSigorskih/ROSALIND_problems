#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <limits>
#include <vector>
#include <unordered_map>

using std::string;
using std::vector;
const char* INFILENAME = "rosalind_ba10b.txt";
const char* OUTFILENAME = "out.txt";
long int STR_MAX = std::numeric_limits<std::streamsize>::max();
using Matrix = std::unordered_map<string, std::unordered_map<string, double>>;

class HMM
{
    string String;
    string Path;
    vector<string> alphabet;
    vector<string> states;
    Matrix Emission;

    vector<string> parseLine(std::istream& ist);
    void readMatrix(std::istream& ist);

public:
    HMM(std::istream& ist);
    double calcPathProb();
};

HMM::HMM(std::istream& ist)
{
    std::getline(ist, String);
    ist.ignore(STR_MAX, '\n'); // ignoring "------"
    alphabet = parseLine(ist);
    ist.ignore(STR_MAX, '\n');
    std::getline(ist, Path);
    ist.ignore(STR_MAX, '\n'); 
    states = parseLine(ist);
    ist.ignore(STR_MAX, '\n'); 
    readMatrix(ist);
}

vector<string> HMM::parseLine(std::istream& ist)
{
    string buf;
    vector<string> result;
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
    for (int i = 0; i < states.size(); ++i)
    {
        vector<string> buf = parseLine(ist);
        for (int j = 1; j <= header.size(); ++j)
        {
            Emission[buf[0]][header[j-1]] = stod(buf[j]);
        }
        buf.clear();
    }
}

double HMM::calcPathProb()
{
    double prob = 1.0;
    for (int i = 0; i < String.size(); ++i)
    {
        prob *= Emission[Path.substr(i,1)][String.substr(i,1)];
    }
    return prob;
}

int main()
{
    std::ifstream ist{INFILENAME};
    if (!ist) 
    {
        std::cerr << "Cannot open input file!\n";
        exit(1);
    }
    std::ofstream ost{OUTFILENAME};
    if (!ost) 
    {
        std::cerr << "Cannot open output file!\n";
        exit(1);
    }

    HMM model(ist);
    ost.precision(12);
    ost << model.calcPathProb() << '\n';

    return 0;
}
