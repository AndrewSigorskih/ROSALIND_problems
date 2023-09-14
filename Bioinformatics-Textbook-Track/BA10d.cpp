#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <limits>
#include <vector>
#include <unordered_map>

using std::string;
using std::vector;
const char* INFILENAME = "rosalind_ba10d.txt";
const char* OUTFILENAME = "out.txt";
long int STR_MAX = std::numeric_limits<std::streamsize>::max();
using Matrix =  std::unordered_map<string, double>;

class HMM
{
    string String;
    vector<string> alphabet;
    vector<string> states;
    Matrix Transition;
    Matrix Emission;

    vector<string> parseLine(std::istream& ist);
    void readMatrix(std::istream& ist, Matrix& mat);

public:
    HMM(std::istream& ist);
    double OutcomeLikelihood();
};

HMM::HMM(std::istream& ist)
{
    std::getline(ist, String);
    ist.ignore(STR_MAX, '\n'); // ignoring "------"
    alphabet = parseLine(ist);
    ist.ignore(STR_MAX, '\n');
    states = parseLine(ist);
    ist.ignore(STR_MAX, '\n'); 
    readMatrix(ist, Transition);
    ist.ignore(STR_MAX, '\n'); 
    readMatrix(ist, Emission);
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

void HMM::readMatrix(std::istream& ist, Matrix& mat)
{
    vector<string> header = parseLine(ist);
    for (int i = 0; i < states.size(); ++i)
    {
        vector<string> buf = parseLine(ist);
        for (int j = 1; j <= header.size(); ++j)
        {
            mat[buf[0] + header[j-1]] = stod(buf[j]);
        }
        buf.clear();
    }
}

double HMM::OutcomeLikelihood()
{
    //vector<vector<double>> TState(states.size(), vector<double>(String.size()));
    vector<vector<double>> TState(String.size(), vector<double>(states.size()));
    // init
    for (int j = 0; j < states.size(); ++j)
    {
        TState[0][j] = (1.0 / states.size()) * Emission[states[j] + String[0]]; 
    }
    // fill tables
    for (int i = 1; i < String.size(); ++i)
    {
        for (int j = 0; j < states.size(); ++j)
        {
            for (int k = 0; k < states.size(); ++k)
            {
                TState[i][j] += TState[i-1][k] *\
                                Transition[states[k] + states[j]] *\
                                Emission[states[j] + String[i]];
            }
        }
    }
    double result = 0.0;
    for (int j = 0; j < states.size(); ++j)
        result += TState[String.size()-1][j];
    return result;
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
    ost << model.OutcomeLikelihood() << '\n';

    return 0;
}