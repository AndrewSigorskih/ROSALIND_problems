#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <limits>
#include <vector>
#include <unordered_map>

using std::string;
using std::vector;
const char* INFILENAME = "rosalind_ba10c.txt";
const char* OUTFILENAME = "out.txt";
long int STR_MAX = std::numeric_limits<std::streamsize>::max();
using Matrix =  std::unordered_map<string, double>;

class HMM
{
    string Path;
    string String;
    vector<string> alphabet;
    vector<string> states;
    Matrix Transition;
    Matrix Emission;

    vector<string> parseLine(std::istream& ist);
    void readMatrix(std::istream& ist, Matrix& mat);

public:
    HMM(std::istream& ist);
    void Viterbi();
    string getPath();
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

// https://en.wikipedia.org/wiki/Viterbi_algorithm
void HMM::Viterbi()
{
    if (Path.size() > 0)
        Path.clear();
    // TState stores the probability of the most likely path so far
    // TIndex stores previous state of the most likely path so far
    vector<vector<int>> TIndex(String.size(), vector<int>(states.size()));
    vector<vector<double>> TState(String.size(), vector<double>(states.size()));
    // init
    for (int j = 0; j < states.size(); ++j)
    {
        TState[0][j] = (1.0 / states.size()) * Emission[states[j] + String[0]]; 
        TIndex[0][j] = 0;
    }

    // fill tables
    for (int i = 1; i < String.size(); ++i)
    {
        for (int j = 0; j < states.size(); ++j)
        {
            double em = Emission[states[j] + String[i]];
            for (int k = 0; k < states.size(); ++k)
            {
                double prob = TState[i-1][k] * em * Transition[states[k] + states[j]];
                if (prob > TState[i][j])
                {
                    TState[i][j] = prob;
                    TIndex[i][j] = k;
                }
            }

        }
    }
    // decode path
    int argmax = 0;
    for (int k = 1; k < states.size(); ++k)
    {
        if (TState[String.size()-1][k] > TState[String.size()-1][argmax])
            argmax = k;
    }
    vector<int> tmp(String.size());
    tmp[String.size()-1] = argmax;
    for (int i = String.size()-1; i >= 1; --i)
    {
        tmp[i-1] = TIndex[i][tmp[i]];
    }
    // create string
    for (int i = 0; i < tmp.size(); ++i)
    {
        Path += states[tmp[i]];
    }
}  

string HMM::getPath()
{
    return Path;
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
    model.Viterbi();
    ost << model.getPath() << '\n';

    return 0;
}
