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
using Matrix = std::unordered_map<string, std::unordered_map<string, double>>;

class HMM
{
    string String;
    string Path;
    vector<string> alphabet;
    vector<string> states;
    Matrix Transition;
    Matrix Emission;

    vector<string> parseLine(std::istream& ist);
    void readMatrix(std::istream& ist, const char mat);

public:
    HMM(std::istream& ist);
    void Viterbi();
    string getPath();
};

HMM::HMM(std::istream& ist)
{
    Path = "";
    std::getline(ist, String);
    ist.ignore(STR_MAX, '\n'); // ignoring "------"
    alphabet = parseLine(ist);
    ist.ignore(STR_MAX, '\n');
    states = parseLine(ist);
    ist.ignore(STR_MAX, '\n'); 
    readMatrix(ist, 'T');
    ist.ignore(STR_MAX, '\n'); 
    readMatrix(ist, 'E');
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

void HMM::readMatrix(std::istream& ist, const char mat)
{
    Matrix* m;
    switch (mat)
    {
        case 'T':
            m = &Transition;
            break;
        case 'E':
            m = &Emission;
            break;
        default:
            std::cerr << "HMM::readMatrix got incorrect matrix type arg: ";
            std::cerr << mat << '\n';
            exit(1);
    }

    vector<string> header = parseLine(ist);
    for (int i = 0; i < states.size(); ++i)
    {
        vector<string> buf = parseLine(ist);
        for (int j = 1; j <= header.size(); ++j)
        {
            (*m)[buf[0]][header[j-1]] = stod(buf[j]);
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
    vector<vector<int>> TIndex(states.size(), vector<int>(String.size()));
    vector<vector<double>> TState(states.size(), vector<double>(String.size()));
    // init
    for (int j = 0; j < states.size(); ++j)
    {
        TState[j][0] = (1.0 / states.size()) * Emission[states[j]][String.substr(0,1)]; 
        TIndex[j][0] = 0;
    }
    // fill tables
    for (int i = 1; i < String.size(); ++i)
    {
        for (int j = 0; j < states.size(); ++j)
        {
            // argmax
            int argmax = 0;
            double prob = 0.0;
            for (int k = 0; k < states.size(); ++k)
            {
                double curprob = TState[k][i-1] *\
                                 Transition[states[k]][states[j]] *\
                                 Emission[states[j]][String.substr(i,1)];
                if (curprob > prob)
                {
                    argmax = k;
                    prob = curprob;
                }
            }
            TIndex[j][i] = argmax;
            TState[j][i] = TState[TIndex[j][i]][i-1] *\
                            Transition[states[argmax]][states[j]] *\
                            Emission[states[j]][String.substr(i,1)];
        }
    }
    // decode path
    int argmax = 0;
    for (int k = 1; k < states.size(); ++k)
    {
        if (TState[k][String.size()-1] > TState[argmax][String.size()-1])
            argmax = k;
    }

    vector<int> tmp(String.size());
    tmp[String.size()-1] = argmax;
    for (int i = String.size()-1; i >= 1; --i)
    {
        tmp[i-1] = TIndex[tmp[i]][i];
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