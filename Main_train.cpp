#include "Baum-Welch.cpp"
#include "Parser.cpp"

using namespace std;

int main(void){

    vector<pair<string, string>> pairs = parsePairs("./train_pairs");
    map <string, double> Pi = parseStartProbs("./start_params/Pi");
    map <string, double> A = parseTransitionProbs("./start_params/A");
    map <string, map<string, double>> E = parseEmissionProbs("./start_params/E");
    
    HMM *hmm = new HMM(Pi, A, E);

    int num_iterations = 5;

    baumwelch(hmm, pairs, num_iterations);

    return 0;
}