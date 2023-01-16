#include "Baum-Welch.cpp"
#include "Parser.cpp"

using namespace std;

int main(void){

    vector<pair<string, string>> pairs = parsePairs("./train_pairs");
    map <string, double> Pi = parseStartProbs("./start_params/Pi");
    map <string, double> A = parseTransitionProbs("./start_params/A");
    map <string, map<string, double>> E = parseEmissionProbs("./start_params/E");
    
    HMM *hmm = new HMM(Pi, A, E);

    baumwelch(hmm, pairs, 10);

    return 0;
}