#include "Parser.cpp"
#include "Viterbi.cpp"

using namespace std;

int main(void){

    vector<pair<string, string>> pairs = parsePairs("./test_pairs");
    map <string, double> Pi = parseStartProbs("./trained_params/Pi_trained.txt");
    map <string, double> A = parseTransitionProbs("./trained_params/A_trained.txt");
    map <string, map<string, double>> E = parseEmissionProbs("./trained_params/E_trained.txt");
    
    HMM *hmm = new HMM(Pi, A, E);

    int num = 0;
    for(pair<string, string> pair : pairs){
        viterbi(hmm, pair, num);
        num++;
    }   

    return 0;
}