#include <vector>
#include <iostream>
#include <fstream> 
#include "HMM.cpp"
#include "Values.cpp"

using namespace std;

//provjera je li u stanju state (0 = M, 1 = Ix, 2 = Iy) moguce emitirati par simbola (o1, o2)
//npr. u stanju Match nije moguce emitirati "-" i "G"
//npr. u stanju Ix nije moguce emitirati "C" i "T"
bool observationPossible(int state, string o1, string o2){
    if(state == 0){
        if(o1 == "-" || o2 == "-"){
            return false;
        }
    }
    else{
        if(o1 != "-" && o2 != "-"){
            return false;
        }
        else if(state == 1 && o1 == "-"){
            return false;
        }
        else if(state == 2 && o2 == "-"){
            return false;
        }
    }
    return true;
}

//pretvara par simbola (o1, o2) u jedan string
//"A" i "A" se pretvara u "AA"
//"A" i "-" se pretvara u "A"
string convertToSymbol(string o1, string o2){
    if(o1 == "-"){
        return o2;
    }
    else if(o2 == "-"){
        return o1;
    }
    else{
        return o1 + o2;
    }
}

//forward algoritam, vraca polje duljine broj_stanja * broj_observacija
double* forward(HMM *hmm, pair <string, string> seqPair){
    string first = seqPair.first;
    string second = seqPair.second;
    string symbol1;
    string symbol2;

    int obs_num = first.length();
    if(second.length() < obs_num){
        obs_num = second.length();
    }

    double* result = (double*) malloc(3 * obs_num * sizeof(double));

    //initialization
    for(int i=0; i<3; i++){
        symbol1 = string(1, first[0]);
        symbol2 = string(1, second[0]);

        if(observationPossible(i, symbol1, symbol2)){
            result[i*obs_num] = (hmm->Pi[states[i]]) + (hmm->E[states[i]][convertToSymbol(symbol1, symbol2)]);
        }
        else{
            result[i*obs_num] = 0.0;
        }
      
    }

    //forward calculation
    for(int t=1; t<obs_num; t++){
        symbol1 = string(1, first[t]);
        symbol2 = string(1, second[t]);

        for(int i=0; i<3; i++){
            result[i*obs_num + t] = 0.0;
            for(int j=0; j<3; j++){
                if(observationPossible(i, symbol1, symbol2)){
                    result[i*obs_num + t] += result[j*obs_num + t - 1] + hmm->A[states[j]+states[i]] + hmm->E[states[i]][convertToSymbol(symbol1, symbol2)];
                }
            }
        }
    }

    return result;         
}

//backward algoritam, vraca polje duljine broj_stanja * broj_observacija
double* backward(HMM *hmm, pair <string, string> seqPair){

    string first = seqPair.first;
    string second = seqPair.second;
    string symbol1;
    string symbol2;

    int obs_num = first.length();
    if(second.length() < obs_num){
        obs_num = second.length();
    }

    double* result = (double*) malloc(3 * obs_num * sizeof(double));

    //initialization
    for(int i=0; i<3; i++){
        result[i*obs_num + obs_num - 1] = 0.0;
    }

    //backward calculation
    for(int t = obs_num-2; t>=0; t--){
        symbol1 = string(1, first[t+1]);
        symbol2 = string(1, second[t+1]);

        for(int i=0; i<3; i++){
            result[i*obs_num + t] = 0.0;
            for(int j=0; j<3; j++){
                if(observationPossible(j, symbol1, symbol2)){
                    result[i*obs_num + t] += result[j*obs_num + t + 1] + hmm->A[states[i]+states[j]] + hmm->E[states[j]][convertToSymbol(symbol1, symbol2)];
                }
            }

        }
    }

    return result;
}

//jedna iteracija Baum-Welch algoritma
//za svaki par sljedova racuna backward i forward vjerojatnosti, racuna parametre ksi i gamma i na temelju toga azurira parametre HMM-a
void run_iteration(HMM *hmm, vector<pair<string, string>> sequences){

    //inicijaliziramo nove parametre
    map <string, double> Pi;
    for(string x : states){
        Pi[x] = 0.0;
    }

    int R = 0;

    map <string, double> A_up;
    for(string x : states){
        for(string y : states){
            A_up[x+y] = 0.0;
        }
    }
    map <string, double> A_down;
    for(string x : states){
        for(string y : states){
            A_down[x+y] = 0.0;
        }
    }

    map <string, double> E_m_up;
    for(string x : symbols){
        for(string y : symbols){
            E_m_up[x+y] = 0.0;
        }
    }
    map <string, double> E_ix_up;
    for(string x : symbols){
        E_ix_up[x] = 0.0;
    }
    map <string, double> E_iy_up;
    for(string x : symbols){
        E_iy_up[x] = 0.0;
    }

    map <string, double> E_m_down;
    for(string x : symbols){
        for(string y : symbols){
            E_m_down[x+y] = 0.0;
        }
    }
    map <string, double> E_ix_down;
    for(string x : symbols){
        E_ix_down[x] = 0.0;
    }
    map <string, double> E_iy_down;
    for(string x : symbols){
        E_iy_down[x] = 0.0;
    }

    for(pair <string,string> seqPair : sequences){

        cout << "Calculating for pair " + to_string(R) << endl;

        double* forward_result = forward(hmm, seqPair);
        double* backward_result = backward(hmm, seqPair);

        string first = seqPair.first;
        string second = seqPair.second;
        string symbol1;
        string symbol2;

        int obs_num = first.length() < second.length() ? first.length() : second.length();

        double *gamma = (double*) malloc(3 * obs_num * sizeof(double));
        double *ksi = (double*) malloc(3 * 3 * (obs_num - 1) * sizeof(double));

        //gamma
        for(int t=0; t<obs_num; t++){
            //nazivnik
            double denominator = 0.0;
            for(int j=0; j<3; j++){
                denominator += forward_result[j*obs_num + t] + backward_result[j*obs_num + t];
            }

            for(int i=0; i<3; i++){
                double product = forward_result[i*obs_num + t] + backward_result[i*obs_num + t];
                if(denominator == 0.0){
                    gamma[i*obs_num + t] = 0.0;
                }
                else{
                    gamma[i*obs_num + t] = product / denominator;
                }                
            }
        }

        //ksi
        for(int t=0; t<obs_num-1; t++){
            //nazivnik
            double denominator = 0.0;
            for(int k=0; k<3; k++){
                for(int w=0; w<3; w++){
                    if(observationPossible(w, string(1, first[t+1]), string(1, second[t+1]))){
                        double factor1 = forward_result[k*obs_num + t];
                        double factor2 = hmm->A[states[k] + states[w]];
                        double factor3 = backward_result[w*obs_num + t+1];
                        double factor4 = hmm->E[states[w]][convertToSymbol(string(1, first[t+1]), string(1, second[t+1]))];
                        denominator += factor1 + factor2 + factor3 + factor4;
                    }                    
                }            
            }

            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    if(observationPossible(j, string(1, first[t+1]), string(1, second[t+1])) && denominator != 0.0){
                        double factor1 = forward_result[i*obs_num + t];
                        double factor2 = hmm->A[states[i] + states[j]];
                        double factor3 = backward_result[j*obs_num + t+1];
                        double factor4 = hmm->E[states[j]][convertToSymbol(string(1, first[t+1]), string(1, second[t+1]))];
                        ksi[t*3*3 + i*3 + j] = (factor1 + factor2 + factor3 + factor4) / denominator;
                    }
                    else{
                        ksi[t*3*3 + i*3 + j] = 0.0;
                    }
                }
            }
        }

        //update values
        R++;
        for(int i=0; i<3; i++){

            //update Pi
            Pi[states[i]] += gamma[i*obs_num];

            //update A
            for(int j=0; j<3; j++){
                for(int t=0; t<obs_num-1; t++){
                    A_up[states[i] + states[j]] += ksi[t*3*3 + i*3 + j];
                    A_down[states[i] + states[j]] += gamma[i*obs_num + t];
                }
            }
        }
        //update E
        for(int t=0; t<obs_num; t++){
            symbol1 = string(1, first[t]);
            symbol2 = string(1, second[t]);
            if(observationPossible(0, symbol1, symbol2)){
                E_m_up[symbol1 + symbol2] += gamma[t];
            }
            else if(observationPossible(1, symbol1, symbol2)){
                E_ix_up[symbol1] += gamma[obs_num + t];
            }
            else if(observationPossible(2, symbol1, symbol2)){
                E_iy_up[symbol2] += gamma[2*obs_num + t];
            }

            for(string x : symbols){
                E_ix_down[x] += gamma[obs_num + t];
                E_iy_down[x] += gamma[2*obs_num + t];
                for(string y : symbols){
                    E_m_down[x+y] += gamma[t];
                }
            }

        }

        free(forward_result);
        free(backward_result);
        free(gamma);
        free(ksi);
    }

    for(string x : states){
        Pi[x] = Pi[x] / R;
    }

    for(string x : states){
        for(string y : states){
            if(A_down[x+y] != 0.0){
                A_up[x+y] = A_up[x+y] / A_down[x+y];
            }            
        }
    }

    for(string x : symbols){
        if(E_ix_down[x] != 0.0){
            E_ix_up[x] = E_ix_up[x] / E_ix_down[x];
        }
        if(E_iy_down[x] != 0.0){
            E_iy_up[x] = E_iy_up[x] / E_iy_down[x];
        }
        for(string y : symbols){
            if(E_m_down[x+y] != 0.0){
                E_m_up[x+y] = E_m_up[x+y] / E_m_down[x+y];
            }            
        }
    }

    hmm->Pi = Pi;
    hmm->A = A_up;

    map<string, map<string, double>> E;
    E["M"] = E_m_up;
    E["Ix"] = E_ix_up;
    E["Iy"] = E_iy_up;
    hmm->E = E;
}

//glavna funkcija za pokretanje algoritma
//prolazi kroz zadan broj iteracija i rezultate sprema u posebne datoteke
void baumwelch(HMM *hmm, vector<pair<string, string>> sequences, int num_iterations){

    //zadani broj iteracija
    for(int i=0; i < num_iterations; i++){
        cout << "RUNNING ITERATION " + to_string(i) << endl;
        run_iteration(hmm, sequences);
    }

    //ispis parametara u datoteku
    ofstream pi_file("./trained_params/Pi_trained.txt");
    for(string x : states){
        pi_file << hmm->Pi[x] << endl;
    }
    pi_file.close();

    ofstream a_file("./trained_params/A_trained.txt");
    for(string x : states){
        for(string y : states){
            a_file << hmm->A[x+y] << endl;
        }
    }
    a_file.close();

    ofstream e_file("./trained_params/E_trained.txt");
    for(string x : symbols){
        for(string y : symbols){
            e_file << hmm->E["M"][x+y] << endl;
        }
    }
    for(string x : symbols){
        e_file << hmm->E["Ix"][x] << endl;
    }
    for(string x : symbols){
        e_file << hmm->E["Iy"][x] << endl;
    }
    e_file.close();
}