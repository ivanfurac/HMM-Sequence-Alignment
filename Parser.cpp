#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>

#include <iostream>
#include <map>

#include "Values.cpp"


using namespace std;


//funkcija prima putanju do direktorija gdje se u zasebnim datotekama nalaze parovi sljedova, cita parove te vraca vektor
vector<pair<string, string>> parsePairs(string path) {

    vector<pair<string, string>> pairs;

    for (const auto &entry : filesystem::directory_iterator(path))
    {
        ifstream file(entry.path());
        string seq1;
        string seq2;
        getline(file, seq1);
        getline(file, seq2);
        file.close();
        
        std::pair <string, string> seqPair = make_pair(seq1, seq2);
        pairs.push_back(seqPair);
    }

    return pairs;
}



//ucitava vjerojatnosti prijelaza stanja iz datoteke
map <string, double> parseTransitionProbs(string path){
    map <string, double> transProbs;
    ifstream file(path);
    string line;

    for(string x : states){
        for(string y : states){
            getline(file, line);
            double prob = stod(line);
            if(prob == 0.0){
                prob = 0.00000001;
            }
            transProbs[x+y] = log(prob);
        }
    }

    file.close();
    return transProbs;
}



//ucitava pocetne vjerojatnosti stanja
map <string, double> parseStartProbs(string path){
    map <string, double> startProbs;
    ifstream file(path);
    string line;

    for(string x : states){
        getline(file, line);
        double prob = stod(line);
        if(prob == 0.0){
            prob = 0.00000001;
        }
        startProbs[x] = log(prob);
    }

    file.close();
    return startProbs;
}



//ucitava emisijske vjerojatnosti
map <string, map<string, double>> parseEmissionProbs(string path){
    ifstream file(path);
    string line;

    //za stanje match
    map <string, double> emissionProbsM;
    for(string x : symbols){
        for(string y : symbols){
            getline(file, line);
            double prob = stod(line);
            if(prob == 0.0){
                prob = 0.00000001;
            }
            emissionProbsM[x+y] = log(prob);
        }
    }

    //za stanje Ix
    map <string, double> emissionProbsX;
    for(string x : symbols){
        getline(file, line);
        double prob = stod(line);
        if(prob == 0.0){
            prob = 0.00000001;
        }
        emissionProbsX[x] = log(prob);
    }

    //za stanje Iy
    map <string, double> emissionProbsY;
    for(string x : symbols){
        getline(file, line);
        double prob = stod(line);
        if(prob == 0.0){
            prob = 0.00000001;
        }
        emissionProbsY[x] = log(prob);
    }

    file.close();

    map <string, map<string, double>> emissionProbs;
    emissionProbs["M"] = emissionProbsM;
    emissionProbs["Ix"] = emissionProbsX;
    emissionProbs["Iy"] = emissionProbsY;

    return emissionProbs;
}