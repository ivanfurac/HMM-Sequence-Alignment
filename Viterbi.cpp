#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include "HMM.cpp"

using namespace std;

const double INF = -1e9;

string reverseString(string first){
    string second = "";
    for(int i = first.length()-1; i >= 0; i--){
        second += first[i];
    }
    return second;
}

void viterbi(HMM *hmm, pair<string, string> observation, int pair_num){

    cout << "Aligning pair " + to_string(pair_num) << endl;

    double epsilon = (hmm->A["IxIx"] + hmm->A["IyIy"]) / 2;
    double delta = (hmm->A["MIx"] + hmm->A["MIy"]) / 2;

    //slijedovi
    string first = observation.first;
    string second = observation.second;
    int n = first.length();
    int m = second.length();
    string symbol1;
    string symbol2;

    //viterbijeve matrice
    double* Vm_old = (double *) malloc((m+1)*sizeof(double));
    double* Vx_old = (double *) malloc((m+1)*sizeof(double));
    double* Vy_old = (double *) malloc((m+1)*sizeof(double));

    double* Vm_new = (double *) malloc((m+1)*sizeof(double));
    double* Vx_new = (double *) malloc((m+1)*sizeof(double));
    double* Vy_new = (double *) malloc((m+1)*sizeof(double));

    //matrice pokazivaca na najbolji put
    int *path_m = (int *) malloc((n+1)*(m+1)*sizeof(int));
    int *path_x = (int *) malloc((n+1)*(m+1)*sizeof(int));
    int *path_y = (int *) malloc((n+1)*(m+1)*sizeof(int));

    //inicijalizacija vrijednosti
    for(int j=0; j<=m; j++){
        Vm_old[j] = INF;
        Vx_old[j] = INF;
        Vy_old[j] = INF;
    }
    Vm_old[0] = 0.0;
    Vm_new[0] = INF;
    Vx_new[0] = INF;
    Vy_new[0] = INF;

    //for petlje
    for(int i=1; i<=n; i++){
        for(int j=1; j<=m; j++){
            symbol1 = string(1, first[i-1]);
            symbol2 = string(1, second[i-1]);

            //stanje M
            double v1 = Vm_old[j-1];
            double v2 = Vx_old[j-1];
            double v3 = Vy_old[j-1];

            if(v1 >= v2 && v1 >= v3){
                Vm_new[j] = hmm->E["M"][symbol1 + symbol2] + v1;
                path_m[i*(m+1) + j] = 0;
            }
            else if(v2 >= v1 && v2 >= v3){
                Vm_new[j] = hmm->E["M"][symbol1 + symbol2] + v2;
                path_m[i*(m+1) + j] = 1;
            }
            else if(v3 >= v1 && v3 >= v2){
                Vm_new[j] = hmm->E["M"][symbol1 + symbol2] + v3;
                path_m[i*(m+1) + j] = 2;
            }

            //stanje Ix
            double v4 = Vm_old[j] + delta;
            double v5 = Vx_old[j] + epsilon;

            if(v4 >= v5){
                Vx_new[j] = hmm->E["Ix"][symbol1] + v4;
                path_x[i*(m+1) + j] = 0;
            }
            else if(v5 >= v4){
                Vx_new[j] = hmm->E["Ix"][symbol1] + v5;
                path_x[i*(m+1) + j] = 1;
            }

            //stanje Iy
            double v6 = Vm_new[j - 1] + delta;
            double v7 = Vy_new[j - 1] + epsilon;

            if(v6 >= v7){
                Vy_new[j] = hmm->E["Iy"][symbol2] + v6;
                path_y[i*(m+1) + j] = 0;
            }
            else if(v7 >= v6){
                Vy_new[j] = hmm->E["Iy"][symbol2] + v7;
                path_y[i*(m+1) + j] = 2;
            }
        }

        memcpy(Vm_old, Vm_new, (m+1)*sizeof(double));
        memcpy(Vx_old, Vx_new, (m+1)*sizeof(double));
        memcpy(Vy_old, Vy_new, (m+1)*sizeof(double));
    }

    //generiranje liste stanja
    int *path;
    double max_m = Vm_old[m];
    double max_x = Vx_old[m];
    double max_y = Vy_old[m];

    string result = "";

    if(max_m >= max_x && max_m >= max_y){
        path = path_m;
        result += "M";
    }
    else if(max_x >= max_m && max_x >= max_y){
        path = path_x;
        result += "X";
    }
    else{
        path = path_y;
        result += "Y";
    }

    int i = n;
    int j = m;
    int state;
    while(i != 0 && j != 0){
        state = path[i*(m+1) + j];
        if(state == 0){
            result += "M";
            i = i-1;
            j = j-1;
            path = path_m;

        }
        else if(state == 1){
            result += "X";
            i = i-1;
            path = path_x;
        }
        else{
            result += "Y";
            j = j-1;
            path = path_y;
        }
    }

    result = reverseString(result);

    //generiranje poravnanja
    string aligned1 = "";
    string aligned2 = "";
    i = 0;
    j = 0;
    for(char state : result){
        if(state == 'M'){
            aligned1 += first[i];
            aligned2 += second[j];
            i++;
            j++;
        }
        else if(state == 'X'){
            aligned1 += first[i];
            aligned2 += "-";
            i++;
        }
        else{
            aligned1 += "-";
            aligned2 += second[j];
            j++;
        }

    }

    string filename = "test_no_" + to_string(pair_num) + ".txt";

    ofstream file("./alignments/" + filename);
    file << aligned1 << endl;
    file << aligned2 << endl;
    file.close();
}