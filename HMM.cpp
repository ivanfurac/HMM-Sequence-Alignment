#include <iostream>
using namespace std;

class HMM {
    private:
        //skrivena stanja
        int n;
        char* states;
        //opažanja
        int m;
        char* symbols;
        //matrice vjerojatnosti
        float* A;
        float* E;
        float* Pi;
    public:
        HMM(int n, char* states, int m, char* symbols, float* A, float* E, float* Pi){
            this->n = n;
            this->states = states;
            this->m = m;
            this->symbols = symbols;
            this->A = A;
            this->E = E;
            this->Pi = Pi;
        }

        void print() {
            //Ispis svih mogućih stanja
            cout << "States:" << endl;;
            for(int i = 0; i < this->n; i++){
                cout << this->states[i] << " ";
            }
            cout << endl << endl;
            
            //Ispis svih mogućih simbola
            cout << "Symbols: " << endl;
            for(int i = 0; i < this->m; i++){
                cout << this->symbols[i] << " ";
            }
            cout << endl << endl;;

            //Ispis početnih vjerojatnosti
            cout << "startProbs: " << endl;
            for(int i=0; i<this->n; i++){
                cout << this->Pi[i] << " ";
            }
            cout << endl << endl;

            //Ispis vjerojatnosti prijelaza
            cout << "transProbs: " << endl;
            for(int i=0; i<this->n; i++){
                for(int j=0; j<this->n; j++){
                    cout << this->A[i*(this->n) + j] << " ";
                }
                cout << endl;
            }
            cout << endl;

            //Ispis emisijskih vjerojatnosti
            cout << "emissionProbs: " << endl;
            for(int i=0; i<this->n; i++){
                for(int j=0; j<this->m; j++){
                    cout << this->E[i*(this->m) + j] << " ";
                }
                cout << endl;
            }
        }

        void forward(char* observations, int obs_num){
            int* indexes = mapping(observations, obs_num);
            float* result = (float*) malloc(n * obs_num * sizeof(float));
            //initialization
            for(int i=0; i<this->n; i++){
                result[i*obs_num] = (this->Pi[i])*(this->E[i*m + indexes[0]]);
            }

            for(int t=1; t<obs_num; t++){
                for(int i=0; i<this->n; i++){
                    result[i*obs_num + t] = 0;
                    for(int j=0; j<this->n; j++){
                        result[i*obs_num + t] += result[j*obs_num + t - 1] * (this->A[j*this->n + i]);
                    }
                    result[i*obs_num + t] *= this->E[i*this->m + indexes[t]];
                }
            }

            cout << "forwardProbs: " << endl;
            for(int i=0; i<this->n; i++){
                for(int j=0; j<obs_num; j++){
                    cout << result[i*(obs_num) + j] << " ";
                }
                cout << endl;
            }
            
        }

        void backward(char* observations, int obs_num){
            int* indexes = mapping(observations, obs_num);
            float* result = (float*) malloc(n * obs_num * sizeof(float));
            //initialization
            for(int i=0; i<this->n; i++){
                result[i*obs_num + obs_num - 1] = 1;
            }

            for(int t=obs_num-2; t>=0; t--){
                for(int i=0; i<this->n; i++){
                    result[i*obs_num + t] = 0;
                    for(int j=0; j<this->n; j++){
                        result[i*obs_num + t] += result[j*obs_num + t + 1] * (this->A[i*this->n + j]) * (this->E[j*this->m + indexes[t+1]]);
                    }
                }
            }

            cout << "backwardProbs: " << endl;
            for(int i=0; i<this->n; i++){
                for(int j=0; j<obs_num; j++){
                    cout << result[i*(obs_num) + j] << " ";
                }
                cout << endl;
            }
        }

        int* mapping(char* observations, int obs_num){
            int* indexes = (int*) malloc(obs_num*sizeof(int));
            for(int i=0; i<obs_num; i++){
                for(int j=0; j<this->m; j++){
                    if(observations[i] == this->symbols[j]){
                        indexes[i] = j;
                        break;
                    }
                }
            }
            return indexes;
        }

};

int main(void){
    int n = 2;
    int m = 6;
    char states[] = {'f', 'b'};
    char symbols[] = {'1', '2', '3', '4', '5', '6'};
    float A[] = {0.8, 0.2, 0.3, 0.7};
    float E[] = {(float) 1/6, (float) 1/6, (float) 1/6, (float) 1/6, (float) 1/6, (float) 1/6, (float) 1/10, (float) 1/10, (float) 1/10, (float) 1/10, (float) 1/10, (float) 1/2};
    float Pi[] = {0.5, 0.5};
    HMM *hmm = new HMM(n, states, m, symbols, A, E, Pi);

    hmm->print();

    char obs[] = {'5', '6', '2', '3', '6', '6'};
    hmm->backward(obs, 6);
}
