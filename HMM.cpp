#include <map>
#include <string>

using namespace std;

class HMM {
    public:
        //pocetne vjerojatnosti stanja
        map <string, double> Pi;

        //vjerojatnosti prijelaza
        map <string, double> A;

        //emisijske vjerojatnosti
        map <string, map <string, double>> E;

        //konstruktor
        HMM(map <string, double> Pi, map <string, double> A, map <string, map <string, double>> E){
            this->Pi = Pi;
            this->A = A;
            this->E = E;
        }
};