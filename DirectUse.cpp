//
// Created by const on 8/4/15.
//
#include "KVOR.h"
#include "TOV.h"
#include "KVDriver.h"
#include <iostream>
#include <sstream>
using namespace std;

KVDriver * parseEos(string eosfile){
    ifstream f(eosfile);
    int num_lines = 0;
    string line;
    while (getline(f, line))
        if (line[0] != '#')
            num_lines++;
    f.clear();
    f.seekg(0, ios::beg);

    double * N = new double[num_lines];
    double * E = new double[num_lines];
    double * P = new double[num_lines];
    int i = 0;
    while (getline(f, line)){
        if (line[0] != '#') {
//            cout << line << endl;
            stringstream line_ss(line);
            line_ss >> N[i] >> E[i] >> P[i];
            i++;
        }
    }
    for (int i = 1; i < num_lines; i++){
        if (E[i] < E[i-1] || N[i] < N[i-1] || (P[i] < P[i-1])){
            cout << i << " " << N[i] << " " << E[i] << " " << P[i] << endl;
            cout << " " << " " << N[i-1] << " " << E[i-1] << " " << P[i-1] << endl;
        }
    }
    return new KVDriver(E, num_lines, P, num_lines, N, num_lines);
}

int main(){
    cout << "HW" << endl;
    KVOR * m = new KVOR();
    cout << m->b << endl;
    double result[3];
    KVDriver * dr = parseEos("/home/const/workspace2/swigEosWrapper/Cuts/EPmpi4_KVOR06.dat");
    star_crust2(1., result, 3, dr, 1e-10);
}