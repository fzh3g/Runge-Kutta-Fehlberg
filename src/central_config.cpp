/*
 * central_config.cpp
 *
 *  Created on: Mon Mar 28 2016
 *      Author: zhengfaxiang
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "../include/rkf78.hpp"
using namespace std;

long double cent_con_phi0    = 1.0;
long double cent_con_phidot0 = 1.0;
long double cent_con_lambda  = 1.0;

template<class T>
T f0(T t, T y[2]) {
    return y[1];
}

template<class T>
T f1(T t, T y[2]) {
    return -cent_con_lambda / (y[0] * y[0]);
}

int main(int argc, char *argv[]) {
    ofstream outfile;
    outfile.open("central_config.dat", ios::out);
    RKF78<long double, 2> RKF;
    RKF.f[0] = &f0<long double>;
    RKF.f[1] = &f1<long double>;
    long double h = 0.01;
    long double t = 0.0;
    long double tend = 20.0;
    long double index[10] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2,
                             1.4, 1.6, 1.8, 2.0};
    for (int i = 0; i < 10; i++) {
        cent_con_phi0 = index[i];
        h = 0.01;
        t = 0.0;
        long double rkf[2] = {cent_con_phi0, cent_con_phidot0};
        for (; t < tend;) {
            try {
                RKF.rkf78(h, t, rkf, 1e-2, 1e-30, 1e-2);
            } catch (invalid_argument& e) {
                cerr << e.what() << endl;
                return -1;
            }
            cout<<setiosflags(ios::scientific)
                <<setprecision(18)<<setw(28)<<t
                <<setw(28)<<h
                <<setw(28)<<rkf[0]
                <<setw(28)<<rkf[1]<<endl;
            outfile<<setiosflags(ios::scientific)
                   <<setprecision(18)<<setw(28)<<t
                   <<setw(28)<<h
                   <<setw(28)<<rkf[0]
                   <<setw(28)<<rkf[1]<<endl;
        }
    }
    outfile.close();
    cout<<"Procedure completed!"<<endl;
    return 0;
}
