/*
 * pcr3b.cpp
 *
 *  Created on: Fri April 1 2016
 *      Author: zhengfaxiang
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "../include/rkf78.hpp"
using namespace std;

const long double pcr3bMiu = 0.1;
const long double pcr3bCj = 3.0996;
const long double pcr3bX0 = 0.65;
const long double pcr3bY0 = 0.0;
const long double pcr3bVx0 = 0.0;
const long double pcr3bTend = 200.0;


template<class T>
T pcr3br1(T x, T y) {
    return sqrt((x + pcr3bMiu) * (x + pcr3bMiu) + y * y);
}

template<class T>
T pcr3br2(T x, T y) {
    return sqrt((x + pcr3bMiu - 1) * (x + pcr3bMiu - 1)
                + y * y);
}

template<class T>
T pcr3bf0(T t, T y[4]) {
    return y[1];
}

template<class T>
T pcr3bf1(T t, T y[4]) {
    T R1 = pcr3br1(y[0], y[2]);
    T R2 = pcr3br2(y[0], y[2]);
    T tmp = (y[0] - (1 - pcr3bMiu) * (y[0] + pcr3bMiu)
             / pow(R1, 3) - pcr3bMiu * (y[0] + pcr3bMiu - 1)
             / pow(R2, 3));
    return 2 * y[3]  + tmp;
}

template<class T>
T pcr3bf2(T t, T y[4]) {
    return y[3];
}

template<class T>
T pcr3bf3(T t, T y[4]) {
    T R1 = pcr3br1(y[0], y[2]);
    T R2 = pcr3br2(y[0], y[2]);
    T tmp = (y[2] * (1 - (1 - pcr3bMiu) / pow(R1, 3) - pcr3bMiu
                     / pow(R2, 3)));
    return -2 * y[1] + tmp;
}

template<class T>
T pcr3bcj(T y[4]) {
    T R1 = pcr3br1(y[0], y[2]);
    T R2 = pcr3br2(y[0], y[2]);
    return (y[0] * y[0] + y[2] * y[2] + 2 * (1 - pcr3bMiu) / R1
            + 2 * pcr3bMiu / R2 - (y[1] * y[1] + y[3] * y[3]));
}


template<class T>
T pcr3bvy(T x, T y, T vx) {
    T R1 = pcr3br1(x, y);
    T R2 = pcr3br2(x, y);
    return sqrt(x * x + y * y + 2 * (1 - pcr3bMiu) / R1
                + 2 * pcr3bMiu / R2 - pcr3bCj - vx * vx);
}


int main(int argc, char *argv[]) {
    ofstream outfile;
    outfile.open("pcr3b.dat", ios::out);
    RKF78<long double, 4> RKF;
    RKF.f[0] = &pcr3bf0<long double>;
    RKF.f[1] = &pcr3bf1<long double>;
    RKF.f[2] = &pcr3bf2<long double>;
    RKF.f[3] = &pcr3bf3<long double>;
    long double Vy = pcr3bvy(pcr3bX0, pcr3bY0, pcr3bVx0);
    long double rkf[4] = {pcr3bX0, pcr3bVx0, pcr3bY0, Vy};
    long double h = 0.1;
    long double t = 0.0;
    for (; t < pcr3bTend;) {
        try {
            RKF.rkf78(h, t, rkf, 1e-1, 1e-6, 1e-12);
        } catch (invalid_argument& e) {
            cerr << e.what() << endl;
            return -1;
        }
        long double Cj = pcr3bcj(rkf);
        cout<<setiosflags(ios::scientific)
            <<setprecision(18)<<setw(28)<<t
            <<setw(28)<<h;
        cout<<setw(28)<<rkf[0]
            <<setw(28)<<rkf[2]
            <<setw(28)<<rkf[1]
            <<setw(28)<<Cj<<endl;
        outfile<<setiosflags(ios::scientific)
               <<setprecision(18)<<setw(28)<<t
               <<setw(28)<<h;
        outfile<<setw(28)<<rkf[0]
               <<setw(28)<<rkf[2]
               <<setw(28)<<rkf[1]
               <<setw(28)<<Cj<<endl;
    }
    long double Cj_end = pcr3bcj(rkf);
    long double Cj_relative_err = (Cj_end - pcr3bCj) / pcr3bCj;
    cout<<"Relative Cj error: "<<Cj_relative_err<<endl;
    outfile.close();
    cout<<"Procedure completed!"<<endl;
    return 0;
}
