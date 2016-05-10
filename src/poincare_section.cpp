/*
 * poincare_section.cpp
 *
 *  Created on: Sun April 2 2016
 *      Author: zhengfaxiang
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "../include/rkf78.hpp"
using namespace std;

const long double pcr3bCj    = 2.999 - 0.1;
const long double pcr3bMiu   = 0.001;
const long double pcr3bCjTOL = 1e-11;
const long double pcr3bTOL   = 1e-12;
const long double pcr3bBegin = 0.45;
const long double pcr3bEnd   = 0.55;
const long double pcr3bTend  = 800.0;
const int pcr3bNum           = 30;


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
T pcr3bvysquare(T x, T y, T vx) {
    T R1 = pcr3br1(x, y);
    T R2 = pcr3br2(x, y);
    return (x * x + y * y + 2 * (1 - pcr3bMiu) / R1
            + 2 * pcr3bMiu / R2 - pcr3bCj - vx * vx);
}

int main(int argc, char *argv[]) {
    ofstream outfile;
    outfile.open("poincare_section.dat", ios::out);
    RKF78<long double, 4> RKF;
    RKF.f[0] = &pcr3bf0<long double>;
    RKF.f[1] = &pcr3bf1<long double>;
    RKF.f[2] = &pcr3bf2<long double>;
    RKF.f[3] = &pcr3bf3<long double>;
    long double gap = (pcr3bEnd - pcr3bBegin) / pcr3bNum;
    long double X = 0.0;
    long double Y = 0.0;
    long double Vx = 0.0;
    long double Vy = 0.0;
    for (int i = 0; i < pcr3bNum; i++) {
        X = pcr3bBegin + gap * i;
        long double VySquare = pcr3bvysquare(X, Y, Vx);
        if (VySquare <= 0)
            continue;
        else
            Vy = sqrt(VySquare);
        long double rkf[4] = {X, Vx, Y, Vy};
        long double h = 0.1;
        long double t = 0.0;
        for (; t < pcr3bTend;) {
            try {
                RKF.rkf78(h, t, rkf, 10, 1e-6, pcr3bTOL);
            } catch (invalid_argument& e) {
                cerr << e.what() << endl;
                break;
            }
            long double Cj = pcr3bcj(rkf);
            long double Cj_relative_err = (Cj - pcr3bCj) / pcr3bCj;
            if (abs(Cj_relative_err) > pcr3bCjTOL) {
                break;
            }
            cout<<i<<setiosflags(ios::scientific)
                <<setprecision(18)<<setw(28)<<t
                <<setw(28)<<h;
            cout<<setw(28)<<rkf[0]
                <<setw(28)<<rkf[2]
                <<setw(28)<<rkf[1]<<endl;
            outfile<<setiosflags(ios::scientific)
                   <<setprecision(18)<<setw(28)<<t
                   <<setw(28)<<h;
            outfile<<setw(28)<<rkf[0]
                   <<setw(28)<<rkf[2]
                   <<setw(28)<<rkf[1]<<endl;
        }
    }
    outfile.close();
    cout<<"Procedure completed!"<<endl;
    return 0;
}
