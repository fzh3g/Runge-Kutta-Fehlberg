/*
 * sitnikov.cpp
 *
 *  Created on: Sat April 3 2016
 *      Author: zhengfaxiang
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "../include/rkf78.hpp"
using namespace std;

long double sitnikovEccentric = 0.6;
long double sitnikovZ0        = 0.2;
long double sitnikovVz0       = 0.0;
long double sitnikovTend      = 60.0;


template<class T>
T sitnikovr0(T angle) {
    return ((-sitnikovEccentric * sitnikovEccentric + 1.0)
            / (sitnikovEccentric * cos(angle) + 1.0));
}

template<class T>
T sitnikovf0(T t, T y[3]) {
    return y[1];
}

template<class T>
T sitnikovf1(T t, T y[3]) {
    T R0 = sitnikovr0(y[2]);
    return -y[0] / pow(R0 * R0 / 4.0 + y[0] * y[0], 1.5);
}

template<class T>
T sitnikovf2(T t, T y[3]) {
    return (pow(sitnikovEccentric * cos(y[2]) + 1.0, 2.0)
            / pow(-pow(sitnikovEccentric, 2.0) + 1.0, 1.5));
}

int main(int argc, char *argv[]) {
    const char *datafile = "sitnikov.dat";
    RKF78<long double, 3> RKF;
    RKF.f[0] = &sitnikovf0<long double>;
    RKF.f[1] = &sitnikovf1<long double>;
    RKF.f[2] = &sitnikovf2<long double>;
    long double rkf[3] = {sitnikovZ0, sitnikovVz0, 0.0};
    try {
        RKF.solve(0.01, 1e-8, rkf, 1e-12, 0.0, sitnikovTend, datafile);
    } catch (invalid_argument& e) {
        cerr << e.what() << endl;
        return -1;
    }
    return 0;
}
