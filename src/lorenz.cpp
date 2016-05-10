/*
 * lorenz.cpp
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


template<class T>
T lorenzf0(T t, T y[3]) {
    return -10 * y[0] + 10 * y[1];
}

template<class T>
T lorenzf1(T t, T y[3]) {
    return 28 * y[0] - y[1] - y[0] * y[2];
}

template<class T>
T lorenzf2(T t, T y[3]) {
    return -8/3 * y[2] + y[0] * y[1];
}

int main(int argc, char *argv[]) {
    const char *datafile = "lorenz.dat";    
    RKF78<long double, 3> RKF;
    RKF.f[0] = &lorenzf0<long double>;
    RKF.f[1] = &lorenzf1<long double>;
    RKF.f[2] = &lorenzf2<long double>;
    long double rkf[3] = {1.0, 1.0, 0.0};
    try {
        RKF.solve(0.1, 1e-8, rkf, 1e-12, 0.0, 50, datafile);
    } catch (invalid_argument& e) {
        cerr << e.what() << endl;
        return -1;
    }
    return 0;
}

