/*
 * twobody.cpp
 *
 *  Created on: Wed Mar 16 2016
 *      Author: zhengfaxiang
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "../include/rkf78.hpp"
#include "../include/orbit_ellipse_2d.hpp"
using namespace std;

// initialize the position and velocity of the orbit
// Miu
const long double twobody_Miu            = 10.0;
// X
const long double twobody_X              = 1.1;
// Y
const long double twobody_Y              = 0.0;
// V_x
const long double twobody_V_x            = 0.0;
// V_y
const long double twobody_V_y            = 1.0;
// number of period
const long double twobody_NumberofPeriod = 2;
// TOL
const long double twobody_TOL            = 1e-12;

// function templates
template<class T>
T f0(T t, T y[4]) {
    return y[2];
}
template<class T>
T f1(T t, T y[4]) {
    return y[3];
}
template<class T>
T f2(T t, T y[4]) {
    return -twobody_Miu * y[0] * pow(y[0] * y[0] + y[1] * y[1], -1.5);
}
template<class T>
T f3(T t, T y[4]) {
    return -twobody_Miu * y[1] * pow(y[0] * y[0] + y[1] * y[1], -1.5);
}

int main() {
    const char *outputfile = "twobody_output.dat"; // output file
    long double rkf[4] = {twobody_X, twobody_Y, twobody_V_x,
                          twobody_V_y}; // variables for rkf
    // initial orbit
    long double energy_begin = 0.0;
    long double period_begin = 0.0;
    try {                       // error handling
        OrbitEllipse2D<long double> OrbElli_begin(rkf[0], rkf[1],
                                                  rkf[2], rkf[3],
                                                  twobody_Miu);
        energy_begin = OrbElli_begin.Ener; // energy
        period_begin = OrbElli_begin.Peri; // period
    } catch (invalid_argument& e) {
        cerr << e.what() << endl;
        return -1;
    }
    RKF78<long double, 4>  RKF; // rkf78-4
    // pass functions addresses to rkf78-4
    RKF.f[0] = &(f0<long double>);
    RKF.f[1] = &(f1<long double>);
    RKF.f[2] = &(f2<long double>);
    RKF.f[3] = &(f3<long double>);
    long double tend = period_begin * twobody_NumberofPeriod; // tend
    try {                       // error handling
        RKF.solve(period_begin / 10.0, period_begin / 1e+8, rkf,
                  twobody_TOL, 0.0, tend, outputfile); // apply rkf78-4
    } catch (invalid_argument& e) {
        cerr << e.what() << endl;
        return -1;
    }
    long double energy_end = 0.0;
    long double period_end = 0.0;
    // final orbit
    try {                       // error handling
        OrbitEllipse2D<long double> OrbElli_end(rkf[0], rkf[1],
                                                rkf[2], rkf[3],
                                                twobody_Miu);
        energy_end = OrbElli_end.Ener; // energy
        period_end = OrbElli_end.Peri; // period
    } catch (invalid_argument& e) {
        cerr << e.what() << endl;
        return -1;
    }
    long double energy_err = energy_end - energy_begin; // error of energy
    long double period_err = period_end - period_begin; // error of period
    // relative errors
    long double energy_relative_err = abs(energy_err / energy_begin);
    long double period_relative_err = abs(period_err / period_begin);
    // write to ascii file
    ofstream outfile;
    outfile.open("energy_period_errs.txt");
    outfile << "Energy begin: " << setiosflags(ios::scientific)
            << setprecision(18) << energy_begin << endl
            << "Energy end: " << energy_end << endl
            << "Energy error: " << energy_err << endl
            << "Energy relative error: " << resetiosflags(ios::scientific)
            << energy_relative_err * 100.0 << "\%" << endl
            << "Period begin: " << setiosflags(ios::scientific)
            << period_begin << endl
            << "Period end: " << period_end << endl
            << "Period error: " << period_err << endl
            << "Period relative error: " << resetiosflags(ios::scientific)
            << period_relative_err * 100.0 << "\%" << endl;
    outfile.close();
    // output to screen
    cout << "Energy begin: " << setiosflags(ios::scientific)
         << setprecision(18) << energy_begin << endl
         << "Energy end: " << energy_end << endl
         << "Energy error: " << energy_err << endl
         << "Energy relative error: " << resetiosflags(ios::scientific)
         << energy_relative_err * 100.0 << "\%" << endl
         << "Period begin: " << setiosflags(ios::scientific)
         << period_begin << endl
         << "Period end: " << period_end << endl
         << "Period error: " << period_err << endl
         << "Period relative error: " << resetiosflags(ios::scientific)
         << period_relative_err * 100.0 << "\%" << endl;
    return 0;
}
