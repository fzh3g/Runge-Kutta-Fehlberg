/*
 * spin_orbit_resonance.cpp
 *
 *  Created on: <2016-05-21 Sat>
 *      Author: zhengfaxiang
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sys/time.h>
#include "../include/rkf78.hpp"
using namespace std;

static const long double pi          = atan(1.0) * 4.0;
static const long double sor_e       = 0.03;
static const long double sor_epsilon = 0.05;
static const long double sor_n       = 1.0;
static const long double sor_begin   = 0.0;
static const long double sor_end     = 2.2;
static const long double sor_tol     = 1e-16;
static const int sor_period_num      = 40;
static const int sor_thetadot_num    = 12;
static const int sor_theta_num       = 10;

double TimeDiff(timeval t1, timeval t2) {
    double t;
    t = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    t += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
    return t;
}

template<class T>
T NearestDiff(T x, T a) {
    T diff = x - floor(x / a) * a;
    if (diff >= a / 2)
        return diff - a;
    else
        return diff;
}

template<class T>
T sor_f0(T t, T y[2]) {
    return y[1];
}

template<class T>
T sor_f1(T t, T y[2]) {
    return (-sor_epsilon * sor_n * sor_n *
            (sin((y[0] - sor_n * t) * 2.0) -
             sor_e * (sin(y[0] * 2.0 - sor_n * t) -
                      sin(y[0] * 2.0 - sor_n * t * 3.0) * 7.0) / 2.0));
}

int main(int argc, char *argv[]) {
    timeval t1, t2;
    gettimeofday(&t1, NULL);
    ofstream outfile;
    outfile.open("spin_orbit_resonance.dat", ios::out);
    RKF78<long double, 2> RKF;
    RKF.f[0] = &sor_f0<long double>;
    RKF.f[1] = &sor_f1<long double>;
    long double gap_thetadot = (sor_end - sor_begin) / sor_thetadot_num;
    long double gap_theta = pi * 2.0 / sor_theta_num;
    long double tend = 2.0 * pi / sor_n * sor_period_num;
    int totstep = sor_theta_num * sor_thetadot_num;
    int step = 0;
    for (int i = 0; i < sor_thetadot_num; i++) {
        for (int j = 0; j < sor_theta_num; j++) {
            step++;
            int process = 100 * step / totstep;
            long double thetadot = sor_begin + gap_thetadot * i;
            long double theta = gap_theta * j;
            long double t = 0.0;
            long double h = 0.1;
            long double rkf[2] = {theta, thetadot};
            long double rkf1[2] = {theta, thetadot};
            for (; t < tend;) {
                rkf1[0] = rkf[0];
                rkf1[1] = rkf[1];
                long double tdiffbefore = NearestDiff<long double>(t, pi * 2.0);
                try {
                    RKF.rkf78(h, t, rkf, 1, 1e-6, sor_tol);
                } catch (invalid_argument& e) {
                    cerr << e.what() << endl;
                    break;
                }
                long double tdiffafter = NearestDiff<long double>(t, pi * 2.0);
                if (tdiffbefore < 0 && tdiffafter > 0) {
                    outfile << setiosflags(ios::scientific)
                            << setprecision(18)
                            << setw(28)
                            << fmodl((rkf1[0] + rkf[0]) / pi / 4, 1)
                            << setw(28) << (rkf1[1] + rkf[1]) / 2
                            << endl;
                }
                else if (tdiffafter == 0 && tdiffbefore < 0) {
                    outfile << setiosflags(ios::scientific)
                            << setprecision(18)
                            << setw(28)
                            << fmodl(rkf[0] / pi / 2, 1)
                            << setw(28) << rkf[1]
                            << endl;
                }
                cout << setw(3) << process << " %" << endl;
            }
        }
    }
    outfile.close();
    gettimeofday(&t2, NULL);
    double d = TimeDiff(t1,t2);
    cout << "Time spent: " << d << " ms" << endl;
    cout << "Procedure completed!" << endl;
    return 0;
}
