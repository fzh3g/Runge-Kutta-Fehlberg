/*
 * gaussperturb.cpp
 *
 *  Created on: <2016-05-04 Wed>
 *      Author: zhengfaxiang
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include "../include/rkf78.hpp"
using namespace std;

long double gaussperturb_f             = 0.0;
long double gaussperturb_e             = 0.5;
long double gaussperturb_a             = 1.0;
long double gaussperturb_c             = 1e-5;
long double gaussperturb_miu           = 1e-2;
long double gaussperturb_tol           = 1e-14;
long double gaussperturb_period_number = 5000;
long double gaussperturb_pace_init     = 1;
long double gaussperturb_pace_max      = 10;
long double gaussperturb_pace_min      = 1e-6;

double TimeDiff(timeval t1, timeval t2) {
    double t;
    t = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    t += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
    return t;
}

template<class T>
T gaussperturb_beta2(T e) {
    T tmp = 1.0 - e * e;
    if (tmp <= 0.0)
        tmp = 1e-12;
    return tmp;
}

template<class T>
T gaussperturb_gamma(T e) {
    return (1.0 / (1.0 + e * cos(gaussperturb_f)));
}

template<class T>
T gaussperturb_f0(T t, T y[4]) {
    return (-4.0 * gaussperturb_c * y[1] * sin(gaussperturb_f) /
            (pow(gaussperturb_miu, 0.5) * pow(y[0], 1.5) *
             pow(gaussperturb_gamma(y[1]), 3.0) * pow(gaussperturb_beta2(y[1]), 3.5)));
}

template<class T>
T gaussperturb_f1(T t, T y[4]) {
    return (-2.0 * gaussperturb_c * sin(gaussperturb_f) /
            (pow(gaussperturb_miu, 0.5) * pow(y[0], 2.5) *
             pow(gaussperturb_gamma(y[1]), 3.0) * pow(gaussperturb_beta2(y[1]), 2.5)));
}

template<class T>
T gaussperturb_f2(T t, T y[4]) {
    return (2.0 * gaussperturb_c * cos(gaussperturb_f) /
            (pow(gaussperturb_miu, 0.5) * pow(y[0], 2.5) * y[1] *
             pow(gaussperturb_gamma(y[1]), 3.0) * pow(gaussperturb_beta2(y[1]), 2.5)));
}

template<class T>
T gaussperturb_f3(T t, T y[4]) {
    return (pow(gaussperturb_miu, 0.5) / pow(y[0], 1.5) - 2.0 *gaussperturb_c /
            (pow(gaussperturb_miu, 0.5) * pow(y[0], 2.5) * y[1] *
             pow(gaussperturb_gamma(y[1]), 3.0) * pow(gaussperturb_beta2(y[1]), 2.0)) *
            (cos(gaussperturb_f) - 2.0 * gaussperturb_gamma(y[1]) * y[1]));
}

template<class T>
T newton_method(T M, T e, T TOL) {
    T En = M;
    T tmp = En;
    for(;;) {
        tmp = En;
        En = En - (En - e * sin(En) - M) /
            (1.0 - e * cos(En));
        T Err = abs(tmp - En);
        if (Err < TOL) 
            break;
    }
    return En;
}

template<class T>
void gaussperturb_get_f(T y[4], T TOL) {
    T E = newton_method<T>(y[3], y[1], TOL);
    long double cosf = ((gaussperturb_beta2(y[1]) / (1.0 - y[1] * cos(E)) - 1.0)
                        / y[1]);
    long double sinf = ((pow(gaussperturb_beta2(y[1]), 0.5) /
                         (1.0 - y[1] * cos(E))) * sin(E));
    gaussperturb_f = atan2(sinf, cosf);
}

int main(int argc, char *argv[]) {
    timeval t1, t2;
    gettimeofday(&t1, NULL);
    ofstream outfile;
    outfile.open("gauss_perturb.dat", ios::out);
    RKF78<long double, 4> RKF;
    RKF.f[0] = &gaussperturb_f0<long double>;
    RKF.f[1] = &gaussperturb_f1<long double>;
    RKF.f[2] = &gaussperturb_f2<long double>;
    RKF.f[3] = &gaussperturb_f3<long double>;
    long double a_e_omega_m[4] = {gaussperturb_a, gaussperturb_e, 0.0, 0.0};
    long double period = 8 * atan(1.0) * pow(a_e_omega_m[0], 1.5) /
        pow(gaussperturb_miu, 0.5);
    long double tend = gaussperturb_period_number * period;
    long double t = 0.0;
    long double h = gaussperturb_pace_init;
    int step = 0;
    for (; t < tend;) {
        gaussperturb_get_f<long double>(a_e_omega_m, gaussperturb_tol);
        try {
            RKF.rkf78(h, t, a_e_omega_m, gaussperturb_pace_max, gaussperturb_pace_min, gaussperturb_tol);
        } catch (invalid_argument& e) {
            cerr << e.what() << endl;
            return -1;
        }
        step++;
        long double r = (a_e_omega_m[0] *
                         gaussperturb_beta2(a_e_omega_m[1]) *
                         gaussperturb_gamma(a_e_omega_m[1]));
        long double x = r * cos(a_e_omega_m[2] + gaussperturb_f);
        long double y = r * sin(a_e_omega_m[2] + gaussperturb_f);
        cout<<setiosflags(ios::scientific)
            <<setprecision(18)<<setw(28)<<t
            <<setw(28)<<h
            <<setw(28)<<a_e_omega_m[0]
            <<setw(28)<<a_e_omega_m[1]
            <<setw(28)<<a_e_omega_m[2]
            <<setw(28)<<a_e_omega_m[3]<<endl;
        outfile<<setiosflags(ios::scientific)
               <<setprecision(18)<<setw(28)<<t
               <<setw(28)<<h
               <<setw(28)<<a_e_omega_m[0]
               <<setw(28)<<a_e_omega_m[1]
               <<setw(28)<<a_e_omega_m[2]
               <<setw(28)<<a_e_omega_m[3]
               <<setw(28)<<x
               <<setw(28)<<y
               <<setw(28)<<gaussperturb_f<<endl;
    }
    outfile.close();
    gettimeofday(&t2, NULL);
    double d = TimeDiff(t1,t2);
    cout<<"Time spent: "<<d<<" ms"<<endl;
    cout<<step<<" steps in total!"<<endl;
    cout<<"Procedure completed!"<<endl;
    return 0;
}

