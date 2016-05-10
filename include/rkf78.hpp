/*
 * rkf78.hpp
 *
 *  Created on: Wed Mar 16 2016
 *      Author: zhengfaxiang
 */

#ifndef RKF78_HPP_
#define RKF78_HPP_

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <sys/time.h>
using namespace std;

template<class T, int dim>
class RKF78{
private:
    static const T a[13];       // rk78 parameters
    static const T b[13][12];   // rk78 parameters
    T K[13][dim];               // runge-kutta parameters
    T R[dim];                   // error
    T z[13][dim];               // make calculation of Ks easier
    double TimeDiff(timeval t1, timeval t2);
    void GetZ(int l, T y[dim]);
    void RungeKuttaParams78(T t, T h, T y[dim]);
    void GetY(T y[dim]);
public:
    T (*f[dim])(T t, T y[dim]); // functions to solve
    void rkf78(T& h, T& t, T ynow[dim], T hmax, T hmin, T TOL);
    void solve(T hinit, T hmin, T y0[dim], T TOL, T begin, T end,
               const char *filename);
};

template<class T, int dim>
const T RKF78<T, dim>::a[13] = {
    0.0, 2.0/27.0, 1.0/9.0, 1.0/6.0, 5.0/12.0, 0.5,
    5.0/6.0, 1.0/6.0, 2.0/3.0, 1.0/3.0, 1.0, 0.0, 1.0
};

template<class T, int dim>
const T RKF78<T, dim>::b[13][12] = {
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {2.0/27.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {1.0/36.0, 1.0/12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {1.0/24.0, 0.0, 1.0/8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {5.0/12.0, 0.0, -25.0/16.0, 25.0/16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, },
    {1.0/20.0, 0.0, 0.0, 1.0/4.0, 1.0/5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0},
    {-25.0/108.0, 0.0, 0.0, 125.0/108.0, -65.0/27.0, 125.0/54.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0},
    {31.0/300, 0.0, 0.0, 0.0, 61.0/225.0, -2.0/9.0, 13.0/900.0, 0.0,
     0.0, 0.0, 0.0, 0.0},
    {2.0, 0.0, 0.0, -53.0/6.0, 704.0/45.0, -107.0/9.0, 67.0/90.0, 3.0,
     0.0, 0.0, 0.0, 0.0},
    {-91.0/108.0, 0.0, 0.0, 23.0/108.0, -976.0/135.0, 311.0/54.0, -19.0/60.0,
     17.0/6.0, -1.0/12.0, 0.0, 0.0, 0.0},
    {2383.0/4100, 0.0, 0.0, -341.0/164.0, 4496.0/1025.0, -301.0/82.0,
     2133.0/4100.0, 45.0/82.0, 45.0/164.0, 18.0/41.0, 0.0, 0.0},
    {3.0/205.0, 0.0, 0.0, 0.0, 0.0, -6.0/41.0, -3.0/205.0, -3.0/41.0,
     3.0/41.0, 6.0/41.0, 0.0, 0.0},
    {-1777.0/4100.0, 0.0, 0.0, -341.0/164.0, 4496.0/1025.0, -289.0/82.0,
     2193.0/4100.0, 51.0/82.0, 33.0/164.0, 12.0/41.0, 0.0, 1.0}
};

template<class T, int dim>
double RKF78<T, dim>::TimeDiff(timeval t1, timeval t2) {
    double t;
    t = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    t += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
    return t;
}

template<class T, int dim>
void RKF78<T, dim>::GetZ(int l, T y[dim]) {
    for (int i=0; i < dim; i++) {
        z[l][i] = y[i];
        if (l != 0) {
            for (int j=0; j < l; j++) {
                z[l][i] += K[j][i] * b[l][j];
            }
        }
    }
}

template<class T, int dim>
void RKF78<T, dim>::RungeKuttaParams78(T t, T h, T y[dim]) {
    // get runge-kutta parameters
    for (int j=0; j < 13; j++) {
        GetZ(j, y);
        for (int i=0; i < dim; i++) {
            K[j][i] = h * (*f[i])(t + h * a[j], z[j]);
        }
    }
}

template<class T, int dim>
void RKF78<T, dim>::GetY(T y[dim]) {
    // get y using runge-kutta parameters
    for (int i=0; i < dim; i++) {
        y[i] += (K[5][i] * 34.0 / 105.0 +
                 (K[6][i] + K[7][i]) * 9.0 / 35.0 +
                 (K[8][i] + K[9][i]) * 9.0 / 280.0 +
                 (K[11][i] + K[12][i]) * 41.0 / 840.0);
    }
}

template<class T, int dim>
void RKF78<T, dim>::rkf78(T& h, T& t, T y[dim], T hmax, T hmin, T TOL) {
    // function to apply the runge-kutta-fehlberg method for one step
    if (hmin == hmax) {
        RungeKuttaParams78(t, h, y);  // get Ks
        GetY(y);                      // get Ys
        t += h;
    }
    else {
        for (;;) {
            RungeKuttaParams78(t, h, y);  // get Ks
            for (int i=0; i < dim; i++) { // finding errors
                R[i] = (abs(K[0][i] + K[10][i] - K[11][i] - K[12][i])
                        / (TOL * h) * 41.0 / 810.0);
            }
            T MaxErr = *max_element(R, R + dim); // maximium value of R
            if (MaxErr < 1) {
                GetY(y);        // get Ys
                t += h;
                if (MaxErr < 0.1) { // steps too small
                    if (MaxErr == 0) {
                        h = hmax;
                    }
                    else {
                        h = min(h * 2.0, hmax);
                    }
                }
                break;
            }
            else {              // error not tolerable
                if (h == hmin){
                    throw invalid_argument("Minimum h exceeded!");
                }
                h = max(h / 2.0, hmin);
            }
        }
    }
}

template<class T, int dim>
void RKF78<T, dim>::solve(T hinit, T hmin, T y[dim], T TOL, T begin, T end,
           const char *filename) {
    // function to apply the runge-kutta-fehlberg method
    // time
    timeval t1, t2;
    gettimeofday(&t1, NULL);
    // open a file in write mode.
    ofstream outfile;
    outfile.open(filename, ios::out);
    T t = begin;                // begin of t
    T h = hinit;                // begin of h
    long step = 0;                   // step
    // output header
    cout<<setw(28)<<"t"<<setw(28)<<"h";
    outfile<<setw(28)<<"t"<<setw(28)<<"h";
    for (int i=0; i < dim; i++) {
        // convert number to string
        string yi = "y" + to_string(i);
        cout<<setw(28)<<yi;
        outfile<<setw(28)<<yi;
    }
    cout<<endl;
    outfile<<endl;
    for (;t < end;) {
        rkf78(h, t, y, hinit, hmin, TOL); // calculate one step
        step++;                           // step plus one
        // output result
        cout<<setiosflags(ios::scientific)
            <<setprecision(18)<<setw(28)<<t
            <<setw(28)<<h;
        outfile<<setiosflags(ios::scientific)
               <<setprecision(18)<<setw(28)<<t
               <<setw(28)<<h;
        for (int i=0; i < dim; i++) {
            cout<<setiosflags(ios::scientific)<<setprecision(18)
                <<setw(28)<<y[i];
            outfile<<setiosflags(ios::scientific)<<setprecision(18)
                   <<setw(28)<<y[i];
        }
        cout<<endl;
        outfile<<endl;
    }
    outfile.close();                    // close the opened file.
    gettimeofday(&t2, NULL);
    double d = TimeDiff(t1,t2);
    cout<<"Procedure completed!"<<endl;
    cout<<"Total steps: "<<step<<endl;
    cout<<"Time spent: "<<d<<" ms"<<endl;
}

#endif  /* RKF78_HPP_ */
