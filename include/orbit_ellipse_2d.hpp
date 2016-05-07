/*
 * orbit_ellipse_2d.hpp
 *
 *  Created on: Fri Mar 18 2016
 *      Author: zhengfaxiang
 */

#ifndef ORBIT_ELLIPSE_2D_HPP_
#define ORBIT_ELLIPSE_2D_HPP_

#include <cmath>
#include <stdexcept>

template<class T>
class OrbitEllipse2D{
private:
    T AngularMomentum(T R, T AngVel);
    T Velocity(T V_r, T AngMom, T R);
    T SemiMajorAxis(T Ener, T miu);
    T Energy(T V, T miu, T R);
    T Period(T SemiMajorAxis, T miu);
    T GetEnergy();
    T GetPeriod();
    void XY2RTheta(T x, T y, T v_x, T v_y);
public:
    OrbitEllipse2D(T x, T y, T v_x, T v_y, T miu);
    T X;                        // x
    T Y;                        // y
    T V_x;                      // V x
    T V_y;                      // V y
    T R;                        // R
    T V_r;                      // radial velocity
    T Theta;                    // true anomaly
    T AngVel;                   // angular velocity
    T V;                        // velocity
    T AngMom;                   // angular momentum
    T Ener;                     // energy
    T Peri;                     // Period
    T Miu;                      // miu
};

template<class T>
T OrbitEllipse2D<T>::AngularMomentum(T R, T AngVel) {
    return AngVel * R * R;
}

template<class T>
T OrbitEllipse2D<T>::Velocity(T V_r, T AngMom, T R) {
    return sqrt(V_r * V_r + AngMom * AngMom / (R * R));
}

template<class T>
T OrbitEllipse2D<T>::SemiMajorAxis(T Ener, T miu) {
    return -miu / (Ener * 2.0);
}

template<class T>
T OrbitEllipse2D<T>::Energy(T V, T miu, T R) {
    T energy = V * V / 2.0 - miu / R;
    if (energy > 0) {
        throw std::invalid_argument("Energy greater than zero!");
    }
    return energy;
}

template<class T>
T OrbitEllipse2D<T>::Period(T SemiMajorAxis, T miu) {
    return sqrt(pow(SemiMajorAxis, 3.0) / miu) * 8.0 * atan(1.0);
}

template<class T>
T OrbitEllipse2D<T>::GetEnergy() {
    T velo = Velocity(V_r, AngMom, R);
    return Energy(velo, Miu, R);
}

template<class T>
T OrbitEllipse2D<T>::GetPeriod() {
    T ener = GetEnergy();
    T a = SemiMajorAxis(ener, Miu);
    return Period(a, Miu);
}

template<class T>
void OrbitEllipse2D<T>::XY2RTheta(T x, T y, T v_x, T v_y) {
    R = sqrt(x * x + y * y);
    V = sqrt(v_x * v_x + v_y * v_y);
    Theta = atan2(y, x);
    T alpha = acos((x * v_x + y * v_y) / (R * V));
    V_r = V * cos(alpha);
    AngVel = V * sin(alpha) / R;
}

template<class T>
OrbitEllipse2D<T>::OrbitEllipse2D(T x, T y, T v_x, T v_y, T miu) {
    XY2RTheta(x, y, v_x, v_y);
    AngMom = AngularMomentum(R, AngVel);
    Miu = miu;
    Ener = GetEnergy();
    Peri = GetPeriod();
}

#endif  /* ORBIT_ELLIPSE_2D_HPP_ */
