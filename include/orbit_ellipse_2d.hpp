/*
 * orbit_ellipse_2d.hpp
 *
 *  Created on: Fri Mar 18 2016
 *      Author: zhengfAxiang
 *
 *  Last Modified: Sat Apr 2 2016
 *      Author: solarqiang
 */

#ifndef ORBIT_ELLIPSE_2D_HPP_
#define ORBIT_ELLIPSE_2D_HPP_

#include <cmath>
#include <stdexcept>

template<class T>
class OrbitEllipse2D{
private:
    T AngularMomentum();
    T Velocity();
    T SemiMajorAxis();
    T GetEnergy();
    T GetPeriod();
    T X;                        // x
    T Y;                        // y
    T V_x;                      // V x
    T V_y;                      // V y
    T R;                        // R
    T Ax;                       // Semimajor Axis
    T V_r;                      // radial velocity
    T Theta;                    // true anomaly
    T AngVel;                   // angular velocity
    T V;                        // velocity
    T AngMom;                   // angular momentum
    T Ener;                     // energy
    T Peri;                     // Period
    T Miu;                      // miu
    void XY2RTheta();
public:
    OrbitEllipse2D(T x, T y, T v_x, T v_y, T miu);
    T Energy(){return Ener;}
    T Period(){return Peri;}
};

// 
// template<class T>
// T OrbitEllipse2D<T>::Velocity() {
//     return sqrt(V_r * V_r + AngMom * AngMom / (R * R));
// }

template<class T>
T OrbitEllipse2D<T>::SemiMajorAxis() {
    return -Miu / (Ener * 2.0);
}

template<class T>
T OrbitEllipse2D<T>::GetEnergy() {
    T energy = V * V / 2.0 - Miu / R;
    if (energy > 0) {
        throw std::invalid_argument("Energy greater than zero!");
    }
    return energy;
}

template<class T>
T OrbitEllipse2D<T>::GetPeriod() {
    return 2 * atan(1) * 4 * sqrt(pow(Ax, 3.0) / Miu);
}

template<class T>
OrbitEllipse2D<T>::OrbitEllipse2D(T x, T y, T V_x, T V_y, T miu) {
    X = x, Y = y, V_x = V_x, V_y = V_y;
    Miu     = miu;

    R       = sqrt(X * X + Y * Y);
    V       = sqrt(V_x * V_x + V_y * V_y);
    Theta   = atan2(Y, X);
    T alpha = acos((X * V_x + Y * V_y) / (R * V));
    V_r     = V * cos(alpha);

    AngVel  = V * sin(alpha) / R;
    AngMom  = V * sin(alpha) * R;

    Ener    = GetEnergy();
    Ax      = SemiMajorAxis();
    Peri    = GetPeriod();
}

#endif  /* ORBIT_ELLIPSE_2D_HPP_ */
