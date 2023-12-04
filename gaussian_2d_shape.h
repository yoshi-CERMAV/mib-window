//
//  gaussian_2d_shape.h
//  
//
//  Created by Yoshiharu Nishiyama on 23/11/2023.
//

#ifndef gaussian_2d_shape_h
#define gaussian_2d_shape_h
#include <math.h>

using namespace std;
const double gauss_width_cof = 2.7725887222397811449070559319807216525077819824219;
const double gauss_area_cof = 0.93943727869965132359908466241904534399509429931641;

class Gaussian2D_shape
{
public:
    Gaussian2D_shape(){}
    ~Gaussian2D_shape(){}
    void set(const double *p){
        pos_x = p[0]; pos_y = p[1]; width_x = p[2]; width_y = p[3];
    }
    void set(const double *p, const double *w){
        pos_x = p[0]; pos_y = p[1]; width_x = w[0]; width_y = w[1];
    }
    void set(const double &px, const double &py, const double &wx, const double &wy){
        pos_x = px;
        pos_y = py;
        width_x = wx;
        width_y = wy;
    }

    inline double jacobian(const double &x, const double &y, double *der)
    {
        double arg_x = (x-pos_x)/width_x;
        double arg2_x = arg_x*arg_x;
        double arg_y = (y-pos_y)/width_y;
        double arg2_y = arg_y*arg_y;
        double arg2 = arg2_x+arg2_y;
        double val = exp(-gauss_width_cof*arg2);
        double fac_x = 2 * arg_x / width_x;
        double fac_y = 2 * arg_y / width_y;
        double dpx = fac_x * gauss_width_cof * val;
        double dpy = fac_y * gauss_width_cof * val;
        der[0] = dpx;
        der[1] = dpy;
        der[2] = dpx * arg_x;
        der[3] = dpy * arg_y;
        return val;
    }
    
    inline double operator ()(const double &x, const double &y){
        double arg_x = (x-pos_x)/width_x;
        double arg2_x = arg_x*arg_x;
        double arg_y = (y-pos_y)/width_y;
        double arg2_y = arg_y*arg_y;
        double arg2 = arg2_x+arg2_y;
        return exp(-gauss_width_cof*arg2);
    }
private:
    double pos_x;
    double pos_y;
    double width_x;
    double width_y;
};

#endif /* gaussian_2d_shape_h */
