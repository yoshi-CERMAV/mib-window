//
//  fit_2d_peak.c
//  
//
//  Created by Yoshiharu Nishiyama on 22/11/2023.
//

#include "fit_2d_peak.h"
#include <Accelerate/Accelerate.h>
#include "gaussian_2d_shape.h"
#include <mpfit.h>
#include <iostream>
using namespace std;
static double residual[2500];
int fit_func(int m, int n, double *p, double *dy, double **dvec, void *vars ){
    double temp[8];
    double *xy = reinterpret_cast<double *>(vars);
    double *I = xy + 2*m;
    //p0-2 for flat bg, p3, 4 for two gaussian height, p5
    Gaussian2D_shape G1, G2;
    G1.set(p+5, p+7);
    G2.set(p[5], p[6], p[7]*2, p[8]*2);
    int n3 = 3;
    int one = 1;
    fill(dy, dy+m, p[0]);

    cblas_dgemv(CblasColMajor, CblasNoTrans,
                m, 2, 1, xy, m, p+1, 1, 1, dy, 1);
    if(dvec){
        if(dvec[0]) fill(dvec[0], dvec[0] + m, 1. );
        if(dvec[1]) memcpy(dvec[1], xy, m*sizeof(double));
        if(dvec[2]) memcpy(dvec[2], xy+m, m*sizeof(double));
        for(int i = 0; i < m; i++){
            double x = xy[i];
            double y = xy[i+m];
            double g1 = G1.jacobian(x,y, temp) ;
            double g2 = G2.jacobian(x,y, temp+4);
            dy[i] += p[3] * g1;
            dy[i] += p[4] * g2;
            if(dvec[3]) dvec[3][i] = g1;
            if(dvec[4]) dvec[4][i] = g2;
            if(dvec[5]) dvec[5][i] = temp[0] + temp[4];
            if(dvec[6]) dvec[6][i] = temp[1] + temp[5];
            if(dvec[7]) dvec[7][i] = temp[2] + temp[6]*0.5;
            if(dvec[8]) dvec[8][i] = temp[3] + temp[7]*0.5;
            dy[i] -= I[i];
        }
        return 0;
    }else{
         for(int i = 0; i < m; i++){
            dy[i] += p[0];
            double x = xy[i];
            double y = xy[i+m];
            double g1 = G1.jacobian(x,y, temp) ;
            double g2 = G2.jacobian(x,y, temp+4);
            dy[i] += p[3] * g1;
            dy[i] += p[4] * g2;
            dy[i] -= I[i];
        }
        return 0;
    }
    return 1;
}

int fit_func1(int m, int n, double *p, double *dy, double **dvec, void *vars ){
    assert(n==8);
    double temp[8];
    double *xy = reinterpret_cast<double *>(vars);
    double *I = xy + 2*m;
    //p0-2 for flat bg, p3, 4 for two gaussian height, p5
    Gaussian2D_shape G1, G2;
    G1.set(p+4, p+6);
    G2.set(p[4], p[5], p[6]*2, p[7]*2);
    int n3 = 3;
    int one = 1;
    fill(dy, dy+m, p[0]);

    cblas_dgemv(CblasColMajor, CblasNoTrans,
                m, 2, 1, xy, m, p+1, 1, 1, dy, 1);
    if(dvec){
        if(dvec[0]) fill(dvec[0], dvec[0] + m, 1. );
        if(dvec[1]) memcpy(dvec[1], xy, m*sizeof(double));
        if(dvec[2]) memcpy(dvec[2], xy+m, m*sizeof(double));
        for(int i = 0; i < m; i++){
            double x = xy[i];
            double y = xy[i+m];
            double g1 = G1.jacobian(x,y, temp) ;
            double g2 = G2.jacobian(x,y, temp+4);
            double g = g1+g2;
            dy[i] += p[3] * g;
  //          dy[i] += p[3] * g2;
            if(dvec[3]) dvec[3][i] = g1;
  //          if(dvec[4]) dvec[4][i] = g2;
            if(dvec[4]) dvec[4][i] = temp[0] + temp[4];
            if(dvec[5]) dvec[5][i] = temp[1] + temp[5];
            if(dvec[6]) dvec[6][i] = temp[2] + temp[6]*0.5;
            if(dvec[7]) dvec[7][i] = temp[3] + temp[7]*0.5;
            dy[i] -= I[i];
        }
        return 0;
    }else{
         for(int i = 0; i < m; i++){
            dy[i] += p[0];
            double x = xy[i];
            double y = xy[i+m];
            double g1 = G1.jacobian(x,y, temp) ;
            double g2 = G2.jacobian(x,y, temp+4);
             double g = g1+g2;
            dy[i] += p[3] * g;
//            dy[i] += p[4] * g2;
            dy[i] -= I[i];
        }
        return 0;
    }
    return 1;
}

void fit(double p[], double p_error[], double *data, int data_size)
{
    mp_par pars[9];
    memset(pars, 0, sizeof(pars));
    mp_result result;
    memset(&result, 0, sizeof(result));
    result.resid = residual;
    result.xerror = p_error;
    int status;
    p[0] = 0;// const
    p[1] = 0;// *x
    p[2] = 0;// *y
    p[3] = 1;// peak_height1
    p[4] = 1;//peak_height2
    p[5] = 0; //pos_x
    p[6] = 0; //pos_y
    p[7] = 1.5; //width_x
    p[8] = 1.5; //width_y
    for(int i = 3; i < 5; i++){
        pars[i].limited[0] = 1;
        pars[i].limits[0] = 0;
    }
    for(int i = 5; i < 9; i++) pars[i].fixed = 1;
//        cout << "start fit" <<" "<< data_size<<endl;
    status = mpfit(fit_func, data_size, 9, p, pars, 0, data, &result);
//        cout << "status 1 :"<< status<<endl;
    for(int i = 5; i < 9; i++)pars[i].fixed = 0;
    for(int i = 5; i < 7; i++){
        pars[i].limited[0] = 1;
        pars[i].limited[1] = 1;
        pars[i].limits[0] = -2;
        pars[i].limits[1] = 2;
    }
    for(int i = 7; i < 9; i++){
        pars[i].limited[0] = 1;
        pars[i].limited[1] = 1;
        pars[i].limits[0] = 1;
        pars[i].limits[1] = 5;
    }
//    cout << "starting second fit "<<endl;
    status = mpfit(fit_func, data_size, 9, p, pars, 0, data, &result);
    for(int i = 0 ; i < 3; i++) pars[i].fixed=1;
    for(int i = 7;i < 9; i++) pars[i].fixed = 1;
    status = mpfit(fit_func, data_size, 9, p, pars, 0, data, &result);
    
    
//    for(int i = 0; i < data_size; i++){
//        cout << data[i]<< " "<<data[i+data_size]<<" "<<data[i+data_size*2]<<" "<<residual[i]<<endl;
//    }
 //   cout << "status 2 :"<< status<<endl;
 //   cout << p[3]<<" "<<p[4] <<" "<<p[5]<<endl;
//   plot_func_2dspot(p);
}

void fit1(double p[], double p_error[], double *data, int data_size)
{
    mp_par pars[8];
    memset(pars, 0, sizeof(pars));
    mp_result result;
    memset(&result, 0, sizeof(result));
    result.resid = residual;
    result.xerror = p_error;
    int status;
    p[0] = 0;// const
    p[1] = 0;// *x
    p[2] = 0;// *y
    p[3] = 1;// peak_height1
    p[4] = 0; //pos_x
    p[5] = 0; //pos_y
    p[6] = 1.5; //width_x
    p[7] = 1.5; //width_y
    pars[3].limited[0] = 1;
    pars[3].limits[0] = 0;
    for(int i = 4; i < 8; i++) pars[i].fixed = 1;
//        cout << "start fit" <<" "<< data_size<<endl;
    status = mpfit(fit_func1, data_size, 8, p, pars, 0, data, &result);
//        cout << "status 1 :"<< status<<endl;
    for(int i = 4; i < 8; i++)pars[i].fixed = 0;
    for(int i = 4; i < 6; i++){
        pars[i].limited[0] = 1;
        pars[i].limited[1] = 1;
        pars[i].limits[0] = -2;
        pars[i].limits[1] = 2;
    }
    for(int i = 6; i < 8; i++){
        pars[i].limited[0] = 1;
        pars[i].limited[1] = 1;
        pars[i].limits[0] = 1;
        pars[i].limits[1] = 5;
    }
//    cout << "starting second fit "<<endl;
    status = mpfit(fit_func1, data_size, 8, p, pars, 0, data, &result);
    for(int i = 0 ; i < 3; i++) pars[i].fixed=1;
    for(int i = 6;i < 8; i++) pars[i].fixed = 1;
    status = mpfit(fit_func1, data_size, 8, p, pars, 0, data, &result);
    
    
//    for(int i = 0; i < data_size; i++){
//        cout << data[i]<< " "<<data[i+data_size]<<" "<<data[i+data_size*2]<<" "<<residual[i]<<endl;
//    }
 //   cout << "status 2 :"<< status<<endl;
 //   cout << p[3]<<" "<<p[4] <<" "<<p[5]<<endl;
//   plot_func_2dspot(p);
}
