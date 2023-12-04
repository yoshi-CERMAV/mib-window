//
//  peak.h
//  
//
//  Created by Yoshiharu Nishiyama on 22/11/2023.
//

#ifndef diffraction_peak_h
#define diffraction_peak_h

typedef struct diffraction_peak{
    int frame;
    int ix, iy; // coordinate on the screen
    double x, y; // coordinate after centering
    double rx, ry, rz;// after rotation
    double wx, wy;//width;
    double I;
    double I_err;
} diffraction_peak;


#endif /* diffraction_peak_h */
