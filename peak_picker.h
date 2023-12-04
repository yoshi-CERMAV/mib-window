//
//  peak_picker.h
//  
//
//  Created by Yoshiharu Nishiyama on 16/11/2023.
//

#ifndef peak_picker_h
#define peak_picker_h

#include <stdio.h>
int get_peak_list(double *data, int w, int h, int line_len, int box,
                  int *p, int max_size,  int &count);
#endif /* peak_picker_h */
