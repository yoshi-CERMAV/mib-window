//
//  fit_all_frames.hpp
//  
//
//  Created by Yoshiharu Nishiyama on 26/11/2023.
//

#ifndef fit_all_frames_hpp
#define fit_all_frames_hpp
#include <Savitzky_Golay_2d.h>
#include <mib-read.h>
#include <stdio.h>
void fit_peaks_all_frames(mib *reader, char filename[], int mask_radius, int fit_radius, int frame0, int frame1, int frame_step, Filter2d *filter, int box_size);

#endif /* fit_all_frames_hpp */
