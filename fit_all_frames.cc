//
//  fit_all_frames.cpp
//  
//
//  Created by Yoshiharu Nishiyama on 26/11/2023.
//

//#include "fit_all_frames.h"
#include "diffraction_data.h"
#include <Savitzky_Golay_2d.h>
#include "peak_picker.h"
#include "fit_2d_peak.h"
#include <mib-read.h>
using namespace std;

static int *peaks = 0;
const static int max_size = 400;
static char index_filename[256];
static char peak_list_filename[256];
static char center_filename[256];
static char mib_filename[256];
static double *fit_data;
static inline double rad2(double &x, double &y){return x*x+y*y;}
static int total_count;
void set_filename(const char file[])
{
    struct stat st;
    if(stat(file, &st) != 0) {
        cerr << "file does not exist";
        exit(0);
    }
    strncpy(mib_filename, file, 255);
    int len = strlen(mib_filename);
    strncpy(index_filename, file, len-3);
    strncat(index_filename, "index", 255);
    strncpy(peak_list_filename, file, len-3);
    strncat(peak_list_filename, "dat", 255);
    strncpy(center_filename, file, len-3);
    strncat(center_filename, "center", 255);
}

void fit_peaks_all_frames(mib *reader, char filename[], int mask_radius, int fit_radius, int frame0, int frame1, int frame_step, Filter2d *filter, int box_size)
{
    fit_data = new double[2500*3];
    total_count = 0;
    set_filename(filename);
    vector<diffraction_data > diff_vec;
    vector<float> center_vec;
    diff_vec.reserve(10000);
    center_vec.reserve(2000);
    vector<int> frame_end;
    if(!peaks) peaks = new int[max_size];
    double p[9], p_error[9];
    double radius2 = mask_radius*mask_radius;
    reader->set_radius(fit_radius);
//    cout <<frame0 << " "<< frame1<<endl;
    for(int current_frame = frame0; current_frame < frame1; current_frame += frame_step){
   //     double angle = step * current_frame;
   //     double ca = cos(angle);
   //     double sa = sin(angle);
//        cout << "current frame "<<current_frame << endl;
        double cx, cy;
        reader->read(current_frame);
        reader->find_center();
        reader->get_center(cx,cy);
        center_vec.push_back(cx);
        center_vec.push_back(cy);
        double C = reader->bg_scale(mask_radius);
        reader->bg_subtract(C);
        reader->apply(*filter);
 //       cout << current_frame <<" "<<cx <<" "<< cy << endl;
        int count = 0;
        get_peak_list(reader->filtered_data(), 514, 514, 514, box_size,
                      peaks, max_size,  count);
        int len = count;
     //   cout << current_frame <<" "<< count<<endl;
        for(int i = 0; i < len; i+=2){
            double fx, fy;
            fx = peaks[i]   - cx;
            fy = peaks[i+1] - cy;
            if(rad2(fx, fy) < radius2) continue;
            int n = reader->copy_data_set(peaks+i, fit_data);
            fit1(p, p_error, fit_data, n);
            double yy = peaks[i+1] + p[5] - cy;
            double xx = peaks[i] + p[4] - cx;
            double r2 = xx*xx + yy *yy;
    //        cout << "fit finished"<<endl;
            diff_vec.push_back(diffraction_data(current_frame, xx, yy, p_error[4], p_error[5], p[3], p_error[3], p[6], p[7], r2));
            //cout<< current_frame <<" "<< xx <<" "<< yy <<endl;
            total_count++;
        }
        frame_end.push_back(total_count);
    }
//    cout << "frame end "<<frame_end.size()<<endl;
    cout << "num data "<<diff_vec.size()<<endl;

    ofstream fo_index(index_filename);
    fo_index.write(reinterpret_cast<char *> (&(frame_end[0])), sizeof(int)*frame_end.size());
    ofstream fo(peak_list_filename);
    fo.write(reinterpret_cast<char *>(&(diff_vec[0])), sizeof (diff_vec[0])*diff_vec.size());
    ofstream foc(center_filename);
    foc.write(reinterpret_cast<char *> (&(center_vec[0])), sizeof(center_vec[0])* center_vec.size());
    delete [] fit_data;
}
