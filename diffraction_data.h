//
//  diffraction_data.h
//  
//
//  Created by Yoshiharu Nishiyama on 25/11/2023.
//

#ifndef diffraction_data_h
#define diffraction_data_h
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sys/stat.h>
#include <vector>
#define DIFF_DATA_ARRAY_SIZE 9
using namespace std;
class diffraction_data{
public:
    diffraction_data()
    {
//        data = new float[DIFF_DATA_ARRAY_SIZE];
    }
    diffraction_data(int in_frame, float in_x, float in_y, float in_x_err, float in_y_err, float in_Intensity, float in_Intensity_err, float in_width_x, float in_width_y, float r2)
    {
 //       data = new float[DIFF_DATA_ARRAY_SIZE];
        frame = in_frame;
        data[0] = in_x;
        data[1] = in_y;
        data[2] = in_x_err;
        data[3] = in_y_err;
        data[4] = in_Intensity;
        data[5] = in_Intensity_err;
        data[6] = in_width_x;
        data[7] = in_width_y;
        data[8] = r2;
    }
    void write(ofstream &fo){
        fo.write(reinterpret_cast <char *>(&frame), sizeof(int) );
        fo.write(reinterpret_cast <char *>(data), sizeof(float)*DIFF_DATA_ARRAY_SIZE);
    }
    void dump()
    {
        cout << frame <<"\t"<< x() <<"\t"<< y()<<endl;
    }
//    int read(ifstream &fi){
//        if(!fi.read(reinterpret_cast<char *>(frame), sizeof(int))) return 0;
//        return fi.read(reinterpret_cast <char *>(data), sizeof(float)* 6);
//    }
    int frame;
    float data[DIFF_DATA_ARRAY_SIZE];
    float &x(){return (data[0]);}
    float &y(){return (data[1]);}
    float x_err(){return data[2];}
    float y_err(){return data[3];}
    float Intensity(){return data[4];}
    float Intensity_err(){return data[5];}
    float width_x(){return data[6];}
    float width_y(){return data[7];}
    float r2(){return data[8];}
} ;

class diffraction_data_set
{
    typedef struct{
        diffraction_data *d0;
        diffraction_data *d1;
        int frame;
    } frame_data;
    void read_diffraction_data()
    {
        struct stat st;
        if(stat(data_filename,&st)!=0) {
            cerr<<"file does not exist "<<data_filename<<endl;
        }
        diffraction_data dat;
        number_of_data = st.st_size / sizeof(dat);
        cout << "number of data "<<number_of_data<<endl;
        if(stat(index_filename,&st)!=0) {
            cerr<<"file does not exist "<<index_filename<<endl;
        }
        
        number_of_frames = st.st_size/sizeof(int);
        cout << "number of frames "<<number_of_frames << endl;
        int *frame_end = new int[number_of_frames];
        ifstream fi_index(index_filename);
        fi_index.read(reinterpret_cast<char *>(frame_end), sizeof(int)*number_of_frames);
        frames = new frame_data[number_of_frames];
        
        data = new diffraction_data [number_of_data];
        ifstream fi(data_filename);
        cout << "size of (diff_dat) "<< sizeof(data[0])<<endl;
        fi.read(reinterpret_cast<char *>(data), sizeof(data[0])*number_of_data);
        frames[0].d0 = data;
        for(int i = 0; i < number_of_frames; i++) {
            frames[i].d1 = data+frame_end[i];
            if(i) frames[i].d0 = frames[i-1].d1;
            frames[i].frame = i;
        }
        delete [] frame_end;
        for(int i = 0; i < 10; i++) data[i].dump();
        cout <<"------------"<<endl;
        data[number_of_data-1].dump();
    }
    
public:
    diffraction_data_set(char *name){
        set_filename(name);
        read_diffraction_data();
        coord.reserve(100000);
    }
    ~diffraction_data_set()
    {
        if(data) delete[]data;
    }
    inline void rotate_z(float &x, float &y, float *out)
    {
        out[0] = cr * x + sr * y;
        out[1] = -sr * x + cr * y;
    }

    inline void rotate_x(float *x){
        x[2] = x[1] * srot;
        x[1] *= crot;
    }

    inline void geom_conv(float *ptr, diffraction_data *diff)
    {
        rotate_z(diff->x(), diff->y(), ptr);
        rotate_x(ptr);
    }
    void head()
    {
        for(int i = 0; i < 10; i++){
            cout << coord[i*3]<<" "<<coord[i*3+1]<<" "<<coord[i*3+2]<<endl;
        }
    }
    void copy_coordinates(float max, int frame_begin, int frame_end, int step, float ax0, float frame_rate){
        cr = cos(ax0);
        sr = sin(ax0);
        float temp[3];
        float radius2 = max*max;
        cout << "copying "<<endl;
        for(int i = frame_begin; i < frame_end; i += step){
    //        cout << i << " "<< coord.size()<<" "<<frames[i].d0->x()<<endl;
             float rot = i * frame_rate;
             srot = sin(rot);
             crot = cos(rot);
             for(diffraction_data *ptr = frames[i].d0; ptr != frames[i].d1; ptr++){
                 if(ptr->r2() < radius2){
                     geom_conv(temp, ptr);
                     coord.push_back(temp[0]);
                     coord.push_back(temp[1]);
                     coord.push_back(temp[2]);
                }
            }
        }
        cout <<"size "<< coord.size()<<endl;
    }
    float *coordinates(){return &(coord[0]);}
    int len(){return coord.size()/3;}
    void rotate_coordinates(float offset_z, float offset_y, float offset_x);
protected:
    char index_filename[256];
    char data_filename[256];
    char center_filename[256];
    int number_of_data;
    int number_of_frames;
    diffraction_data *data;
    frame_data *frames;
    vector<float> coord;
    
    float cr, sr, srot, crot;
    void set_filename(const char file[])
    {
        strncpy(index_filename, file, 255);
        strncat(index_filename, ".index", 255);
        strncpy(data_filename, file, 255);
        strncat(data_filename, ".dat", 255);
        strncpy(center_filename, file, 255);
        strncat(center_filename, ".center", 255);
    }

};

#endif /* diffraction_data_h */
