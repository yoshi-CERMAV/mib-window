//
//  MIB_view.h
//  
//
//  Created by Yoshiharu Nishiyama on 20/11/2023.
//

#ifndef MIB_view_h
#define MIB_view_h
#include <FL/Fl_Widget.H>

#define MAX_NUMPEAK 400
class MIB_view:public Fl_Widget{
public:
    MIB_view(int x, int y, int w, int h, const char *label = 0L):Fl_Widget(x, y, w, h, label){
        data = new uint8_t [w*h];
        set_dummy();
        radius = 30;
    }
    ~MIB_view(){if(data) delete []data;}
    unsigned char* image_ptr(){return data;}
    void set_peaks(double *xy, int len){
        if(len > MAX_NUMPEAK) len = MAX_NUMPEAK;
        peak_len = len;
        int *ptr = peaks;
        double *ptr1 = xy;
        for(int i = 0; i < len; i++, ptr+=2, ptr1+=2){
            ptr[0] = ptr1[0]+x();
            ptr[1] = ptr1[1]+y();
        }
    }
    
    void draw_cross(int *ptr){
        fl_line(ptr[0]-10+x(), ptr[1]+y(), ptr[0]+10+x(), ptr[1]+y());
        fl_line(ptr[0]+x(), ptr[1]-10+y(), ptr[0]+x(), ptr[1]+10+y());
    }
    void set_dummy()
    {
        int *ptr = peaks;
        for(int i = 0; i < 10; i++, ptr+=2){
            ptr[0] = i * 30;
            ptr[1] = i * 30;
        }
        peak_len = 10;
    }
    int* peaks_ptr(){return peaks;}
    void set_peak_len(int i){peak_len = i;}
    void draw(){
        fl_draw_image(data, x(), y(), w(), h(), 1);
        fl_begin_line();
        fl_color(FL_GREEN);
        if(show_peakpos){
            int *ptr = peaks;
            for(int i = 0; i < peak_len/2; i++, ptr+=2){
                draw_cross(ptr);
            }
        }
        if(show_box){
            cout <<"show box"<<endl;
            fl_rect(10, 10, box_size, box_size);
            //       fl_arc(0., 0., 20.0, 0.0, 360.0);
        }
        if(show_center_circle){
            fl_arc(icx, icy, radius, 0.0, 360.0);
            cout << "show center"<<endl;
        }
        fl_end_line();

    }
    void set_show_peakpos(){show_peakpos = true;}
    void unset_show_peakpos(){show_peakpos = false;}
    void set_show_box(){show_box = true;}
    void unset_show_box(){show_box = false;}
    void set_center_circle(){show_center_circle = true;}
    void unset_center_circle(){show_center_circle = false;}
    void set_box(int size){box_size = size;}
    void set_radius(double rad){radius = rad;}
    void set_center(double cx, double cy, double rad)
    {
        icx = cx+x();
        icy = cy+y();
        radius = rad;
        cout << "set center "<<icx <<" "<<icy<<" "<<radius<<endl;
    }
    int get_peak_len(){return peak_len;}
protected:
    uint8_t *data;
    int peaks[MAX_NUMPEAK];
    int peak_len = 0;
    int box_size = 30;
    bool show_peakpos=false;
    bool show_box=false;
    bool show_center_circle=false;
    double icx, icy, radius;
};
#endif /* MIB_view_h */
