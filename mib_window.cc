#include <FL/Fl.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Image.H>
#include <FL/Fl_Round_Button.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Progress.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <Savitzky_Golay_2d.h>
#include <Savitzky_Golay_3d.h>
#include <mib-read.h>
#include <MIB_view.h>
#include <YN_input.h>
#include "diffraction_data.h"
#include "fit_all_frames.h"
void fit1(double p[], double p_error[], double *data, int data_size);
void fit(double p[], double p_error[], double *data, int data_size);

static int max_num_data = 200*1000;
static int total_count = 0;
diffraction_data *diffdata;
#define MAINWIN_W       700                             // main window w()
#define MAINWIN_H       700                             // main window h()
#define BROWSER_X       10                              // browser x()
#define BROWSER_Y       10                              // browser y()
#define BROWSER_W       100                             // browser w()
#define BROWSER_H       MAINWIN_H-35                    // browser h()
#define DATA_AREA_X      (BROWSER_W + 20)                // test area x()
#define DATA_AREA_Y      10                              // test area y()
#define DATA_AREA_W      (MAINWIN_W - BROWSER_W - 30)    // test area w()
#define DATA_AREA_H      BROWSER_H                       // test area h()
int get_peak_list(double *data, int w, int h, int line_len, int box,
                  int *peaklist, int max_size,  int &count);
int box_size= 30;
void set_mul_threshold(double a) ;
void set_threshold(double a);
char mib_filename[255];

static int data_mode = 0;
static double *fit_data;
int *peaklist;
Filter2d *filter;
Filter3d *filter3d;
static int current_frame;
static int current_peak;
static int number_of_frames;
static int frame_step= 1;
static double C;
double rotation_per_frame = 0.165*M_PI/180;
double st = sin(M_PI/180.*154);
double ct = cos(M_PI/180.*154);
static void rotate(double &x, double &y){
    double temp = ct * x + st *y;
    y = -st * x + ct * y;
    x = temp;
}
double step = 0.0165*M_PI/180;

double cx, cy, radius = 30;
class MIB_Window: public Fl_Double_Window{
public:
    MIB_Window(int w, int h, const char *l=0L) :
    Fl_Double_Window(w, h, l){}
};

YN_Int_Input *input, *sg_size_input, *sg_order_input;
YN_Int_Input *box_size_input;
YN_Int_Input *frame0, *frame1;
int iframe0, iframe1;
YN_Int_Input *frame_step_input;
YN_Input *max_input, *sigma_x_input, *sigma_input;
YN_Input *scale_input, *center_radius_input;
Fl_Output *cb_info=(Fl_Output *)0;
YN_Input *axe_input, *rot_step_input;
Fl_Progress *progress;
class rotation_series{
public:
    
protected:
    double axe;
    double c_axe, s_axe;
};

//YN_Input *check_peakpos, *check_box, *check_center;

mib *reader;
MIB_view *view;

void set_progress(int i)
{
    progress->value(i);
    progress->redraw();
}

void copy_mode()
{
    switch(data_mode){
        case 0:
            reader->copy(view->image_ptr());
            break;
        case 1:
            reader->copy_d(view->image_ptr());
            break;
        case 2:
            reader->copy_o(view->image_ptr());
            break;
    }
}

void sg3d_param_changed(Fl_Widget *w, void *data)
{
    int n = sg_size_input->val();
    filter3d->init_polynom(n, n, n, 0, sg_order_input->val());
    filter3d->init_filter();
    reader->apply(*filter);
    copy_mode();
    int count = 0;
    get_peak_list(reader->filtered_data(), 514, 514, 514, 30,
                  view->peaks_ptr(), 200,  count);
    view->set_peak_len(count);
    cout <<"num peaks"<< count/2 <<endl;
    view->redraw();
}

#ifdef AAA
void get_peak_list(Fl_Widget *w, void *data)
{
    current_peak = 0;
    current_frame = 0;
    while(reader->read(current_frame)){
        reader->bg_subtract(C);
        reader->apply(*filter);
        get_peak_list(reader->filtered_data(), 514, 514, 514, box_size,
                      diff_peak);
        current_frame++;
    }
}
#endif

inline double rad2(double &x, double &y){return x*x+y*y;}

void fit_peaks_all_frames(Fl_Widget *w, void *data)
{
    fit_peaks_all_frames(reader, mib_filename, radius, 6,  iframe0, iframe1, frame_step, filter,box_size);
}

void fit_peaks(Fl_Widget *w, void *data)
{
    reader->set_radius(6);
    int *peaks = view->peaks_ptr();
    int len = view->get_peak_len();
    double p[9], p_error[9];
    ofstream fo("temp");
    double radius2 = radius*radius;
    for(int i = 0; i < len; i+=2){
        double fx, fy;
        fx  = peaks[i]-cx;
        fy = peaks[i+1]-cy;
        if(rad2(fx, fy) < radius2) continue;
        int n = reader->copy_data_set(peaks+i, fit_data);
        //        for(int j = 0; j < n; j++){
        //            fo << j <<" "<< fit_data[j] << " "<< fit_data[j + n]<<" "<< fit_data[j+2*n]<<endl;
        //        }
        fit1(p, p_error, fit_data, n);
        cout << current_frame<<" ";
        cout <<peaks[i] + p[4] - cx<<" "<< p_error[4] <<" "<<peaks[i+1] + p[5] - cy<<" ";
        for(int i = 3; i < 4; i ++){
        cout << p[i] <<" "<<p_error[i]<<" ";
        }
    
        cout << endl;
        //      if (i==2) exit(0);
    }
}

void change_bg()
{
    C = reader->bg_scale(radius);
    reader->bg_subtract(C);
    reader->apply(*filter);
    int count = 0;
    get_peak_list(reader->filtered_data(), 514, 514, 514, box_size,
                  view->peaks_ptr(), 200,  count);
    view->set_peak_len(count);
    copy_mode();
}

void iframe_changed(Fl_Widget *w, void *data)
{
    iframe0 = frame0->val();
    iframe1 = frame1->val();
}

void frame_step_changed(Fl_Widget *w, void *data)
{
    frame_step = frame_step_input->val();
}
void file_num_changed(Fl_Widget *w, void *data)
{
    current_frame  = input->val();
    //    cout << "current frame = "<< current_frame <<endl;
    reader->read(current_frame);
    //    cout << "read "<< current_frame <<endl;
    reader->set_max(max_input->get());
    reader->find_center();
    //    cout << "center found"<<endl;
    reader->get_center(cx,cy);
    view->set_center(cx,cy,radius);
    //    cout << "change bg"<<endl;
    change_bg();
    view->redraw();
}

void sg_param_changed(Fl_Widget *w, void *data)
{
    int n = sg_size_input->val();
    filter->init_polynom(n, n, 0, sg_order_input->val());
    filter->init_filter();
    reader->apply(*filter);
    copy_mode();
    int count = 0;
    get_peak_list(reader->filtered_data(), 514, 514, 514, 30,
                  view->peaks_ptr(), 200,  count);
    view->set_peak_len(count);
    cout <<"num peaks"<< count/2 <<endl;
    view->redraw();
}

void max_changed(Fl_Widget *w, void *data)
{
    reader->set_max(max_input->get());
    copy_mode();
    view->redraw();
}
void scale_changed(Fl_Widget *w, void *data)
{
    reader->set_scale(scale_input->get());
    file_num_changed(NULL, NULL);
}

void box_size_changed(Fl_Widget *w, void *data)
{
    box_size = box_size_input->val();
    int count = 0;
    get_peak_list(reader->filtered_data(), 514, 514, 514, box_size,
                  view->peaks_ptr(), 200,  count);
    view->set_box(box_size);
    view->set_peak_len(count);
    //   cout <<"num peaks"<< count/2 <<endl;
    char msg[256];
    snprintf(msg,256, "%d peaks found", count/2 );
    cb_info->value(msg);
    cb_info->redraw();
    view->redraw();
}

void center_radius_changed(Fl_Widget *w, void *data)
{
    radius = center_radius_input->get();
    double C = reader->bg_scale(radius);
    reader->bg_subtract(C);
    reader->apply(*filter);
    view->set_radius(radius);
    view->redraw();
}

void peak_pick_conditions_changed(Fl_Widget *w, void *data)
{
    set_mul_threshold(sigma_x_input->get());
    set_threshold(sigma_input->get());
    box_size_changed(NULL, NULL);
}


static void button_cb(Fl_Button *b, void *) {
    char msg[256];
    //  sprintf(msg, "Label: '%s'\nValue: %d", b->label(),b->value());
    if(! strncmp(b->label(), "subtracted", 6)) data_mode = 1;
    if(! strncmp(b->label(), "raw", 3)) data_mode = 0;
    if(! strncmp(b->label(), "filtered", 6)) data_mode = 2;
    file_num_changed(NULL, NULL);
}

static void check_button_cb(Fl_Check_Button *b, void *) {
    char msg[256];
    //  sprintf(msg, "Label: '%s'\nValue: %d", b->label(),b->value());
    if(! strncmp(b->label(), "peak", 3)) {\
        cout << "peak "<<(int)(b->value())<<endl;
        if( b->value())view->set_show_peakpos(); else view->unset_show_peakpos();
    }
    if(! strncmp(b->label(), "box", 3)){
        cout << "box "<<(int)(b->value())<<endl;
        if( b->value())view->set_show_box(); else view->unset_show_box();
    }
    if(! strncmp(b->label(), "center", 3)) {
        if( b->value())view->set_center_circle(); else view->unset_center_circle();
        cout << "center "<<(int)(b->value())<<endl;
    }
    view->redraw();
    //  file_num_changed(NULL, NULL);
}
using namespace std;


int main(int argc, char *argv[])
{
 //   diff_data = new diffraction_data[max_num_data];
//    set_filename(argv[1]);
    strncpy(mib_filename, argv[1], 255);
    reader = new mib(argv[1], .8);
    filter = new Filter2d (2, 2, 0, 2, 514, 514);
    peaklist = new int[200];
    fit_data = new double[2500*3];
    //    diff_peaklist = new diffraction_peak[200*1000];
    current_peak = 0;
    number_of_frames = reader->number_of_frames();
    Fl_Window window(900, 600, "diffraction frame");
    window.begin();
    view = new MIB_view(10, 10, 514, 514, "diffraction");
    int vpos = 10;
    int hpos = 624;
    input = new YN_Int_Input(hpos, vpos, 50, 25, "frame num", 1);
    vpos +=30;
    input->callback(file_num_changed, NULL);
    max_input = new YN_Input(hpos, vpos, 50, 25, "max", 255);
    vpos +=30;
    max_input->tooltip("maximum pixel value");
    max_input->callback(max_changed, NULL);
    scale_input =new YN_Input(hpos, vpos, 50, 25, "scale", 0.8);
    vpos +=30;
    scale_input->tooltip("scaling factor for the pixels atvthe border of chips");
    scale_input->callback(scale_changed, NULL);
    center_radius_input = new YN_Input(hpos, vpos, 50, 25, "radius", radius);
    vpos +=50;
    center_radius_input->tooltip("radius for fitting the central diffusion against 1/r2");
    center_radius_input->callback(center_radius_changed, NULL);
    int hpos1 = 600;
    Fl_Group* o = new Fl_Group(hpos1, vpos, 90, 90, "chose display");
    o->box(FL_THIN_UP_FRAME);
    { Fl_Round_Button* o = new Fl_Round_Button(hpos1, vpos, 90, 30, "raw");
        vpos +=30;
        o->tooltip("showing detector image.");
        o->type(102);
        o->set();
        o->down_box(FL_ROUND_DOWN_BOX);
        o->callback((Fl_Callback*) button_cb);
    } // Fl_Round_Button* o
    { Fl_Round_Button* o = new Fl_Round_Button(hpos1, vpos, 90, 30, "subtrated"); vpos += 30;
        o->tooltip("showing image after direct beam subtraction");
        o->type(102);
        o->down_box(FL_ROUND_DOWN_BOX);
        o->callback((Fl_Callback*) button_cb);
    }
    { Fl_Round_Button* o = new Fl_Round_Button(hpos1, vpos, 90, 30, "filtered");
        vpos+=30;
        o->tooltip("showing image after Savitzky-Golay filter application");
        o->type(102);
        o->down_box(FL_ROUND_DOWN_BOX);
        o->callback((Fl_Callback*) button_cb);
    }
    { cb_info = new Fl_Output(534, vpos, 160, 50, "");
        vpos +=70;
        cb_info->type(12);
        cb_info->textsize(12);
        cb_info->align(Fl_Align(133));
    } // Fl_Output* cb_info
    
    o->end();
    // Fl_Group* o
    {Fl_Group* o = new Fl_Group(550, vpos, 150, 60, "savitzky golay parameters");
        o->box(FL_THIN_UP_FRAME);
        sg_size_input = new YN_Int_Input(624, vpos, 50, 25, "width", 2);
        vpos+=30;
        sg_size_input->tooltip("the width of the windows on both sides");
        sg_size_input->callback((Fl_Callback *) sg_param_changed);
        sg_order_input = new YN_Int_Input(624, vpos, 50, 25, "order", 2);
        vpos+=70;
        sg_order_input->callback((Fl_Callback *) sg_param_changed);
        o->end();
    }
    {
        Fl_Group* o = new Fl_Group(550, vpos, 150, 80, "peak search") ;
        o->box(FL_THIN_UP_FRAME);
        o->tooltip("only peaks with intensity above x times the local standard deviation are chosen");
        box_size_input = new YN_Int_Input(624, vpos, 50, 25, "box size", box_size);
        vpos+=30;
        box_size_input->callback((Fl_Callback *)box_size_changed);
        box_size_input->tooltip("size of the tile in which we try to find a peak. Should be smaller than the minimum distance between two peaks");
        sigma_x_input = new YN_Input(624, vpos, 50, 25, "sigma x", 15.);
        vpos += 30;
        sigma_x_input->tooltip("only peaks with intensity above x times the local standard deviation are chosen");
        sigma_x_input->callback((Fl_Callback *)peak_pick_conditions_changed );
        sigma_input = new YN_Input (624, vpos, 50, 25, "cutoff", 5.);
        vpos+=30;
        sigma_input->tooltip("only peaks with intensity above this value are chosen");
        sigma_input->callback((Fl_Callback *)peak_pick_conditions_changed );
        //        input = new YN_Int_Input(624, 307, 50, 25, "order", 2);
        o->end();
    }
    // cout<< "file num changed"<<endl;
    file_num_changed(NULL, NULL);
    //  cout<< "file num changed"<<endl;
    
    //    hpos += 100;
    //    check_peakpos = Check_Button(hpos, vpos, 5, 5));
    //    vpos += 30;
    hpos = 750;
    vpos = 10;
    Fl_Check_Button * show_cross = new Fl_Check_Button(hpos, vpos, 50, 15, "peak");
    show_cross->down_box(FL_DOWN_BOX);
    show_cross->callback((Fl_Callback*)  check_button_cb);
    vpos += 25;
    Fl_Check_Button * show_center = new Fl_Check_Button(hpos, vpos, 50, 15, "center");
    show_center->callback((Fl_Callback*) check_button_cb);
    vpos += 25;
    Fl_Check_Button * show_box= new Fl_Check_Button(hpos, vpos, 50, 15, "box");
    show_box->callback((Fl_Callback*) check_button_cb);
    vpos += 25;
    Fl_Button *b4 = new Fl_Button(hpos, vpos, 70, 15, "fit peaks");
    b4->callback((Fl_Callback *)fit_peaks);
    vpos += 25;
    Fl_Button *b5 = new Fl_Button(hpos, vpos, 70, 15, "fit-all frames");
    b5->callback((Fl_Callback *)fit_peaks_all_frames);
    vpos += 25;
    progress = new Fl_Progress(hpos, vpos, 100, 30);
    progress->minimum(0);
    progress->maximum(number_of_frames);
    vpos +=25;
    frame0 = new YN_Int_Input(hpos, vpos, 100, 30, "1st frame", 0);
    vpos +=25;
    iframe0 = 0;
    iframe1 = number_of_frames;
    frame1 = new YN_Int_Input(hpos, vpos, 100, 30, "last frame", number_of_frames);
    frame0->callback((Fl_Callback*) iframe_changed);
    frame1->callback((Fl_Callback*) iframe_changed);
    vpos +=25;
    frame_step_input = new YN_Int_Input(hpos, vpos, 100, 30, "frame step", 1);
    frame_step_input->callback((Fl_Callback*) frame_step_changed);
    
    window.end();
    window.show(argc, argv);
    return Fl::run();
}
