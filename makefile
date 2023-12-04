mib-window: mib_window.cc fit_all_frames.cc
	g++ -o mib_window mib_window.cc fit_all_frames.cc /usr/local/lib/libfltk.a -lpthread -framework Cocoa -I/Users/yoshi/include  -std=c++11 -framework Accelerate /Users/yoshi/lib/Savitzky_Golay_2d.o /Users/yoshi/lib/Savitzky_Golay_3d.o -lfftw3 peak_picker.cc fit_2d_peak.c -lmpfit -L/Users/yoshi/lib
