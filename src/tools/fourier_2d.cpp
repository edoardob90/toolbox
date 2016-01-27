#include "clparser.hpp"
#include "fftw3.h"

using namespace std;
using namespace toolbox;

void banner()
{
   cerr<< " USAGE: fourier_2d [options] < INPUT [ > OUTPUT ]                                  \n"
            << " reads a 2D data series from stdin and computes its Fourier transform by FFT.      \n"
            << " Input should be just the data, nothing else.            \n"           
            << " Data is interpreted as a sequence of points corresponding to an uniform grid with equal spacing DX. \n"
            << " -dx [dx]    sets the grid spacing. {def: 1.0}           \n"
            << " -pad [npad] appends npad zeroes before doing the FT. increases the resolution. \n"
            << endl;
}

int main(int argc, char **argv)
{
   CLParser clp(argc, argv);
   
   double dx=1.0; 
   unsigned long npad=0;
   string wnd;
   bool fhelp;
   bool fok=clp.getoption(dx, "dx", 1.0) &&
            clp.getoption(fhelp, "h", false) &&
            clp.getoption(npad, "pad", (unsigned long) 0);   
   
   if (!fok || fhelp) { banner(); return -1; }

   // fill vector with data
   vector<double> data(0); double y;
   while (cin.good()) { cin >> y; data.push_back(y); }

   // prepare FFT
   unsigned long ndata=data.size(), nfft=ndata+npad;
   valarray<double> vvt(nfft), vvw(nfft);

   vvt=0.0; for (unsigned long i=0; i<nfft; ++i) vvt[i]=data[i];
}
