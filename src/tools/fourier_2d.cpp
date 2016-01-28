/* Simple tool to compute a Discrete Fourier Transform of a 2D set of data
 * and the time average of the square moduli of the coefficients.
 * ------------------------------------------------------------------------
 *  ISSUES & LIMITATIONS:
 *  >>> data MUST be formatted as XYZ file: X (or Y or Z) should be the data and the other columns contain the coordinates of the regular grid underlying the data
 *  >>> planar, equally spaced grid only
 *  >>> number of points in the grid = number of atoms of the XYZ input data
 */
#include "clparser.hpp"
#include "ioatoms.hpp"
#include "fftw3.h"
#include <cmath>

using namespace std;
using namespace toolbox;

void banner()
{
   cerr<< " USAGE: fourier_2d [options] < INPUT [ > OUTPUT ]                                  \n"
            << " reads a 2D data series from stdin and computes its Fourier transform by FFT.      \n"
            << " Input should be just a XYZ-formatted file with:\n"
            << "    X column: data to be Fourier-transformer (by default) \n"
            << "    Y/Z columns: coordinates of the (equally spaced) grid point corresponding to the value in X column. \n"
            << " OPTIONS: \n"
            << "    -nframes [nf]   number of frames in the trajectory \n"
            << "    -data [x/y/z]   column of data to apply FFT to {def: x} \n"
            << "    -dx [dx]        sets the grid spacing. {def: 1.0}           \n"
            << endl;
}

int main(int argc, char **argv)
{
   CLParser clp(argc, argv);
   
   unsigned long nframes;
   string datacol;
   double dx=1.0; 
   bool fhelp;
   bool fok=clp.getoption(nframes, "nframes") &&
            clp.getoption(dx, "dx", 1.0) &&
            clp.getoption(fhelp, "h", false) &&
            clp.getoption(datacol, "data", string("x"));   
   
   if (!fok || fhelp) { banner(); return -1; }

   // data of type "frame"
   AtomFrame data_frame;
    
   // read in trajectory
   for ( unsigned long i=0; i<nframes; ++i ) {
        unsigned long npoints=0; // # grid points. This should not change, but needs to be here
        ReadXYZFrame(cin, data_frame, npoints);
        
        // FFT array and allocation
        size_t ndata=sqrt(npoints); size_t ndatah=(ndata/2)+1;
        double *data_in=fftw_alloc_real(ndata*ndata); 
        fftw_complex *data_fft=fftw_alloc_complex(ndata*ndatah);
        
        // fill in FFT arrays properly
        for ( int iy=0; iy<ndata; ++iy ) {
            for ( int ix=0; ix<ndata; ++ix ) {
                //printf("DEBUG: %i %i %i \n", ix, iy, iy*ndata+ix);
                if ( datacol=="y" ) data_in[iy*ndata+ix] = data_frame.ats[iy*ndata+ix].y;
                else if ( datacol=="z" ) data_in[iy*ndata+ix] = data_frame.ats[iy*ndata+ix].z;
                else data_in[iy*ndata+ix]=data_frame.ats[iy*ndata+ix].x;
            }
        }

        // FFT plan & execute
        fftw_plan plan_forward=fftw_plan_dft_r2c_2d(ndata, ndata, data_in, data_fft, FFTW_ESTIMATE);
        fftw_execute(plan_forward);

        // CHECK: print output
        printf(" Output Fourier coefficients:\n");
        printf("\n");
        for (int i=0; i<ndata; ++i) {
            for (int j=0; j<ndatah; ++j) {
                printf(" %4d %4d    %12f  %12f\n", i,j,data_fft[i*ndatah+j][0], data_fft[i*ndatah+j][1]);
            }
        }
   } // FRAMES
   
}
