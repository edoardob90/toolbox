/* Simple tool to compute a Discrete Fourier Transform of a 2D set of data
 * and the time average of the square moduli of the coefficients.
 * ------------------------------------------------------------------------
 *  ISSUES & LIMITATIONS:
 *  >>> data MUST be formatted as XYZ file: X (or Y or Z) should be the data and the other columns contain the coordinates of the regular grid underlying the data
 *  >>> planar, equally spaced grid only
 *  >>> number of points in the grid = number of atoms of the XYZ input data
 */
#include <cmath>
#include "clparser.hpp"
#include "ioatoms.hpp"
#ifdef _OMP
#include <omp.h>
#endif
#include "fftw3.h"

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

   // global variables
   AtomFrame data_frame;
   unsigned long iframe, npoints;
   size_t ndata, ndatah;
   valarray<double> fft_average;
#ifdef _OMP
   int nthreads=omp_get_num_threads();
   cerr<<"Multi-threaded FFTW in use. FFTW will be executed using "<<nthreads<<" threads."<<endl;
   int omp_error=fftw_init_threads();
   if (omp_error==0) cerr<<"**WARNING**: some error occurred while initializing multi-threaded FFTW!"<<endl;
#endif
    
   // read in trajectory
   for ( iframe=0; iframe < nframes; ++iframe ) {
        // read in current frame
        npoints=0;
        ReadXYZFrame(cin, data_frame, npoints);
        ndata=sqrt(npoints); ndatah=(ndata/2)+1;
        
        // at the beginning, resize FFT average vector to accomodate data
        if (iframe==0) { fft_average.resize(ndata*ndatah); fft_average=0.0;} 

        // FFT array and allocation
        double *data_in = fftw_alloc_real(ndata*ndata); 
        fftw_complex *data_fft = fftw_alloc_complex(ndata*ndatah);
        
        // fill in FFT arrays properly
        for ( int i=0; i<ndata*ndata; ++i ) {
            if ( datacol=="y" ) data_in[i] = data_frame.ats[i].y;
            else if ( datacol=="z" ) data_in[i] = data_frame.ats[i].z;
            else data_in[i]=data_frame.ats[i].x;
        }

        // FFT plan & execute
#ifdef _OMP
        fftw_plan_with_nthreads(nthreads);
        fftw_plan plan_forward=fftw_plan_dft_r2c_2d(ndata, ndata, data_in, data_fft, FFTW_ESTIMATE);
#else
        fftw_plan plan_forward=fftw_plan_dft_r2c_2d(ndata, ndata, data_in, data_fft, FFTW_ESTIMATE);
#endif
        
        fftw_execute(plan_forward);

        // CHECK: print output
        //printf(" Output Fourier coefficients:\n");
        //printf("\n");
        //for (int i=0; i<ndata; ++i) {
        //    for (int j=0; j<ndatah; ++j) {
        //        printf(" %4d %4d    %12f  %12f\n", i,j,data_fft[i*ndatah+j][0], data_fft[i*ndatah+j][1]);
        //    }
        //}

        // compute square moduli and accumulate time average
        for ( int i=0; i<ndata*ndatah; ++i ) {
            fft_average[i] += data_fft[i][0]*data_fft[i][0]+data_fft[i][1]*data_fft[i][1]; //)*creal(data_fft[i])+cimag(data_fft[i])*cimag(data_fft[i]);
        }

        // FFT destroy
        fftw_free(data_in) ; fftw_free(data_fft) ; fftw_destroy_plan(plan_forward);
   } // FRAMES

   // compute the time average
   fft_average/=nframes;

   // print output
   //cout<<"# FIELDS"<<endl<<"# k_x  k_y  fft_average[k_y,k_x]"<<endl<<"#"<<endl;
   cout.precision(8); cout.width(15); cout.setf(ios::scientific);
   // we are printing HALF the data in each k-direction, beacause we know that the interesting behaviour is at small k-vectors
   for (int i=0; i<ndatah; ++i) {
       for (int j=0; j<ndatah; ++j) {
           cout<<(2.*constant::pi)*j/(dx*ndata)<<"   "<<(2.*constant::pi)*i/(dx*ndata)<<"   "<<fft_average[i*ndatah+j]<<endl;
       }
       //cout<<endl;
   }
}
