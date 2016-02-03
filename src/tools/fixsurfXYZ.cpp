/* This tool MUST BE USED in combination with a PLUMED Willard-Chandler surface trajectory
 * fixsurfXYZ check that at each frame, for each of the grid points, there is JUST one
 * surface point.
 * WARNING: this tool may perform something that is NOT PHYSICAL at all! Its purpose is only
 * to provide a "fixed" XYZ in which the surface can be represented as a real-valued function
 * It is based upon the "toolbox" infrastructure
 */
#include <cmath>
#include "ioatoms.hpp"
#include "clparser.hpp"

using namespace std;
using namespace toolbox;

void banner()
{
   cerr<< " USAGE: fixsurfXYZ [options] < INPUT [ > OUTPUT ]                                  \n"
            << " reads a PLUMED Willard-Chandler surface XYZ ... \n"
            << " and check the surface can be represented as a function of (X,Y).      \n"
            << " Input should be just a XYZ-formatted file with:\n"
            << "    X column: surface height (by default) \n"
            << "    Y/Z columns: coordinates of the (equally spaced) grid point corresponding to the value in X column. \n"
            << " OPTIONS: \n"
            << "    -data [x/y/z]  column that contains data {def: X; Y-Z are grid points} \n"
            << endl;
}

int main(int argc, char **argv)
{
    CLparser clp(argc, argv);

    string datacol;
    bool fhelp;
    bool fok=clp.getoption(fhelp,"h","false") &&
            //clp.getoption(nframes,"N") &&
            clp.getoption(datacol,"data", string("x"));

    if (!fok || fhelp) { banner(); return -1; }

    vector<AtomFrame> trajectory; // the surface trajectory
    vector<unsigned int> natoms; // atom per each frame 
    
    // read all the trajectory at once 
    ReadXYZ(cin, trajectory, natoms);

    // loop over frames
    unsigned nframes=trajectory.size();
    for (int i=0; i < nframes; ++i) {
        
        int igrid=0, iatom=0;
        double gx=pow(10,10), gy=pow(10,10), znew=pow(10,10), zold;
        vector<double> grid; grid.resize(nframes*natoms[i]); // vector that contains height at each grid point, at each frame
        
        // loop on the grid
        while (igrid < natoms[i]) {

            // get the values from current grid point
            double ax, ay, az; // az is ALWAYS the height of the surface. BE CAREFUL!
            if (datacol=="y") { ax=trajectory[i].ats[igrid].z; ay=trajectory[i].ats[igrid].x; az=trajectory[i].ats[igrid].y; }
            else if (datacol=="z") { ax=trajectory[i].ats[igrid].x; ay=trajectory[i].ats[igrid].y;az=trajectory[i].ats[igrid].z; }
            else { ax=trajectory[i].ats[igrid].y; ay=trajectory[i].ats[igrid].z; az=trajectory[i].ats[igrid].x;}

            // check if the current grid point is different from the previous one (i.e. previous frame)
            if (abs(gx-ax) > pow(10,-5) || abs(gy-ay)>pow(10,-5)) {
                // here we found that the grid point is different, so...
                // we set the height of the grid to znew
                grid[i*natoms[i]+igrid]=znew;
                // incr grid index
                ++igrid;
                // we set zold: if we are on the first frame, there's no previous frame
                if (i==0) zold=grid[i*natoms[i]+igrid];
                else zold=grid[(i-1)*natoms[i-1]+igrid];
                // lastly, we set znew and update the grid extents
                znew=az; gx=ax; gy=ay;
            
            } else { // we found a duplicate point
                if (abs(znew-zold)>abs(az-zold)) znew=az;
            
            }

            grid.clear(); grid.resize(0);
        } // grid loop
    
    } // frame loop


} // main
