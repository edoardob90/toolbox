#include <math.h>
#include <fstream>
#include <string>
#include "clparser.hpp"
#include "ioatoms.hpp"
#include "ioparser.hpp"

using namespace toolbox;

int main (int argc, char **argv)
{
    CLParser clp(argc, argv);

    AtomFrame af;
    int nfr=0;

    while ( ReadXYZFrame(std::cin, af)){
        ++nfr;
        for (std::vector<AtomData>::const_iterator iat=af.ats.begin(); iat!=af.ats.end(); ++iat) {
            std::cout<<"Prop: "<<iat->props[0]<<std::endl;
        }
    }

    return 0;
}
