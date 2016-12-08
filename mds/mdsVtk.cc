#include "apf.h"
#include "apfMDS.h"
#include "apfMesh2.h"
#include "gmi.h" /* this is for gmi_getline... */

//#include <cstdio>
//#include <cstring>
//#include <cassert>
//#include <cstdlib>

namespace apf {

Mesh2* loadMdsFromVtk(gmi_model* g, const char* filename)
{
    Mesh2* m = makeEmptyMdsMesh(g, 3, false);
    readVtkFiles(filename, m);
    m->verify();
    gmi_write_dmg(m->getModel(), "model_made.dmg");
    return m;
}

}
