#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <cassert>
#include <algorithm>

#include "mthQR.h"

class Linear : public ma::IsotropicFunction
{
  public:
    Linear(ma::Mesh* m, double mesh_sz) : sz(mesh_sz), mesh(m) {}
    virtual double getValue(ma::Entity* v) { return sz; }
  private:
    double sz;
    ma::Mesh* mesh;
};

int main(int argc, char** argv)
{
  assert(argc==4);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  const double mesh_sz = atof(argv[3]);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  ma::Mesh* m = apf::loadMdsMesh(modelFile,meshFile);
  m->verify();
  Linear sf(m, mesh_sz);
  ma::Input* in = ma::configure(m, &sf);
  in->maximumIterations = 10;
  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->shouldRefineLayer = true;
  ma::adapt(in);

//  ma::runUniformRefinement(m,6);

  m->verify();
  apf::writeVtkFiles("after",m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

