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
    Linear(ma::Mesh* m, double mesh_sz)
    {
        sz=mesh_sz;
      mesh = m;
      average = ma::getAverageEdgeLength(m);
      ma::getBoundingBox(m,lower,upper);
    }
    virtual double getValue(ma::Entity* v)
    {
    //  ma::Vector p = ma::getPosition(mesh,v);
    //  double x = p[0], y = p[1], z = p[2];
    //  double Pe = 10;
    //  double uxx = -Pe*Pe * exp(Pe*x) * (1-exp(Pe*y)) / ((1-exp(Pe))*(1-exp(Pe))) ;
    //  double uyy = -Pe*Pe * exp(Pe*y) * (1-exp(Pe*x)) / ((1-exp(Pe))*(1-exp(Pe))) ;
    //  double uxy = Pe*Pe * exp(Pe*y) * exp(Pe*x) / ((1-exp(Pe))*(1-exp(Pe))) ;

    //  ma::Matrix r(uxx, uxy, 0.,
    //               uxy, uyy, 0.,
    //                    0., 0., 1.);
    //  
    //  mth::Matrix<double,3,3> A, Q, Lambda;
    //  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) A(i,j) = r[i][j];

    //  mth::eigenQR(A, Lambda, Q, 100);
    //    
    //  std::vector<double> H(3);
    //  for (int i=0; i<3; ++i) H[i] = std::max(0.01, 0.1 / abs(Lambda(i,i))); 

    //  return std::max(0.01, std::abs(-Pe*Pe*exp(Pe*x)/(1-exp(Pe))));//
//      return *std::max_element(&H[0], &H[0]+3);

       return sz;
    }
  private:
    double sz;
    ma::Mesh* mesh;
    double average;
    ma::Vector lower;
    ma::Vector upper;
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
  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->shouldRefineLayer = true;
  ma::adapt(in);
  m->verify();
  apf::writeVtkFiles("after",m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

