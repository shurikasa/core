#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <cassert>
#include <iostream>
#include <cmath>
#include <cstdlib>

#include "mthQR.h"

class AnIso : public ma::AnisotropicFunction
{
  public:
    AnIso(ma::Mesh* m)
    {
      mesh = m;
      average = ma::getAverageEdgeLength(m);
      ma::getBoundingBox(m,lower,upper);
    }
    virtual void getValue(ma::Entity* v, ma::Matrix& R, ma::Vector& H)
    {
      ma::Vector p = ma::getPosition(mesh,v);
      double x = (p[0] - lower[0])/(upper[0] - lower[0]);
      double sizeFactor = x + 2.0;
      ma::Vector h(0.1, 0.1, 1.);
      double pi = std::acos(-1.);
      double theta = pi/4.;
//      ma::Matrix r(std::cos(theta), -std::sin(theta), 0.0,
//		   std::sin(theta), std::cos(theta), 0.0,
//		   0.0, 0.0, 1.0);
      ma::Matrix r(6., 9., 0., 3., 5., 0., 0., 0., 1.);
      
      std::cout << r << std::endl << std::endl; 
      mth::Matrix<double,3,3> A, Q, Lambda;
      for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
              A(i,j) = r[i][j];
      mth::eigenQR(A, Q, Lambda, 100);

      std::cout << A << std::endl << std::endl; 
      std::cout << Lambda << std::endl << std::endl; 
      std::cout << Q << std::endl << std::endl; 

      exit(-1);
      H = h;
      R = r;
    }
  private:
    ma::Mesh* mesh;
    double average;
    ma::Vector lower;
    ma::Vector upper;
};

int main(int argc, char** argv)
{
  assert(argc==3);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  ma::Mesh* m = apf::loadMdsMesh(modelFile,meshFile);
  m->verify();
  apf::writeVtkFiles("aniso_before",m);
  AnIso sf(m);
  ma::Input* in = ma::configure(m, &sf);
  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->shouldRefineLayer = true;
  ma::adapt(in);
  m->verify();
  apf::writeVtkFiles("aniso_after",m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

