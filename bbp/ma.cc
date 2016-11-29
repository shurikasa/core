#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <cassert>


double getDistance(ma::Mesh* mesh, ma::Entity** verts)
{
  return (ma::getPosition(mesh, verts[1]) - ma::getPosition(mesh, verts[0])).getLength();
}

/// Needs adjustments for parallel, i.e. calculation of average edge on the part boundary
class Iso : public ma::IsotropicFunction
{
  public:
    Iso(ma::Mesh* m)
    {
      mesh = m;
    }
    virtual double getValue(ma::Entity* v)
    {
        ma::Entity* e;
        double esize = 0.;

        int nent = mesh->countUpward(v);
        for (int i = 0; i < nent; ++i) {
            e = mesh->getUpward(v, i);
            apf::Downward verts;
            mesh->getDownward(e, 0, verts);
            esize += getDistance(mesh, verts);

        }
        esize /= nent;
        esize /= 1.5;
        return esize;
    }
  private:
    ma::Mesh* mesh;
};


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
      ma::Vector h(average, average/sizeFactor, average/sizeFactor);
      ma::Matrix r(1.0, 0.0, 0.0,
		   0.0, 1.0, 0.0,
		   0.0, 0.0, 1.0);
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
  Iso sf(m);
//  AnIso sf(m);
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

