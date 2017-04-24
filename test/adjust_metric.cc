#include <apf.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <cassert>

int main(int argc, char** argv)
{
  assert(argc==3);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(modelFile, meshFile);

  /// Iterate over vertices and update the metric size
  apf::MeshEntity* v;
  apf::MeshIterator* it = m->begin(0);
  while ((v = m->iterate(it))) {
      apf::Vector3 p;
      m->getPoint(v, 0, p);
      for (int i = 0; i < 3; ++i)
          p[i] = p[i] * 1e-6;
      m->setPoint(v, 0, p);
  }
  m->end(it);

  m->verify();
  m->writeNative("new_mesh.smb");
  apf::writeVtkFiles("out",m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

